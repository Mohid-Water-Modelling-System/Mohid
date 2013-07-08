!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Assimilation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to compute data assimilation
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

Module ModuleAssimilation

    use ModuleGlobalData
    use ModuleTime
    use ModuleStopWatch
    use ModuleTimeSerie,        only: StartTimeSerie,                                   &
                                      GetTimeSerieLocation, CorrectsCellsTimeSerie,     &
                                      GetNumberOfTimeSeries, TryIgnoreTimeSerie,        &
                                      WriteTimeSerie, KillTimeSerie         

    use ModuleHDF5,             only: ConstructHDF5, GetHDF5FileAccess, HDF5SetLimits,      &
                                      HDF5WriteData, HDF5FlushMemory, KillHDF5
    use ModuleEnterData,        only: ReadFileName, ConstructEnterData, GetData,            &
                                      ExtractBlockFromBuffer, Block_Unlock, GetOutPutTime,  &
                                      ExtractBlockFromBlock, KillEnterData
    use ModuleGridData,         only: GetGridData, UngetGridData        
    use ModuleHorizontalGrid,   only: WriteHorizontalGrid, UnGetHorizontalGrid,             &
                                      GetGridCellArea, GetXYCellZ
    use ModuleHorizontalMap,    only: GetWaterPoints2D, UngetHorizontalMap, GetWaterFaces2D 
    use ModuleGeometry,         only: GetGeometrySize, GetGeometryDistances, UngetGeometry, &
                                      GetGeometryKFloor
    use ModuleMap,              only: GetWaterPoints3D, GetOpenPoints3D, GetWaterFaces3D, UngetMap
    use ModuleFillMatrix,       only: ConstructFillMatrix, GetDefaultValue, KillFillMatrix, &
                                      ModifyFillMatrix
    use ModuleFunctions,        only: ConstructPropertyID, Density, Sigma, CHUNK_J

    implicit none

    private 

    !Subroutines & Functions---------------------------------------------------

    !Constructor
    public  ::  StartAssimilation                                          
    private ::      ReadAssimilationFilesName
    private ::      ConstructAssimilationList
    private ::          ConstructProperty  
    private ::              ConstructAssimilationField
    private ::              ConstructPropertyCoefficients
    private ::                  FindMinimumCoef
    private ::                  TranslateTypeZUV
    private ::      ConstructAltimetricAssimilation
    private ::      Open_HDF5_OutPut_File
    private ::      ConstructTimeSerie

    !Selector
    public  ::  GetAssimilationSize
    public  ::  GetAssimilationField
    private ::      SearchProperty
    private ::          AssimilationFromFile
    public  ::  GetAssimilationCoef
    public  ::  GetAssimilationAltimetry
    public  ::  GetAssimilationAltimetryDT
    public  ::  GetAltimetryDecayTime
    public  ::  GetAltimSigmaDensAnalyzed

    public  ::  UngetAssimilation

    !Modifier
    public  ::  ModifyAssimilation
    private ::      ReadLockExternalVar
    private ::      AltimetricAssimilation
    private ::          LoweringLifting   
    private ::              ChangeTemperatureSalinity
    private ::                  Spline_Alt
    private ::                  Splint
    private ::      ReadUnLockExternalVar
    private ::  OutPutResultsHDF
    private ::  OutPut_TimeSeries



    !Destructor
    public  ::  KillAssimilation                                              
    private ::      DeallocateVariables  
    private ::      KillAltimetricAssimilation
    private ::      WriteFinalHDFOutPut
    

    !Management
    private ::      Ready
    private ::      LocateObjAssimilation


    !Interfaces----------------------------------------------------------------
    private :: UngetAssimilation2D
    private :: UngetAssimilation3D
    private :: UngetAssimilation3D8
    interface  UngetAssimilation
        module procedure UngetAssimilation2D
        module procedure UngetAssimilation3D
        module procedure UngetAssimilation3D8
    end interface UngetAssimilation

    !Parameter-----------------------------------------------------------------

    !Property Dimensions 
    integer, parameter                          :: Dim_2D               = 2
    integer, parameter                          :: Dim_3D               = 3
                                                                        
    character(LEN = StringLength), parameter    :: block_begin          = '<beginproperty>'
    character(LEN = StringLength), parameter    :: block_end            = '<endproperty>'
    
    character(LEN = StringLength), parameter    :: begin_field          = '<<begin_field>>'
    character(LEN = StringLength), parameter    :: end_field            = '<<end_field>>'

    character(LEN = StringLength), parameter    :: begin_coef           = '<<begin_coef>>'
    character(LEN = StringLength), parameter    :: end_coef             = '<<end_coef>>'

    !Types---------------------------------------------------------------------

    private :: T_Field
    type       T_Field
        real, dimension(:,:  ), pointer         :: R2D
        real, dimension(:,:,:), pointer         :: R3D
        real                                    :: Minimum
        real                                    :: DefaultValue
        integer                                 :: TypeZUV              = null_int
    end type   T_Field

    private :: T_OutPut
    type       T_OutPut
         type (T_Time), pointer, dimension(:)   :: OutTime
         integer                                :: Next
         logical                                :: ON
         integer                                :: Number
    end type T_OutPut


    private :: T_Property
    type       T_Property
        type (T_PropertyID)                     :: ID
        type (T_PropertyID)                     :: CoefID
        integer                                 :: Dim                  = null_int
        type (T_Field)                          :: Field
        type (T_Field)                          :: CoefField
        real                                    :: ColdRelaxPeriod, ColdOrder
        type (T_Time)                           :: LastActualization
        logical                                 :: TimeSerie            = .false.
        logical                                 :: OutputHDF            = .false.
        type (T_Time)                           :: LastTimeSerieOutput
        type (T_Property), pointer              :: Next, Prev
    end type T_Property

    private :: T_Files
    type       T_Files
         character(len=PathLength)              :: ConstructData
         character(len=PathLength)              :: Results
    end type T_Files

    ! guillaume nogueira : aqui define-se a estrutura de assimilacao de altimetria

    type       T_External                                                   ! J. Nogueira Assim
        type(T_Time)                            :: Now
        real,    pointer, dimension(:,:,:)      :: DWZ
        real,    pointer, dimension(:,:,:)      :: SZZ
        real,    pointer, dimension(:,:  )      :: GridCellArea
        real,    pointer, dimension(:,:,:)      :: ZCellCenter
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D
        integer, pointer, dimension(:,:  )      :: KFloor_Z
    end type T_External

    private :: T_Altimetric_Assim
    type       T_Altimetric_Assim                                           ! J. Nogueira Assim
        real                                    :: DtAltimAssimilation
        real                                    :: AltimetricDepth
        real                                    :: AltimDecayTime
        logical                                 :: UseVarianceField         ! Guillaume
        real, dimension(:    ), pointer         :: ColumnTemperature
        real, dimension(:    ), pointer         :: ColumnSalinity
        real, dimension(:    ), pointer         :: ColumnDensity
        real, dimension(:    ), pointer         :: DeltaTemperature
        real, dimension(:    ), pointer         :: DeltaSalinity
        real, dimension(:    ), pointer         :: DeltaDensity
        real, dimension(:    ), pointer         :: DeltaPressure
        real, dimension(:    ), pointer         :: AuxDepth
        real, dimension(:    ), pointer         :: D_GradTemp
        real, dimension(:    ), pointer         :: D_GradSalin
        real, dimension(:    ), pointer         :: AuxSpline
        real, dimension(:,:  ), pointer         :: WaterLevelToAssimilate
        real, dimension(:,:  ), pointer         :: VarianceFieldToAssimilate
        real, dimension(:,:  ), pointer         :: Delta_WaterLevel
        real, dimension(:,:  ), pointer         :: Delta_WaterLevelToAssimilate
        real, dimension(:,:  ), pointer         :: Delta_WaterLevelColumn
        real, dimension(:,:  ), pointer         :: PressureAnomaly
        real, dimension(:,:  ), pointer         :: WaterLevelAnalized
        real, dimension(:,:  ), pointer         :: Observation_Error
        real, dimension(:,:  ), pointer         :: Model_Error
        real, dimension(:,:  ), pointer         :: Gain
        real, dimension(:,:,:), pointer         :: SigmaDensAnalyzed
        real, dimension(:,:,:), pointer         :: TemperatureAnalyzed
        real, dimension(:,:,:), pointer         :: SalinityAnalyzed
        type(T_Time)                            :: NextCompute
    end type T_Altimetric_Assim


    private :: T_Assimilation
    type       T_Assimilation
        private
        logical                                 :: AltimetricAssimilation   ! J. Nogueira Assim
        integer                                 :: InstanceID       
        type(T_Time)                            :: LastCall
        type(T_Time)                            :: EndTime      
        type(T_Time)                            :: ActualTime
        type(T_Size3D)                          :: Size
        type(T_Size3D)                          :: WorkSize
        type(T_Files)                           :: Files
        type(T_OutPut)                          :: OutPut
        type(T_External )                       :: ExternalVar              ! J. Nogueira Assim
        type(T_Altimetric_Assim)                :: Altimetric_Assim         ! J. Nogueira Assim
        type(T_Property), pointer               :: FirstAssimilationProp
        type(T_Property), pointer               :: LastAssimilationProp
        integer                                 :: PropertiesNumber   = FillValueInt

        !Instance of ModuleHDF5     
        integer                                 :: ObjHDF5            = 0

        !Instance of Module_EnterData
        integer                                 :: ObjEnterData       = 0
       
        !Instance of ModuleBathymetry
        integer                                 :: ObjGridData        = 0
      
        !Instance of ModuleHorizontalGrid
        integer                                 :: ObjHorizontalGrid  = 0

        !Instance of ModuleGeometry
        integer                                 :: ObjGeometry        = 0

        !Instance of ModuleHorizontalMap
        integer                                 :: ObjHorizontalMap   = 0

        !Instance of ModuleMap
        integer                                 :: ObjMap             = 0

        !Instance of ModuleTime
        integer                                 :: ObjTime            = 0

        !Instance of ModuleTimeSerie
        integer                                 :: ObjTimeSerie       = 0

        ! guillaume : Required by the Cooper-Haines method
        !Instance of ModuleWaterProperties     
        integer                                 :: ObjWaterProperties = 0

        type(T_Assimilation), pointer           :: Next
    
    end type T_Assimilation

    !Global Module Variables
    type (T_Assimilation), pointer              :: FirstObjAssimilation
    type (T_Assimilation), pointer              :: Me



    !--------------------------------------------------------------------------
    
    contains



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine StartAssimilation(AssimilationID,                            &
                                 TimeID,                                    &
                                 GridDataID,                                &
                                 HorizontalGridID,                          &
                                 HorizontalMapID,                           &
                                 MapID,                                     &
                                 GeometryID,                                &
                                 STAT)

        !Arguments--------------------------------------------------------------
        integer                                 :: AssimilationID
        integer                                 :: TimeID
        integer                                 :: GridDataID
        integer                                 :: HorizontalGridID
        integer                                 :: HorizontalMapID
        integer                                 :: MapID
        integer                                 :: GeometryID
        integer, optional, intent(OUT)          :: STAT     

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: ready_         

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_
 
        !Begin-----------------------------------------------------------------
        
        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mAssimilation_)) then
            nullify (FirstObjAssimilation)
            call RegisterModule (mAssimilation_) 
        endif

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            !Associates External Instances
            Me%ObjTime            = AssociateInstance (mTIME_,           TimeID            )
            Me%ObjGridData        = AssociateInstance (mGRIDDATA_,       GridDataID        )
            Me%ObjHorizontalMap   = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID   )
            Me%ObjHorizontalGrid  = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID  )
            Me%ObjGeometry        = AssociateInstance (mGEOMETRY_,       GeometryID        )
            Me%ObjMap             = AssociateInstance (mMAP_,            MapID             ) 


            call GetComputeTimeLimits(Me%ObjTime, BeginTime = Me%ActualTime, &
                                      EndTime = Me%EndTime, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'StartAssimilation - ModuleAssimilation - ERR01'


            call GetGeometrySize(Me%ObjGeometry, Size = Me%Size, WorkSize = Me%WorkSize, &
                                 STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'StartAssimilation - ModuleAssimilation - ERR02'

            !Read the name file of the Assimilation module
            call ReadAssimilationFilesName

            !Construct enter data 
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL) 
            if(STAT_CALL .ne. SUCCESS_)stop 'StartAssimilation - ModuleAssimilation - ERR03'

            ! guillaume : aqui le-se os parametros do assimilation.dat
            ! Constructs the property list 
            call ConstructAssimilationList
                     
            ! guillaume : aqui chama-se a rotina de alocacoes.
            if(Me%AltimetricAssimilation)                                   &
                    call ConstructAltimetricAssimilation

            !Constructs Time Serie
            call ConstructTimeSerie

            !Constructs HDF output time
            call ConstructOutPut

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'StartAssimilation - ModuleAssimilation - ERR04'

            ! By default a output file is always open in the construction phase
            if (Me%OutPut%ON) call Open_HDF5_OutPut_File

            call null_time(Me%ActualTime)

            AssimilationID = Me%InstanceID

            STAT_ = SUCCESS_
        
        else 
            
            stop 'StartAssimilation - ModuleAssimilation - ERR99' 

        end if 

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartAssimilation

    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance
                                                    
        !Local-----------------------------------------------------------------
        type (T_Assimilation), pointer                  :: NewObjAssimilation
        type (T_Assimilation), pointer                  :: PreviousObjAssimilation

        !Begin-----------------------------------------------------------------

        !Allocates new instance
        allocate (NewObjAssimilation)
        nullify  (NewObjAssimilation%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjAssimilation)) then
            FirstObjAssimilation         => NewObjAssimilation
            Me                           => NewObjAssimilation
        else
            PreviousObjAssimilation      => FirstObjAssimilation
            Me                           => FirstObjAssimilation%Next
            do while (associated(Me))
                PreviousObjAssimilation  => Me
                Me                       => Me%Next
            enddo
            Me                           => NewObjAssimilation
            PreviousObjAssimilation%Next => NewObjAssimilation
        endif

        Me%InstanceID = RegisterNewInstance (mASSIMILATION_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------
    !Read the name of the files need to construct and modify
    ! the Assimilation properties 

    subroutine ReadAssimilationFilesName

        !External--------------------------------------------------------------
        character(len = StringLength)   :: Message
        type(T_Time)                    :: TIME_END
        integer                         :: STAT_CALL

        !Begin------------------------------------------------------------------

        TIME_END = Me%EndTime

        !Opens the Assimilation data file 
        ! ASCII file used to construct new properties
        Message  ='Assimilation Data Properties.'
        Message  = trim(Message)

        call ReadFileName('ASSIMILA_DAT', Me%Files%ConstructData,               &
                           Message = Message, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadAssimilationFilesName - ModuleAssimilation - ERR01'

        ! ---> File in HDF format where is written instant fields of Assimilation properties
        Message   ='Instant fields of Assimilation properties in HDF format.'
        Message   = trim(Message)

        call ReadFileName('ASSIMILA_HDF', Me%Files%Results, Message = Message,  &
                           TIME_END = TIME_END, Extension = 'ass', STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadAssimilationFilesName - ModuleAssimilation - ERR02'


    end subroutine ReadAssimilationFilesName

    !--------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !External--------------------------------------------------------------
        real,    pointer, dimension(:,:  )      :: Bathymetry
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D
        integer                                 :: STAT_CALL, HDF5_CREATE

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, ILB, JUB, JLB, KLB, KUB
        character(len=StringLength)             :: FileName

        !----------------------------------------------------------------------

        IUB  = Me%WorkSize%IUB 
        JUB  = Me%WorkSize%JUB 
        KUB  = Me%WorkSize%KUB 

        ILB  = Me%WorkSize%ILB 
        JLB  = Me%WorkSize%JLB 
        KLB  = Me%WorkSize%KLB 


        FileName = trim(adjustl(Me%Files%Results))

        call GetGridData(Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR00'

        !Gets WaterPoints3D
        call GetWaterPoints3D   (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR01'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5      (Me%ObjHDF5,                                             &
                                 trim(Me%Files%Results)//"5",                            &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR02'


        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR03'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR04'

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = Bathymetry,                                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR05'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",                 &
                              Array3D = WaterPoints3D,                                   &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR06'


        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR07'

        !Ungets the WaterPoints
        call UnGetMap        (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR08'
        
        call UngetGridData(Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAssimilation - ERR09'


    end subroutine Open_HDF5_OutPut_File


    !--------------------------------------------------------------------------

    
    subroutine ConstructTimeSerie

        !External--------------------------------------------------------------
        character(len=StringLength)                         :: TimeSerieLocationFile
        integer, dimension(:,:,:), pointer                  :: WaterPoints3D
        integer                                             :: STAT_CALL, iflag
       
        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber        
        integer                                             :: nProperties
        type(T_Property), pointer                           :: PropertyX
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Begin-----------------------------------------------------------------

        call GetWaterPoints3D(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAssimilation - ERR10'


        !First checks out how many properties will have time series
        PropertyX   => Me%FirstAssimilationProp
        nProperties =  0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) nProperties = nProperties + 1
            PropertyX=>PropertyX%Next
        enddo

        if (nProperties > 0) then

            !Allocates PropertyList
            allocate(PropertyList(nProperties))

            !Fills up PropertyList
            PropertyX   => Me%FirstAssimilationProp
            nProperties =  0
            do while (associated(PropertyX))
                if (PropertyX%TimeSerie) then
                    nProperties = nProperties + 1
                    PropertyList(nProperties) = trim(adjustl(PropertyX%ID%name))
                endif
                PropertyX=>PropertyX%Next
            enddo

            call GetData(TimeSerieLocationFile,                                             &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromFile,                                           &
                         keyword      = 'TIME_SERIE_LOCATION',                              &
                         ClientModule = 'ModuleWaterProperties',                            &
                         Default      = Me%Files%ConstructData,                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAssimilation - ERR20'


            !Constructs TimeSerie
            call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                                trim(TimeSerieLocationFile),                                &
                                PropertyList, "sra",                                        &
                                WaterPoints3D = WaterPoints3D,                              &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAssimilation - ERR30'

            !Deallocates PropertyList
            deallocate(PropertyList)

            !Corrects if necessary the cell of the time serie based in the time serie coordinates
            call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAssimilation - ERR40'

            do dn = 1, TimeSerieNumber

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          CoordX   = CoordX,                                &
                                          CoordY   = CoordY,                                & 
                                          CoordON  = CoordON,                               &
                                          STAT     = STAT_CALL)

                if (CoordON) then
                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= OUT_OF_BOUNDS_ERR_) then
                        stop 'ConstructTimeSerie - ModuleAssimilation - ERR50'
                    endif                            

                    if (STAT_CALL == OUT_OF_BOUNDS_ERR_ .or. Id < 0 .or. Jd < 0) then
                                    
                        call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAssimilation - ERR60'

                        if (IgnoreOK) then
                            cycle
                        else
                            stop 'ConstructTimeSerie - ModuleAssimilation - ERR70'
                        endif

                    endif

                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAssimilation - ERR80'

                endif

            enddo

        endif

        call UnGetMap(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAssimilation - ERR70'

    end subroutine ConstructTimeSerie

    !--------------------------------------------------------------------------

    subroutine ConstructAssimilationList

        !Local-----------------------------------------------------------------
        integer                                 :: iflag
        type (T_Property),    pointer                       :: NewProperty
        integer                                             :: ClientNumber
        integer                                             :: STAT_CALL
        logical                                             :: BlockFound

        !----------------------------------------------------------------------
        ! Initialize the Assimilation properties number   

        Me%PropertiesNumber = 0

        ! Initialize the Assimilation properties list   
        nullify (Me%FirstAssimilationProp)
        nullify (Me%LastAssimilationProp)

        ! J. Nogueira Assim

        call GetData(Me%AltimetricAssimilation,                                &
                     Me%ObjEnterData, iflag,                                                    &
                     SearchType = FromFile,                                                     &
                     keyword    = 'ALTIMETRIC_ASSIMILATION',                                    &
                     Default    = .false.,                                                      &
                     ClientModule = 'ModuleAssimilation',                                       &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructAssimilationList - ModuleAssimilation - ERR01'

        ! J. Nogueira Assim

        if (Me%AltimetricAssimilation)then

            call GetData(Me%Altimetric_Assim%DtAltimAssimilation,                               &
                         Me%ObjEnterData, iflag,                                                &
                         SearchType = FromFile,                                                 &
                         keyword    = 'ALTIMETRIC_DT',                                  &
                         Default    = 864000.,                                                    &
                         ClientModule = 'ModuleAssimilation',                                   &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructAssimilationList - ModuleAssimilation - ERR02'

            call GetData(Me%Altimetric_Assim%AltimetricDepth,                               &
                         Me%ObjEnterData, iflag,                                                &
                         SearchType = FromFile,                                                 &
                         keyword    = 'ALTIMETRIC_DEPTH',                                  &
                         Default    = 50.,                                                    &
                         ClientModule = 'ModuleAssimilation',                                   &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructAssimilationList - ModuleAssimilation - ERR03'

            !Valor estimado com L = 50 km e c_onda_rossby = 0.05 m/s e fazendo tau=L/c
            call GetData(Me%Altimetric_Assim%AltimDecayTime,                               &
                         Me%ObjEnterData, iflag,                                                &
                         SearchType = FromFile,                                                 &
                         keyword    = 'ALTIMETRIC_DECAYTIME',                                   &
                         Default    = 1.e6,                                                     &
                         ClientModule = 'ModuleAssimilation',                                   &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructAssimilationList - ModuleAssimilation - ERR04'

            ! Guillaume
            call GetData(Me%Altimetric_Assim%UseVarianceField,                               &
                         Me%ObjEnterData, iflag,                                                &
                         SearchType = FromFile,                                                 &
                         keyword    = 'USE_VARIANCE_FIELD',                                   &
                         Default    = .false.,                                                     &
                         ClientModule = 'ModuleAssimilation',                                   &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructAssimilationList - ModuleAssimilation - ERR04'

            if(Me%Altimetric_Assim%AltimDecayTime .le. 0.) then
                write(*,*) 'Error: negative altimetry decay time.'
                stop 'ConstructAssimilationList - ModuleAssimilation - ERR05'
            end if

         end if

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,      &
                                        block_begin, block_end, BlockFound, &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  
                    
                    ! Construct a New Property 
                    Call ConstructProperty(NewProperty, ClientNumber)

                    ! Add new Property to the Assimilation List 
                    Call Add_Property     (NewProperty)
                else
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) &
                        stop 'ConstructAssimilationList - ModuleAssimilation - ERR01'
                        
                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_) &
                    stop 'ConstructAssimilationList - ModuleAssimilation - ERR02'
            end if cd1
        end do do1

        !----------------------------------------------------------------------

    end subroutine ConstructAssimilationList

    !--------------------------------------------------------------------------
    ! nogueira e guillaume : alocacao de variaveis para Cooper-Haines
    subroutine ConstructAltimetricAssimilation

        !----------------------------------------------------------------------
        integer                                     :: IUB, JUB, KUB, ILB, JLB, KLB

         !Begin Shorten variable names

        IUB  = Me%WorkSize%IUB 
        JUB  = Me%WorkSize%JUB 
        KUB  = Me%WorkSize%KUB 

        ILB  = Me%WorkSize%ILB 
        JLB  = Me%WorkSize%JLB 
        KLB  = Me%WorkSize%KLB 

       !Alocation of internal variables--------------------------------------------

        allocate (  Me%Altimetric_Assim%ColumnTemperature(KLB:KUB),                                               &
                    Me%Altimetric_Assim%ColumnSalinity(KLB:KUB),                                                   &
                    Me%Altimetric_Assim%ColumnDensity(KLB:KUB),                                                   &
                    Me%Altimetric_Assim%DeltaTemperature(KLB:KUB),                                                &
                    Me%Altimetric_Assim%DeltaSalinity(KLB:KUB),                                                   &
                    Me%Altimetric_Assim%DeltaDensity(KLB:KUB),                                                    &
                    Me%Altimetric_Assim%DeltaPressure(KLB:KUB))

        allocate (  Me%Altimetric_Assim%Delta_WaterLevel(ILB:IUB, JLB:JUB),                                       &
                    Me%Altimetric_Assim%Delta_WaterLevelToAssimilate(ILB:IUB, JLB:JUB),                           &
                    Me%Altimetric_Assim%PressureAnomaly(ILB:IUB, JLB:JUB))

        allocate (  Me%Altimetric_Assim%Observation_Error(ILB:IUB, JLB:JUB),                                      &
                    Me%Altimetric_Assim%Model_Error(ILB:IUB, JLB:JUB),                                            &
                    Me%Altimetric_Assim%Gain(ILB:IUB, JLB:JUB))

        allocate (  Me%Altimetric_Assim%AuxDepth(KLB:KUB))
        allocate (  Me%Altimetric_Assim%D_GradTemp(KLB:KUB))
        allocate (  Me%Altimetric_Assim%D_GradSalin(KLB:KUB))
        allocate (  Me%Altimetric_Assim%AuxSpline(KLB:KUB)) 
        
        allocate (  Me%Altimetric_Assim%SigmaDensAnalyzed(ILB:IUB, JLB:JUB, KLB:KUB))


        !----------------------------------------------------------------------

    end subroutine ConstructAltimetricAssimilation

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine ConstructProperty(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer           :: NewProperty
        integer                             :: ClientNumber

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------

        !Allocates new property
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProperty - ModuleAssimilation - ERR01'

        nullify(NewProperty%Field%R2D)
        nullify(NewProperty%Field%R3D)

        nullify(NewProperty%CoefField%R2D)
        nullify(NewProperty%CoefField%R3D)
        
        call null_time(NewProperty%LastActualization)


        call ConstructPropertyID            (NewProperty%ID, Me%ObjEnterData, FromBlock)


        call ConstructAssimilationField     (NewProperty, ClientNumber)


        call ConstructPropertyCoefficients  (NewProperty, ClientNumber)

        !----------------------------------------------------------------------

    end subroutine ConstructProperty


    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct the property values       
    ! in the domain and in the boundaries            

    subroutine ConstructAssimilationField(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property),          pointer      :: NewProperty
        integer                                 :: ClientNumber

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, i, j, k
        integer, dimension(:,:  ), pointer      :: WaterPoints2D, WaterFaces2D_U, WaterFaces2D_V, PointsToFill2D
        integer, dimension(:,:,:), pointer      :: WaterPoints3D, WaterFaces3D_U, WaterFaces3D_V, PointsToFill3D 
        real,    dimension(:,:  ), pointer      :: Matrix2D
        real,    dimension(:,:,:), pointer      :: Matrix3D


        !Local-----------------------------------------------------------------
        integer                                 :: iflag
        integer                                 :: SizeILB, SizeIUB, SizeJLB
        integer                                 :: SizeJUB, SizeKLB, SizeKUB
        character(len=StringLength)             :: Char_TypeZUV
        logical                                 :: BlockFound

        !----------------------------------------------------------------------
 
        SizeILB = Me%Size%ILB
        SizeIUB = Me%Size%IUB
        SizeJLB = Me%Size%JLB
        SizeJUB = Me%Size%JUB
        SizeKLB = Me%Size%KLB
        SizeKUB = Me%Size%KUB

        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR10'

        call GetWaterFaces2D(Me%ObjHorizontalMap,                                               &
                               WaterFaces2DU = WaterFaces2D_U,      &
                               WaterFaces2DV = WaterFaces2D_V,      &
                               STAT            = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR20'

        call GetWaterFaces3D(Me%ObjMap,                                               &
                               WaterFacesU3D = WaterFaces3D_U,                      &
                               WaterFacesV3D = WaterFaces3D_V,                      &
                               STAT            = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR30'

        call GetWaterPoints3D(Me%ObjMap,WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR40'


        !Initial value of NewProperty%LastActualization = =-9999999
        call null_Time(NewProperty%LastActualization)


        !By default the property dimension is R3D (scalar) 
        call GetData(NewProperty%Dim, Me%ObjEnterData, iflag,                           &
                     keyword        = 'DIMENSION',                                      &  
                     default        = Dim_3D,                                           &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAssimilation',                             &
                     STAT           = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR50'


        call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                       &
                                   begin_field, end_field, BlockFound,                  &
                                   STAT = STAT_CALL)

        if(STAT_CALL .EQ. SUCCESS_)then

            if (BlockFound) then

                call GetData(Char_TypeZUV, Me%ObjEnterData, iflag,                      &
                             keyword        = 'TYPE_ZUV',                               &  
                             SearchType     = FromBlockInBlock,                         &
                             ClientModule   = 'ModuleAssimilation',                     &
                             default        = "Z",                                      &
                             STAT           = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR60'

                if (NewProperty%Dim == Dim_2D) then 

                    allocate(NewProperty%Field%R2D (SizeILB:SizeIUB, SizeJLB:SizeJUB), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR70'

                    NewProperty%Field%R2D(:,:) = FillValueReal
            
                    NewProperty%Field%TypeZUV  = TranslateTypeZUV(Char_TypeZUV)

                    if (NewProperty%Field%TypeZUV == TypeZ_) then

                        PointsToFill2D => WaterPoints2D

                        if (GetPropertyIDNumber(NewProperty%ID%Name) == BarotropicVelocityU_ .or. &
                            GetPropertyIDNumber(NewProperty%ID%Name) == BarotropicVelocityV_) then

                            allocate (Matrix2D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                            Matrix2D(:,:) = FillValueReal

                        else

                            Matrix2D => NewProperty%Field%R2D

                        endif

                    else if (NewProperty%Field%TypeZUV == TypeU_) then

                        PointsToFill2D => WaterFaces2D_U
                        Matrix2D       => NewProperty%Field%R2D

                    else if (NewProperty%Field%TypeZUV == TypeV_) then

                        PointsToFill2D => WaterFaces2D_V
                        Matrix2D       => NewProperty%Field%R2D

                    else

                        stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR80'

                    endif



                    call ConstructFillMatrix  (PropertyID           = NewProperty%ID,           &
                                               EnterDataID          = Me%ObjEnterData,          &
                                               TimeID               = Me%ObjTime,               &
                                               HorizontalGridID     = Me%ObjHorizontalGrid,     &
                                               ExtractType          = FromBlockInBlock,         &
                                               PointsToFill2D       = PointsToFill2D,           &
                                               Matrix2D             = Matrix2D,                 &
                                               TypeZUV              = NewProperty%Field%TypeZUV,&
                                               STAT                 = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR90'

                    call GetDefaultValue(NewProperty%ID%ObjFillMatrix, NewProperty%CoefField%DefaultValue, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR100'

                    if (NewProperty%Field%TypeZUV == TypeZ_) then
                        if (GetPropertyIDNumber(NewProperty%ID%Name) == BarotropicVelocityU_) then

                            do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                    
                                NewProperty%Field%R2D(i,j) = 0.

                                if (PointsToFill2D(i,j-1) == OpenPoint .and.          &
                                    PointsToFill2D(i,j  ) == OpenPoint) then
                                    NewProperty%Field%R2D(i,j) = (Matrix2D(i,j-1)+Matrix2D(i,j))/2.
                                endif

                            enddo
                            enddo

                            deallocate(Matrix2D)

                        endif

                        if (GetPropertyIDNumber(NewProperty%ID%Name) == BarotropicVelocityV_) then

                            do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                    
                                NewProperty%Field%R2D(i,j) = 0.

                                if (PointsToFill2D(i-1,j) == OpenPoint .and.          &
                                    PointsToFill2D(i  ,j) == OpenPoint) then
                                    NewProperty%Field%R2D(i,j) = (Matrix2D(i-1,j)+Matrix2D(i,j))/2.
                                endif

                            enddo
                            enddo

                            deallocate(Matrix2D)

                        endif
                    endif

                    nullify   (PointsToFill2D)
                    nullify   (Matrix2D)


                    if(.not. NewProperty%ID%SolutionFromFile)then

                        call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR110'
            
                    end if


                else if (NewProperty%Dim == Dim_3D) then

                    allocate(NewProperty%Field%R3D   (SizeILB:SizeIUB, SizeJLB:SizeJUB, SizeKLB:SizeKUB), &
                             STAT = STAT_CALL)            
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR120'

                    NewProperty%Field%R3D(:,:,:) = FillValueReal

                    NewProperty%Field%TypeZUV   = TranslateTypeZUV(Char_TypeZUV)

                    if (NewProperty%Field%TypeZUV == TypeZ_) then

                        PointsToFill3D => WaterPoints3D

                        if (GetPropertyIDNumber(NewProperty%ID%Name) == VelocityU_ .or. &
                            GetPropertyIDNumber(NewProperty%ID%Name) == VelocityV_) then

                            allocate (Matrix3D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
                            Matrix3D(:,:,:) = FillValueReal

                        else

                            Matrix3D => NewProperty%Field%R3D

                        endif

                    else if (NewProperty%Field%TypeZUV == TypeU_) then

                        PointsToFill3D => WaterFaces3D_U
                        Matrix3D       => NewProperty%Field%R3D

                    else if (NewProperty%Field%TypeZUV == TypeV_) then

                        PointsToFill3D => WaterFaces3D_V
                        Matrix3D       => NewProperty%Field%R3D

                    else

                        stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR130'

                    endif


                    call ConstructFillMatrix  (PropertyID           = NewProperty%ID,           &
                                               EnterDataID          = Me%ObjEnterData,          &
                                               TimeID               = Me%ObjTime,               &
                                               HorizontalGridID     = Me%ObjHorizontalGrid,     &
                                               GeometryID           = Me%ObjGeometry,           &
                                               ExtractType          = FromBlockInBlock,         &
                                               PointsToFill3D       = PointsToFill3D,           &
                                               Matrix3D             = Matrix3D,                 &
                                               TypeZUV              = NewProperty%Field%TypeZUV,&
                                               STAT                 = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR140'

                    call GetDefaultValue(NewProperty%ID%ObjFillMatrix, NewProperty%Field%DefaultValue, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR150'

                    if (NewProperty%Field%TypeZUV == TypeZ_) then
                        if (GetPropertyIDNumber(NewProperty%ID%Name) == VelocityU_) then

                            do k = Me%WorkSize%KLB,Me%WorkSize%KUB
                            do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                    
                                NewProperty%Field%R3D(i,j,k) = 0.

                                if (PointsToFill3D(i,j-1,k) == OpenPoint .and.          &
                                    PointsToFill3D(i,j  ,k) == OpenPoint) then
                                    NewProperty%Field%R3D(i,j,k) = (Matrix3D(i,j-1,k)+Matrix3D(i,j,k))/2.
                                endif

                            enddo
                            enddo
                            enddo

                            deallocate(Matrix3D)

                        endif

                        if (GetPropertyIDNumber(NewProperty%ID%Name) == VelocityV_) then

                            do k = Me%WorkSize%KLB,Me%WorkSize%KUB
                            do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                    
                                NewProperty%Field%R3D(i,j,k) = 0.

                                if (PointsToFill3D(i-1,j,k) == OpenPoint .and.          &
                                    PointsToFill3D(i  ,j,k) == OpenPoint) then
                                    NewProperty%Field%R3D(i,j,k) = (Matrix3D(i-1,j,k)+Matrix3D(i,j,k))/2.
                                endif

                            enddo
                            enddo
                            enddo

                            deallocate(Matrix3D)

                        endif
                    endif

                    nullify   (PointsToFill3D)
                    nullify   (Matrix3D)


                    if(.not. NewProperty%ID%SolutionFromFile)then

                        call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR160'
            
                    end if
               
                else

                    stop 'ConstructAssimilationField - ModuleAssimilation - ERR170'

                end if

            else

                stop 'ConstructAssimilationField - ModuleAssimilation - ERR180'

            end if

        else

            stop 'ConstructAssimilationField - ModuleAssimilation - ERR190'

        end if

        !Along this period the relaxation term grows lineary until the value specified
        !<BeginKeyword>
            !Keyword          : COLD_RELAX_PERIOD
            !<BeginDescription>       
               ! 
               !The user specify the period along which wants the relaxation have a linear growth
               !
            !<EndDescription>
            !Type             : Real 
            !Default          : 0.

            !File keyword     : IN_DAD3D 
            !Multiple Options : 0 (.false.), 1 (.true.)
            !Search Type      : From File
        !<EndKeyword>


        !By default property do not have a cold relaxation period
        call GetData(NewProperty%ColdRelaxPeriod,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'COLD_RELAX_PERIOD',                              &
                     Default        = 0.,                                               &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAssimilation',                             &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR200'

        !By default cold coefficient follows a x^4 evolution
        call GetData(NewProperty%ColdOrder,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'COLD_ORDER',                                     &
                     Default        = 4.,                                               &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAssimilation',                             &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR210'

        !By default property don't have time serie
        call GetData(NewProperty%TimeSerie,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TIME_SERIE',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAssimilation',                             &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR220'

        if (NewProperty%TimeSerie) then
            NewProperty%LastTimeSerieOutput = Me%ActualTime
        endif


        !By default property don't have time serie
        call GetData(NewProperty%OutputHDF,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'OUTPUT_HDF',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAssimilation',                             &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR230'

        call UngetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR240'

        call UngetHorizontalMap(Me%ObjHorizontalMap, WaterFaces2D_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR250'

        call UngetHorizontalMap(Me%ObjHorizontalMap, WaterFaces2D_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAssimilationField - ModuleAssimilation - ERR260'

        
        call UnGetMap(Me%ObjMap, WaterFaces3D_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR270'

        call UnGetMap(Me%ObjMap, WaterFaces3D_V, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR280'

        call UnGetMap(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR290'

        !----------------------------------------------------------------------

    end subroutine ConstructAssimilationField




    
    subroutine ConstructPropertyCoefficients(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property),     pointer           :: NewProperty
        integer                                 :: ClientNumber

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer, dimension(:,:  ), pointer      :: WaterPoints2D, PointsToFill2D 
        integer, dimension(:,:,:), pointer      :: WaterPoints3D, PointsToFill3D
        integer                                 :: iflag
        logical                                 :: BlockFound

        !Local-----------------------------------------------------------------
        integer                                 :: SizeILB, SizeIUB, SizeJLB
        integer                                 :: SizeJUB, SizeKLB, SizeKUB
        integer                                 :: di, dj, i, j, k
        character(len=StringLength)             :: Char_TypeZUV

        !----------------------------------------------------------------------
 
        SizeILB = Me%Size%ILB
        SizeIUB = Me%Size%IUB
        SizeJLB = Me%Size%JLB
        SizeJUB = Me%Size%JUB
        SizeKLB = Me%Size%KLB
        SizeKUB = Me%Size%KUB


        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR10'
        
        call GetWaterPoints3D(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR20'


        call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,      &
                                   begin_coef, end_coef, BlockFound,   &
                                   STAT = STAT_CALL)

        if(STAT_CALL .EQ. SUCCESS_)then

cd0:        if (BlockFound) then

                call GetData(Char_TypeZUV, Me%ObjEnterData, iflag,                          &
                             keyword        = 'TYPE_ZUV',                                   &  
                             SearchType     = FromBlockInBlock,                             &
                             ClientModule   = 'ModuleAssimilation',                         &
                             default        = "Z",                                          &
                             STAT           = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR03'

                NewProperty%CoefField%TypeZUV = TranslateTypeZUV(Char_TypeZUV)

                if (NewProperty%Dim == Dim_2D) then 

                    allocate(NewProperty%CoefField%R2D (SizeILB:SizeIUB, SizeJLB:SizeJUB))
                    allocate(PointsToFill2D            (SizeILB:SizeIUB, SizeJLB:SizeJUB))
                    
                    if (NewProperty%CoefField%TypeZUV == TypeZ_) then
                    
                        di=0;dj=0
                        
                    else if (NewProperty%CoefField%TypeZUV == TypeU_) then
                    
                        di=0;dj=1                   
                    
                    else if (NewProperty%CoefField%TypeZUV == TypeV_) then
                    
                        di=1;dj=0
                    
                    endif

                    do j = Me%WorkSize%JLB,Me%WorkSize%JUB + 1
                    do i = Me%WorkSize%ILB,Me%WorkSize%IUB + 1
                        
                        if (WaterPoints2D(i-di,j-dj) == WaterPoint .or.                 &
                            WaterPoints2D(i   ,j   ) == WaterPoint) then

                            PointsToFill2D(i,j) = 1

                        else
                            
                            PointsToFill2D(i,j) = 0

                            
                        endif

                    enddo
                    enddo 
                    

                    NewProperty%CoefField%R2D(:,:) = FillValueReal


                    call ConstructFillMatrix  (PropertyID           = NewProperty%CoefID,           &
                                               EnterDataID          = Me%ObjEnterData,              &
                                               TimeID               = Me%ObjTime,                   &
                                               HorizontalGridID     = Me%ObjHorizontalGrid,         &
                                               ExtractType          = FromBlockInBlock,             &
                                               PointsToFill2D       = WaterPoints2D,                &
                                               Matrix2D             = NewProperty%CoefField%R2D,    &
                                               TypeZUV              = NewProperty%CoefField%TypeZUV,&
                                               STAT                 = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR05'

                    call GetDefaultValue(NewProperty%CoefID%ObjFillMatrix, NewProperty%Field%DefaultValue, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR06'
                    
                    call FindMinimumCoef(NewProperty%CoefField%Minimum,                 &
                                         Field2D       = NewProperty%CoefField%R2D,     &
                                         WaterPoints2D = PointsToFill2D)

                    if(.not. NewProperty%CoefID%SolutionFromFile)then

                        call KillFillMatrix(NewProperty%CoefID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR07'
                    
                    end if
                    
                    deallocate(PointsToFill2D)


                else if (NewProperty%Dim == Dim_3D) then

                    allocate(NewProperty%CoefField%R3D (SizeILB:SizeIUB, SizeJLB:SizeJUB, SizeKLB:SizeKUB))

                    allocate(PointsToFill3D            (SizeILB:SizeIUB, SizeJLB:SizeJUB, SizeKLB:SizeKUB))
                    
                    if (NewProperty%CoefField%TypeZUV == TypeZ_) then
                    
                        di=0;dj=0
                        
                    else if (NewProperty%CoefField%TypeZUV == TypeU_) then
                    
                        di=0;dj=1                   
                    
                    else if (NewProperty%CoefField%TypeZUV == TypeV_) then
                    
                        di=1;dj=0
                    
                    endif

                    do k = Me%WorkSize%KLB,Me%WorkSize%KUB
                    do j = Me%WorkSize%JLB,Me%WorkSize%JUB + 1
                    do i = Me%WorkSize%ILB,Me%WorkSize%IUB + 1
                        
                        if (WaterPoints3D(i-di,j-dj,k) == WaterPoint .or.               &
                            WaterPoints3D(i   ,j   ,k) == WaterPoint) then

                            PointsToFill3D(i,j,k) = 1

                        else
                            
                            PointsToFill3D(i,j,k) = 0

                            
                        endif

                    enddo
                    enddo 
                    enddo

                    NewProperty%CoefField%R3D(:,:,:) = FillValueReal



                    call ConstructFillMatrix  (PropertyID           = NewProperty%CoefID,           &
                                               EnterDataID          = Me%ObjEnterData,              &
                                               TimeID               = Me%ObjTime,                   &
                                               HorizontalGridID     = Me%ObjHorizontalGrid,         &
                                               GeometryID           = Me%ObjGeometry,               &
                                               ExtractType          = FromBlockInBlock,             &
                                               PointsToFill3D       = PointsToFill3D,               &
                                               Matrix3D             = NewProperty%CoefField%R3D,    &
                                               TypeZUV              = NewProperty%CoefField%TypeZUV,&
                                               STAT                 = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR09'

                    call GetDefaultValue(NewProperty%CoefID%ObjFillMatrix, NewProperty%CoefField%DefaultValue, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR10'

                    call FindMinimumCoef(NewProperty%CoefField%Minimum,                 &
                                         Field3D       = NewProperty%CoefField%R3D,     &
                                         WaterPoints3D = PointsToFill3D)

                    if(.not. NewProperty%CoefID%SolutionFromFile)then

                        call KillFillMatrix(NewProperty%CoefID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR11'

                    end if
                    
                    deallocate(PointsToFill3D)
               
                else

                    stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR12'

                end if

            else

                stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR13'

            end if cd0

        else

            stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR14'

        end if

        call UngetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR15'

        call UngetMap(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - ModuleAssimilation - ERR16'


    end subroutine ConstructPropertyCoefficients


    !--------------------------------------------------------------------------


    subroutine FindMinimumCoef(Minimum, Field2D, Field3D, WaterPoints2D, WaterPoints3D)

        !Arguments-------------------------------------------------------------
        real                                                :: Minimum
        real,       dimension(:,:  ), pointer, optional     :: Field2D
        real,       dimension(:,:,:), pointer, optional     :: Field3D
        integer,    dimension(:,:  ), pointer, optional     :: WaterPoints2D
        integer,    dimension(:,:,:), pointer, optional     :: WaterPoints3D

        !External--------------------------------------------------------------
        integer                                             :: i, j, k
        integer                                             :: ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 

        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 

        KLB = Me%WorkSize%KLB 
        KUB = Me%WorkSize%KUB

        Minimum = - FillValueReal


        !A constant value in all the domain is consider by default
cd1:    if (present(Field2D) .and. present(WaterPoints2D)) then 

            do j = JLB, JUB 
            do i = ILB, IUB 
                if (WaterPoints2D(i,j) == WaterPoint .and. Field2D(i, j) < Minimum) then

                    Minimum = Field2D(i, j)
                    
                endif
            enddo
            enddo

        else if (present(Field3D) .and. present(WaterPoints3D)) then 

            do k = KLB, KUB 
            do j = JLB, JUB 
            do i = ILB, IUB 

                if (WaterPoints3D(i,j,k) == WaterPoint .and. Field3D(i, j, k) < Minimum) then

                    Minimum = Field3D(i, j, k)
                    
                endif

            enddo
            enddo
            enddo

        else 

            stop 'FindMinimumCoef - ModuleAssimilation - ERR01'

        endif cd1

    end subroutine FindMinimumCoef


    !--------------------------------------------------------------------------

    subroutine ConstructOutPut

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------


        call GetOutPutTime(Me%ObjEnterData,                                              &
                           CurrentTime   = Me%ActualTime,                                &
                           EndTime       = Me%EndTime,                                   &
                           keyword       = 'OUTPUT_TIME',                                &
                           SearchType    = FromFile,                                     &
                           OutPutsTime   = Me%OutPut%OutTime,                            &
                           OutPutsOn     = Me%OutPut%ON,                                 &
                           OutPutsNumber = Me%OutPut%Number,                             &
                           STAT          = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructOutPut - ModuleAssimilation - ERR01'


        if (Me%OutPut%ON) Me%OutPut%Next = 1


    end subroutine ConstructOutPut

    !--------------------------------------------------------------------------
    ! This subroutine adds a new property to the Water Property List  

    subroutine Add_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),   pointer         :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstAssimilationProp)) then
            Me%PropertiesNumber = 1
            Me%FirstAssimilationProp     => NewProperty
            Me%LastAssimilationProp      => NewProperty
        else
            NewProperty%Prev             => Me%LastAssimilationProp
            Me%LastAssimilationProp%Next => NewProperty
            Me%LastAssimilationProp      => NewProperty
            Me%PropertiesNumber          = Me%PropertiesNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Property 


    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetAssimilationSize(AssimilationID, ILB, IUB, JLB, JUB, KLB, KUB, STAT)

        !Arguments--------------------------------------------------------------
        integer                         :: AssimilationID
        integer, optional, intent(OUT)  :: ILB, IUB 
        integer, optional, intent(OUT)  :: JLB, JUB 
        integer, optional, intent(OUT)  :: KLB, KUB
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_            

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(ILB)) ILB = Me%Size%ILB
            if (present(IUB)) IUB = Me%Size%IUB

            if (present(JLB)) JLB = Me%Size%JLB
            if (present(JUB)) JUB = Me%Size%JUB

            if (present(KLB)) KLB = Me%Size%KLB
            if (present(KUB)) KUB = Me%Size%KUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetAssimilationSize


    !--------------------------------------------------------------------------

    
    subroutine GetAssimilationField(AssimilationID, ID,                         & 
                                       Field2D, Field3D,                        &
                                       STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: AssimilationID
        real, dimension(:,:  ), pointer, optional   :: Field2D
        real, dimension(:,:,:), pointer, optional   :: Field3D
        integer                                     :: ID
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        
        type (T_Property), pointer                  :: PropertyX    
        integer                                     :: STAT_CALL, STAT_CALL1

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nullify(PropertyX)
            call SearchProperty(PropertyX = PropertyX, PropertyXIDNumber = ID, STAT = STAT_CALL) 
                                         
                           
cd2:        if (STAT_CALL == SUCCESS_) then


                STAT_CALL1 = SUCCESS_

cd3:            if (PropertyX%Dim == Dim_2D) then

                    if (present(Field2D)) then 

                        call Read_Lock(mASSIMILATION_, Me%InstanceID)

                        Field2D     => PropertyX%Field%R2D

                    endif


                    if (present(Field3D    )) STAT_CALL1 =  NOT_ASSOCIATE_


                else if (PropertyX%Dim == Dim_3D) then

                    if (present(Field3D)) then 

                        call Read_Lock(mASSIMILATION_, Me%InstanceID)

                        Field3D     => PropertyX%Field%R3D

                    endif


                    if (present(Field2D    )) STAT_CALL1 =  NOT_ASSOCIATE_

                endif cd3

                STAT_ = STAT_CALL1

            else  cd2

                STAT_ = STAT_CALL

            endif cd2

        else cd1
         
            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetAssimilationField

    !--------------------------------------------------------------------------

    subroutine GetAssimilationCoef(AssimilationID, ID, CoefField2D, CoefField3D,        &
                                   ColdRelaxPeriod, ColdOrder, Minimum, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: AssimilationID
        real, dimension(:,:  ), pointer, optional   :: CoefField2D
        real, dimension(:,:,:), pointer, optional   :: CoefField3D
        real,    optional, intent(OUT)              :: Minimum, ColdRelaxPeriod, ColdOrder
        integer          , intent(IN )              :: ID
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        
        type (T_Property), pointer                  :: PropertyX    
        integer                                     :: STAT_CALL, STAT_CALL1

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nullify(PropertyX)
            call SearchProperty(PropertyX = PropertyX, PropertyXIDNumber = ID , STAT = STAT_CALL)
                                            

cd2:        if (STAT_CALL == SUCCESS_) then

                if (present(ColdRelaxPeriod))                                           &
                    ColdRelaxPeriod = PropertyX%ColdRelaxPeriod

                if (present(ColdOrder))                                                 &
                    ColdOrder = PropertyX%ColdOrder

                STAT_CALL1 = SUCCESS_

cd3:            if (PropertyX%Dim == Dim_2D) then

                    if (present(CoefField2D)) then 

                        call Read_Lock(mASSIMILATION_, Me%InstanceID)

                        CoefField2D => PropertyX%CoefField%R2D

                    endif

                    if (present(CoefField3D)) STAT_CALL1 =  NOT_ASSOCIATE_


                else if (PropertyX%Dim == Dim_3D) then

                    if (present(CoefField3D)) then 

                        call Read_Lock(mASSIMILATION_, Me%InstanceID)

                        CoefField3D => PropertyX%CoefField%R3D

                    endif

                    if (present(CoefField2D)) STAT_CALL1 =  NOT_ASSOCIATE_

                endif cd3

                if (present(Minimum)) Minimum = PropertyX%CoefField%Minimum

                STAT_ = STAT_CALL1

            else  cd2

                STAT_ = STAT_CALL

            endif cd2

        else cd1
         
            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetAssimilationCoef

    !--------------------------------------------------------------------------
    ! guillaume nogueira

    subroutine GetAssimilationAltimetry (AssimilationID, AltimetricAssimilation, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: AssimilationID
        logical,        intent(OUT)                 :: AltimetricAssimilation
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            AltimetricAssimilation = Me%AltimetricAssimilation
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetAssimilationAltimetry

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! guillaume nogueira

    subroutine GetAssimilationAltimetryDT (AssimilationID, AltimetryDT, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: AssimilationID
        real,        intent(OUT)                    :: AltimetryDT
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            AltimetryDT = Me%Altimetric_Assim%DtAltimAssimilation
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetAssimilationAltimetryDT 

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    ! guillaume nogueira

    subroutine GetAltimetryDecayTime (AssimilationID, AltimetryDecayTime, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: AssimilationID
        real,        intent(OUT)                    :: AltimetryDecayTime
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            AltimetryDecayTime = Me%Altimetric_Assim%AltimDecayTime
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetAltimetryDecayTime 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! guillaume nogueira

    subroutine GetAltimSigmaDensAnalyzed (AssimilationID, SigmaDensAnalyzed, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: AssimilationID
        real, dimension(:,:,:), pointer             :: SigmaDensAnalyzed
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            call Read_Lock(mASSIMILATION_, Me%InstanceID)

            SigmaDensAnalyzed => Me%Altimetric_Assim%SigmaDensAnalyzed
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetAltimSigmaDensAnalyzed 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine UngetAssimilation2D(AssimilationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: AssimilationID
        real, pointer, dimension(:,:)       :: Array
        integer, optional, intent (OUT)     :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_   

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mASSIMILATION_, Me%InstanceID, "UngetAssimilation2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetAssimilation2D

    !--------------------------------------------------------------------------

    subroutine UngetAssimilation3D(AssimilationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: AssimilationID
        real(4), pointer, dimension(:,:,:)  :: Array
        integer, optional, intent (OUT)     :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_   

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mASSIMILATION_, Me%InstanceID, "UngetAssimilation3D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetAssimilation3D

    !--------------------------------------------------------------------------

    subroutine UngetAssimilation3D8(AssimilationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: AssimilationID
        real(8), pointer, dimension(:,:,:)  :: Array
        integer, optional, intent (OUT)     :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_   

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mASSIMILATION_, Me%InstanceID, "UngetAssimilation3D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetAssimilation3D8

    !--------------------------------------------------------------------------


    subroutine SearchProperty(PropertyX, PropertyXIDNumber, STAT)

        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer             :: PropertyX
        integer         ,           intent (IN)         :: PropertyXIDNumber
        integer         , optional, intent (OUT)        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_ 

        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstAssimilationProp

        do while (associated(PropertyX)) 
            if (PropertyX%ID%IDNumber==PropertyXIDNumber) then
                exit        
            else
                PropertyX => PropertyX%Next                 
            end if    
        end do    

        if (associated(PropertyX)) then

            call AssimilationFromFile(PropertyX)

            STAT_ = SUCCESS_  

        else
            STAT_  = NOT_FOUND_ERR_  
        end if


        if (present(STAT)) STAT = STAT_

    end subroutine SearchProperty


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!   Esta rotina é nova e calcula os campos de assimilação

    subroutine ModifyAssimilation(AssimilationID, Temperature, Salinity,                        &
                                    WaterLevelModel,                                            &
                                    DensMethod,                                                 &
                                    CorrecPress,                                                &
                                    STAT)                         !J. Nogueira Guillaume

        !Arguments-------------------------------------------------------------
        integer                                     :: AssimilationID
        real, dimension(:,:  ), pointer             :: WaterLevelModel
        real, dimension(:,:,:), pointer             :: Temperature, Salinity
        integer, optional, intent(OUT)              :: STAT
        integer                                     :: DensMethod
        Logical                                     :: CorrecPress
   
        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        real, dimension(:,:  ), pointer             :: WaterLevelToAssimilate, WaterLevelAnalyzed, &
                                                       VarianceFieldToAssimilate
        real, dimension(:,:,:), pointer             :: TemperatureAnalyzed, SalinityAnalyzed
        type(T_Property), pointer                   :: Property
        
        !Begin-----------------------------------------------------------------

       
        STAT_ = UNKNOWN_

        call Ready(AssimilationID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (MonitorPerformance) call StartWatch ("ModuleAssimilation", "ModifyAssimilation")

            call ReadLockExternalVar
         
            !Fetch Temperature Analyzed pointer
            call SearchProperty(PropertyX = Property,                                    &
                                PropertyXIDNumber = AltimTemperatureAnalyzed_,           &
                                STAT = STAT)
            if (STAT /= SUCCESS_) stop 'ModifyAssimilation - ModuleAssimilation - ERR01'
            TemperatureAnalyzed => Property%Field%R3D

            !Fetch Salinity Analyzed pointer
            call SearchProperty(PropertyX = Property,                                    &
                                PropertyXIDNumber = AltimSalinityAnalyzed_,              &
                                STAT = STAT)
            if (STAT /= SUCCESS_) stop 'ModifyAssimilation - ModuleAssimilation - ERR02'
            SalinityAnalyzed => Property%Field%R3D

            !Fetch Water Level Analyzed pointer
            call SearchProperty(PropertyX = Property,                                    &
                                PropertyXIDNumber = AltimLevelAnalyzed_,                 &
                                STAT = STAT)
            if (STAT /= SUCCESS_) stop 'ModifyAssimilation - ModuleAssimilation - ERR03'
            WaterLevelAnalyzed => Property%Field%R2D

            !Fetch Water Level to assimilate pointer
            call SearchProperty(PropertyX = Property,                                    &
                                PropertyXIDNumber = AltimLevelToAssimilate_,             &
                                STAT = STAT)
            if (STAT /= SUCCESS_) stop 'ModifyAssimilation - ModuleAssimilation - ERR03'
            WaterLevelToAssimilate => Property%Field%R2D

            !Guillaume: Testar se e para se usar o campo de variancias
            if (Me%Altimetric_Assim%UseVarianceField) then

                !Fetch Variance Field to assimilate pointer
                call SearchProperty(PropertyX = Property,                                    &
                                PropertyXIDNumber = VarianceFieldToAssimilate_,          &
                                STAT = STAT)
                if (STAT /= SUCCESS_) stop 'ModifyAssimilation - ModuleAssimilation - ERR04'
                VarianceFieldToAssimilate => Property%Field%R2D
            
            end if

            !Cooper-Haines method
            call AltimetricAssimilation (   Temperature, Salinity, WaterLevelModel,      &
                                            WaterLevelToAssimilate,                      &
                                            VarianceFieldToAssimilate,                   &
                                            WaterLevelAnalyzed,                          &
                                            TemperatureAnalyzed,                         &
                                            SalinityAnalyzed,                            &
                                            DensMethod,                                  &
                                            CorrecPress,                                 &
                                            STAT)
           if (STAT /= SUCCESS_) stop 'ModifyAssimilation - ModuleAssimilation - ERR01'

            call ReadUnlockExternalVar

            if (MonitorPerformance) call StopWatch ("ModuleAssimilation", "ModifyAssimilation")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

    end subroutine ModifyAssimilation

! ----------------------------------------------------------------------------------!
! ----------------------------------------------------------------------------------!
! ----------------------------------------------------------------------------------!
! ----------------------------------------------------------------------------------!
!                                                                                   !
!      This routine assimilate altimetric data with Cooper and Haines scheme           !
!                                                                                   !
!   Cooper,M. and K,Haines, 1996.  "Altimetric assimilation with water property     !
!                                conservation". Journal of Geophysical Research,    !
!                                   101, C1, 1059-1077.                                ! 
!                                                                                    !
!   The method is modified according to Drakopoutos et.al.                            !
!                                                                                    !
!   Drakopoulos, P.et.al. 1997. "Altimetric assimilation in a mediterranean general !
!                              circulation model". Journal of Geophysical Research, !
!                                 102, C1, 10509-10523.                                ! 
!                                                                                    !
!    inputs: WaterLevelModel (SZZ(KUB))                                              !
!           Temperature                                                             !
!           Salinity                                                                !
!                                                                                   !
!           WaterLevelToAssimilate (from file)                                      !
!           VarianceFieldToAssimilate (from file)                                   !
!                                                                                   !
!                                                                                    !
!    output: WaterLevelAnalized                                                      !
!           TemperatureAnalyzed                                                     !
!           SalinityAnalyzed                                                        !
!                                                                                    !
!   2005 -   João Nogueira    joaonogueira@ist.utl.pt                                 !
!                                                                                    !
! ----------------------------------------------------------------------------------!
! ----------------------------------------------------------------------------------!
! ----------------------------------------------------------------------------------!
! ----------------------------------------------------------------------------------!

    Subroutine AltimetricAssimilation ( Temperature, Salinity, WaterLevelModel, WaterLevelToAssimilate,      &
                                        VarianceFieldToAssimilate, WaterLevelAnalyzed, TemperatureAnalyzed,  &
                                        SalinityAnalyzed, DensMethod, CorrecPress, STAT)


        Implicit none


       !Arguments-------------------------------------------------------------
        real,       dimension(:,:,:),   pointer       :: Temperature, Salinity
        real,       dimension(:,:),     pointer       :: WaterLevelToAssimilate
        real,       dimension(:,:),     pointer       :: VarianceFieldToAssimilate
        real,       dimension(:,:),     pointer       :: WaterLevelAnalyzed, WaterLevelModel
        real,       dimension(:,:,:),   pointer       :: TemperatureAnalyzed, SalinityAnalyzed
        integer, optional                             :: STAT
        integer                                       :: DensMethod
        Logical                                       :: CorrecPress

        !External--------------------------------------------------------------

        real,       dimension(:,:,:), pointer       :: ZCellCenter
        integer,    dimension(:,:,:), pointer       :: OpenPoints3D
        real,       dimension(:,:,:), pointer       :: DWZ
        real,       dimension(:,:,:), pointer       :: SZZ
        integer,    dimension(:,:), pointer         :: KFloor_Z
        real,       dimension(:,:), pointer         :: CellArea
        integer                                     :: IUB, JUB, KUB, ILB, JLB, KLB

        !Local-----------------------------------------------------------------

        Integer i, j, k, kspl, kmin, kbottom


        Real    Mean_WaterLevel_Diference,                                                  &
                DeltaPressureSurface, ColumnPressure, DepthLevel, Total_Area,               &
                Mean_WaterLevelToAssimilate 
        Real    DepthLimit

        Real,Pointer,Dimension (:)::     ColumnTemperature,                                 &
                                         ColumnSalinity,                                    &
                                         ColumnDensity,                                     &
                                         DeltaTemperature,                                  &
                                         DeltaSalinity,                                     & 
                                         DeltaDensity,                                      &
                                         DeltaPressure

        Real,Pointer,Dimension (:,:)::   Delta_WaterLevel,                                  &
                                         Delta_WaterLevelToAssimilate,                      &
                                         PressureAnomaly

        Real,Pointer,Dimension (:,:,:):: SigmaDensAnalyzed

        Real,Pointer,Dimension (:,:)::   Observation_Error,                                 &
                                         Model_Error,                                       &
                                         Gain

        !Parameter-----------------------------------------------------------------

        Real, Parameter                             ::ModelError        =    1.0
        Real, Parameter                             ::ErrorMax          =    0.8
        Real, Parameter                             ::ReferenceDensity  = 1023.0
        Real, Parameter                             ::Gravity           =    9.806


        !Begin Shorten variable names

        IUB  = Me%WorkSize%IUB 
        JUB  = Me%WorkSize%JUB 
        KUB  = Me%WorkSize%KUB 

        ILB  = Me%WorkSize%ILB 
        JLB  = Me%WorkSize%JLB 
        KLB  = Me%WorkSize%KLB 


        ZCellCenter                     => Me%ExternalVar%ZCellCenter
        OpenPoints3D                    => Me%ExternalVar%OpenPoints3D
        KFloor_Z                        => Me%ExternalVar%KFloor_Z
        CellArea                        => Me%ExternalVar%GridCellArea
        DWZ                             => Me%ExternalVar%DWZ
        SZZ                             => Me%ExternalVar%SZZ

        ColumnTemperature               => Me%Altimetric_Assim%ColumnTemperature
        ColumnSalinity                  => Me%Altimetric_Assim%ColumnSalinity
        ColumnDensity                   => Me%Altimetric_Assim%ColumnDensity
        DeltaTemperature                => Me%Altimetric_Assim%DeltaTemperature
        DeltaSalinity                   => Me%Altimetric_Assim%DeltaSalinity
        DeltaDensity                    => Me%Altimetric_Assim%DeltaDensity
        DeltaPressure                   => Me%Altimetric_Assim%DeltaPressure
        Delta_WaterLevel                => Me%Altimetric_Assim%Delta_WaterLevel
        Delta_WaterLevelToAssimilate    => Me%Altimetric_Assim%Delta_WaterLevelToAssimilate
        PressureAnomaly                 => Me%Altimetric_Assim%PressureAnomaly
        SigmaDensAnalyzed               => Me%Altimetric_Assim%SigmaDensAnalyzed

        Observation_Error               => Me%Altimetric_Assim%Observation_Error
        Model_Error                     => Me%Altimetric_Assim%Model_Error
        Gain                            => Me%Altimetric_Assim%Gain

        DepthLimit                      = Me%Altimetric_Assim%AltimetricDepth



!   2D variables inicialization

        do i=ILB, IUB
        do j=JLB, JUB

            WaterLevelAnalyzed(i,j) = WaterLevelModel(i,j)

        enddo
        enddo

        do i=ILB, IUB
        do j=JLB, JUB

    !        if (present(VarianceFieldToAssimilate)) then
    !       Guillaume :
           if(Me%Altimetric_Assim%UseVarianceField) then
                Observation_Error(i,j) = VarianceFieldToAssimilate(i,j)
                Model_Error(i,j) = ModelError
            else
                Observation_Error(i,j) = 0.
                Model_Error(i,j) = ModelError
           endif

        enddo
        enddo
    
!   3D variables inicialization

        do k=KLB, KUB
        do j=JLB, JUB
        do i=ILB, IUB

            TemperatureAnalyzed(i,j,k) = Temperature(i,j,k)

        enddo
        enddo
        enddo
    
        do k=KLB, KUB
        do j=JLB, JUB
        do i=ILB, IUB

            SalinityAnalyzed(i,j,k) = Salinity(i,j,k)

        enddo
        enddo
        enddo




!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!
!                   MODEL COMPUTE SURFACE HEIGTH (PRESSURE) ANOMALIE FIELD
!
!        The water level correction field result from: 
!
!        Gain * [(WaterLevelToAssimilate - WaterLevelModel) - Mean_WaterLevel_Diference]
!
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------



!
!       The Gain matrix represents the weigth given too diference between model and 
!   altimetric data. Take values between 0 and 1. 
!
!     If Observation_Error=0: Analized level = Assimilated level
!     If Model_Error=0: Analized level = Model level (no correction is made)
!
!

!            do j=JLB, JUB
!            do i=ILB, IUB
!              Gain(i,j) = Model_Error(i,j)
!            enddo
!            enddo

            do j=JLB, JUB
            do i=ILB, IUB
                Gain(i,j) = Model_Error(i,j)/(Model_Error(i,j) + Observation_Error(i,j))
            enddo
            enddo


!   Compute diferences between mesure and model taking into account the Gain

        do j=JLB, JUB
        do i=ILB, IUB
            Delta_WaterLevel(i,j) = WaterLevelToAssimilate(i,j) - WaterLevelModel(i,j)
            Delta_WaterLevel(i,j) = Gain(i,j) * Delta_WaterLevel(i,j)
        enddo
        enddo

!
!       Compute average diference between model and mesures including only points with  
!   reliable data.The criteria is: open points (points with water) and data with small 
!   error.
!

        Total_Area = 0.
        Mean_WaterLevelToAssimilate = 0.
        do i=ILB, IUB
        do j=JLB, JUB
            if(OpenPoints3D(i,j,KUB).eq.OpenPoint .and. Observation_Error(i,j).lt.ErrorMax)then
                Total_Area = Total_Area + CellArea(i,j)
                Mean_WaterLevel_Diference = Mean_WaterLevel_Diference                       &
                                           + Delta_WaterLevel(i,j) * CellArea(i,j) 
            endif
        enddo
        enddo

        if (Total_Area .ne. 0.) then
            Mean_WaterLevel_Diference = Mean_WaterLevel_Diference / Total_Area
        else
            write(*,*) 'Error: There are no waterpoints in your domain!'
            stop 'AltimetricAssimilation - ModuleAssimilation - ERR06'
        endif


!
!   The anomaly value (Delta_WaterLevel - Mean_WaterLevel_Diference) is converted to 
!   pressure anomaly
!

        do i=ILB, IUB
        do j=JLB, JUB
            if(OpenPoints3D(i,j,KUB).eq.OpenPoint .and. Observation_Error(i,j).lt.ErrorMax)then
                PressureAnomaly(i,j) = ReferenceDensity * Gravity *                         &
                                      (Delta_WaterLevel(i,j) - Mean_WaterLevel_Diference)
            else
                PressureAnomaly(i,j) = 0.
            endif
        enddo
        enddo

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!
!                           MODEL COMPUTE VERTICAL PROJECTION
!
!          In order to conserve mass and potencial vorticity, the model changes T, S,  
!      distribution according to pressure anomaly at the surface (the bottom pressure 
!      is unchanged). This changes are computed only below the thermocline 
!      (DepthLimit = 50m):
!
!       Fox et. al., 2000.  "Altimeter assimilation in the OCCAM global model Part I:  
!                           A twin experiment". Journal of Marine Systems, 26, 303-322
!
!       Fox et. al., 2000.  "Altimeter assimilation in the OCCAM global model Part II:  
!                           TOPEX/POSEIDON and ERS-1 assimilation". Journal of Marine
!                           Systems, 26, 323-347
!
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------


!
!   ZCellCenter is the center level for cell k (negative downward)
!
!   SZZ(i,j,k) is the up face level for cell k (positive downward):
!
!           --->  SZZ(i,j,0) is the bathimetry.
!           ---> -SZZ(i,j,KUB) is the water level.
!
!   Depth in the center of cell k is defined by:
!
!         Depth(i,j,k) = -1.0*[SZZ(i,j,KUB)+ZCellCenter(i,j,k)]
!


do1 :   do i=ILB, IUB
do2 :   do j=JLB, JUB

            kbottom = KFloor_Z(i,j)      ! This assumption enssure that bottom is not moved

            kmin = kbottom + 1        

            kspl = kbottom
    !
    !       kspl define the k value until where correction is made at each (i,j).
    !
            if(OpenPoints3D(i,j,KUB) .eq. OpenPoint)then

                do k=kbottom, KUB
                    DepthLevel = -1.0*(SZZ(i,j,KUB)+ZCellCenter(i,j,k))
                    if (DepthLevel.gt.DepthLimit) then
                        kspl = k
                    else
                        exit
                    endif
                enddo

            endif

    !       The model makes vertical projection if there are more than one level (with
    !       water) below the thermocline.

            if(OpenPoints3D(i,j,KUB) .eq. OpenPoint .and. kmin .lt. kspl)then

    !           Inicialization values for pressure

                DeltaPressureSurface = PressureAnomaly(i,j)
                ColumnPressure = PressureAnomaly(i,j)

    !           Set the tracer values for the kspl-1 top levels to value at kspl. This 
    !           is to define the structure of the vertjcaly displaced water column. Set
    !           the rest equal to original values.

                do k=kbottom, KUB
                    if (k .le. kspl) then
                        ColumnTemperature(k) = Temperature(i,j,k)
                        ColumnSalinity(k) = Salinity(i,j,k)
                    else
                        ColumnTemperature(k) = Temperature(i,j,kspl)
                        ColumnSalinity(k) = Salinity(i,j,kspl)
                    endif
                enddo

                Call LoweringLifting(ColumnPressure, DeltaPressureSurface,                  &
                                     ColumnTemperature, ColumnSalinity, ColumnDensity,      &
                                     DeltaTemperature,                                      &
                                     DeltaSalinity, DeltaDensity, DeltaPressure,            & 
                                     SZZ, ZCellCenter, DWZ, kbottom,                        &
                                     kspl, kmin, i, j,                                      &
                                     DensMethod,                                            &
                                     CorrecPress)

                PressureAnomaly(i,j) = ColumnPressure

    !           Update temperature and salinity fields

                do k=kbottom+1, kspl

                    TemperatureAnalyzed(i,j,k) = Temperature(i,j,k) +              &
                                                   DeltaTemperature(k)
                    SalinityAnalyzed(i,j,k) = Salinity(i,j,k) + DeltaSalinity(k)

                    DepthLevel = -1.0*(SZZ(i,j,KUB) + ZCellCenter(i,j,k))
                    
                    SigmaDensAnalyzed(i,j,k) = Sigma (DensMethod,                   &
                                                      CorrecPress,                  &
                                                      TemperatureAnalyzed(i,j,k),   &
                                                      SalinityAnalyzed(i,j,k),      &
                                                      DepthLevel)
                enddo

            else  ! For this grid point (i,j) model keep the original value

                PressureAnomaly(i,j) = 0.
                Delta_WaterLevelToAssimilate(i,j) = 0.

            endif

        enddo do2
        enddo do1

    !   Update Surface Height

        do j=JLB, JUB
        do i=ILB, IUB

            WaterLevelAnalyzed(i,j) = WaterLevelModel(i,j) + (Delta_WaterLevel(i,j)     &
                                      - Mean_WaterLevel_Diference)

        enddo
        enddo

        STAT = SUCCESS_

    End Subroutine AltimetricAssimilation
    !----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!
!                           MODEL COMPUTE VERTICAL DISPLACMENT
!
!
!          "Updates to the surface height field are projected into the interior using a   
!       simple lifting/lowering scheme. The entire water column is displaced vertically, 
!       conserving water properties temperature, salinity and potential vorticity in 
!       isopycnal layers."
!
!       Fox et. al., 2000. "Altimeter assimilation in the OCCAM global model Part II:  
!                           TOPEX/POSEIDON and ERS-1 assimilation". Journal of Marine
!                           Systems, 26, 323-347
!
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------



    Subroutine LoweringLifting (ColumnPressure, DeltaPressureSurface,                   &
                            ColumnTemperature, ColumnSalinity, ColumnDensity,           &
                            DeltaTemperature,                                           &
                            DeltaSalinity, DeltaDensity, DeltaPressure,                 & 
                            SZZ, ZCellCenter, DWZ, kbottom,                             &
                            kspl, kmin, i, j,                                           &
                            DensMethod,                                                 &
                            CorrecPress)

        Implicit none

        !Arguments--------------------------------------------------------------
        Real    DeltaPressureBottom, DeltaPressureSurface,                                  &
                ColumnPressure, DepthLevel,                                                 &
                D_DeltaThermocline,                                                         &
                MaxDisplacement

        logical last                            

        Real,Pointer,Dimension (:)::    ColumnTemperature,                                  &
                                        ColumnSalinity,                                     &
                                        ColumnDensity,                                      &
                                        DeltaTemperature,                                   &
                                        DeltaSalinity,                                      & 
                                        DeltaDensity,                                       &
                                        DeltaPressure

        real,       dimension(:,:,:), pointer           :: ZCellCenter
        real,       dimension(:,:,:), pointer           :: DWZ
        real,       dimension(:,:,:), pointer           :: SZZ
        integer                                         :: DensMethod
        Logical                                         :: CorrecPress


        !External--------------------------------------------------------------
        Integer IUB, JUB, KUB, ILB, JLB, KLB

        !Local-----------------------------------------------------------------


        Integer niter, i, j, k, kspl, kmin, kbottom
        Real  DeltaThermocline, DeltaPressure_correction

        !Parameter-----------------------------------------------------------------

        Real, Parameter                             ::PressureError     =    1.0
        Real, Parameter                             ::ReferenceDensity  = 1023.0
        Real, Parameter                             ::Gravity           =    9.806


        !Begin Shorten Variables---------------------------------------------------

        IUB  = Me%WorkSize%IUB 
        JUB  = Me%WorkSize%JUB 
        KUB  = Me%WorkSize%KUB 

        ILB  = Me%WorkSize%ILB 
        JLB  = Me%WorkSize%JLB 
        KLB  = Me%WorkSize%KLB 


    !   Inicialize updates

        do k=kmin, KUB

            DeltaPressure(k) = 0.
            DeltaSalinity(k) = 0.
            DeltaTemperature(k) = 0.
            DeltaDensity(k) = 0.

        enddo

    !   Find in-situ density

        do k=kmin, kspl

            DepthLevel = -1.0*(SZZ(i,j,KUB) + ZCellCenter(i,j,k))
!            ColumnDensity(k) = Density (UNESCOState_, .true., ColumnTemperature(k), ColumnSalinity(k), DepthLevel)
            ColumnDensity(k) = Density (DensMethod,                                             &
                                        CorrecPress,                                            &
                                        ColumnTemperature(k),                                   &
                                        ColumnSalinity(k),                                      &
                                        DepthLevel)

        enddo

    !   Initially, assume DeltaThermocline equals zero, and so surface pressure
    !   anomaly penetrates straight to bottom

        DeltaThermocline = 0.
        DeltaPressureBottom = DeltaPressureSurface

        if (kspl .eq. KUB) then
            D_DeltaThermocline = (abs(SZZ(i,j,kmin)))/10.
        else
            D_DeltaThermocline = (SZZ(i,j,kmin) - SZZ(i,j,kspl+1))/10.
        endif

        D_DeltaThermocline = min(D_DeltaThermocline,350.0)

!-------------------------------------------------------------------------------   
!-------------------------------------------------------------------------------   
!   
!                       START ITERATIVE PROCEDURE. 
!   
!   
!       The aim is to have:
!   
!                       DeltaPressureBottom = 0
!   
!       The iteration stops if (DeltaPressureBottom .le. PressureError)
!   
!            D_DeltaThermocline is the increment made to DeltaThermocline at each  
!      iteration.It has initial value equal to the ocean depth. On each iteration, 
!      D_DeltaThermocline is either added or subtracted, depending on the sign of 
!      the botom pressure variation, and then halved for the next iteration.  
!      So after 20 iterations, the accuracy of the D_DeltaThermocline is: 
!   
!           (ocean depth)/(2**20)
!   
!      about 3 mm for a 4 km column.
!      This routine never diverges.
!   
!-------------------------------------------------------------------------------   
!-------------------------------------------------------------------------------   


        last = .false.
        niter = 0


        if (ColumnPressure .eq. 0.) then

            do k=kmin, KUB

                DeltaPressure(k) = 0.
                DeltaSalinity(k) = 0.
                DeltaTemperature(k) = 0.

            enddo

        else

            do while (niter .le. 25)

                niter = niter + 1

                if (DeltaPressureBottom .gt. 0.)                                            &
                    DeltaThermocline = DeltaThermocline - D_DeltaThermocline

                if (DeltaPressureBottom .lt. 0.)                                            &
                    DeltaThermocline = DeltaThermocline + D_DeltaThermocline

    !           Model enssure that no modification is made at surface

                if (kspl .eq. KUB) then
                    MaxDisplacement = SZZ(i,j,kbottom) + DWZ(i,j,kbottom)
                else
                    MaxDisplacement = SZZ(i,j,kbottom) + DWZ(i,j,kbottom) - SZZ(i,j,kspl)
                endif
                    MaxDisplacement = min(MaxDisplacement,350.001)

                if (abs(DeltaThermocline).gt.MaxDisplacement .or.                           &
                    DeltaPressureBottom.eq.0.) then

                    last = .true.
                    DeltaThermocline = 0.

                    do k=kmin, KUB

                        DeltaPressure(k) = 0.
                        DeltaSalinity(k) = 0.
                        DeltaTemperature(k) = 0.

                    enddo

                else

                    Call ChangeTemperatureSalinity(ColumnPressure,    &
                                         ColumnTemperature, ColumnSalinity, ColumnDensity,  &
                                         DeltaTemperature, DeltaSalinity, DeltaDensity,     & 
                                         DeltaPressure, DeltaThermocline,                   &
                                         SZZ, ZCellCenter, DWZ, kspl, kmin, i, j)

                endif

                DeltaPressureBottom = DeltaPressure(kmin)
                D_DeltaThermocline = 0.5 * D_DeltaThermocline

                if  (abs(DeltaPressureBottom) .lt. PressureError .or. last) exit

            enddo

    !       If iterative procedure do not converge no modification is make at this grid point

            if (DeltaPressureBottom .ge. PressureError) then

                DeltaThermocline = 0.
                ColumnPressure = 0.

                do k=kmin, kspl

                    DeltaPressure(k) = 0.
                    DeltaSalinity(k) = 0.
                    DeltaTemperature(k) = 0.

                enddo

            endif

    !       If botom pressure variation is not zero, model set it to zero and adjust all 
    !       pressures above similarly.

            DeltaPressure_correction = DeltaPressure(kmin)

            do k=kspl,kmin,-1
              DeltaPressure(k) = DeltaPressure(k) - DeltaPressure_correction
            enddo

            ColumnPressure = ColumnPressure - DeltaPressure_correction

        endif


    End Subroutine LoweringLifting

    !--------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!
!                           MODEL COMPUTE (T,S) PROFILE CORRECTION
!
!
!       Now we have to modified the profile of temperature and salinity in agremment with
!       the vertical displacement. This procedure is made with an interpolation of (T,S) 
!       in all column.
!
!
!       Cooper,M. and Haines,K, 1996.  "Altimetric assimilation with water property
!                                    conservation". Journal of Geophysical Research,
!                                       101, C1, 1059-1077.             
!
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------


    Subroutine ChangeTemperatureSalinity(ColumnPressure,           &
                                     ColumnTemperature, ColumnSalinity, ColumnDensity,  &
                                     DeltaTemperature, DeltaSalinity, DeltaDensity,     & 
                                     DeltaPressure, DeltaThermocline,                   &
                                     SZZ, ZCellCenter, DWZ, kspl, kmin, i, j)



        !Arguments-------------------------------------------------------------------
        Real                                    ::  ColumnPressure,    &
                                                    DeltaThermocline

        Real,Pointer,Dimension (:)              ::  ColumnTemperature, ColumnSalinity,              &
                                                    ColumnDensity, DeltaTemperature, DeltaSalinity, &
                                                    DeltaDensity, DeltaPressure

        real,       dimension(:,:,:), pointer       :: ZCellCenter
        real,       dimension(:,:,:), pointer       :: DWZ
        real,       dimension(:,:,:), pointer       :: SZZ

        Integer k, kmin, kspl, i, j

        !Local-------------------------------------------------------------------------
        Real,Pointer,Dimension (:)                  ::  Depth,                          &
                                                        D_GradTemp, D_GradSalin
        Real, Parameter                             ::  ReferenceDensity  = 1023.0
        Real, Parameter                             ::  Gravity           =    9.806
        integer                                     ::  KUB 
        real                                        ::  GradTemp_L, GradTemp_U,         &
                                                        GradSalin_L, GradSalin_U,       &
                                                        MinDepth, MaxDepth, Z_n,        &
                                                        Temperature_n, Salinity_n,      &
                                                        DepthLevel, ColumnDensity_n

        KUB             = Me%WorkSize%KUB 

        Depth           => Me%Altimetric_Assim%AuxDepth
        D_GradTemp      => Me%Altimetric_Assim%D_GradTemp
        D_GradSalin     => Me%Altimetric_Assim%D_GradSalin


        do k=kmin, KUB

            Depth(k) = -1.0*(SZZ(i,j,KUB) + ZCellCenter(i,j,k))

        enddo


        GradTemp_L = 1.e31
        GradTemp_U = 1.e31

        GradTemp_L = (ColumnTemperature(kmin+1) - ColumnTemperature(kmin))/                 &
                     (Depth(kmin+1) - Depth(kmin))

        GradTemp_U = (ColumnTemperature(KUB-2) - ColumnTemperature(KUB-1))/                   &
                     (Depth(KUB-2) - Depth(KUB-1))

        Call Spline_Alt (Depth, ColumnTemperature, kmin, GradTemp_U, GradTemp_L, D_GradTemp)

        GradSalin_L = 1.e31
        GradSalin_U = 1.e31

        GradSalin_L = (ColumnSalinity(kmin+1) - ColumnSalinity(kmin))/                      &
                     (Depth(kmin+1) - Depth(kmin))

        GradSalin_U = (ColumnSalinity(KUB-2) - ColumnSalinity(KUB-1))/                        &
                     (Depth(KUB-2) - Depth(KUB-1))
    
        Call Spline_Alt (Depth, ColumnSalinity, kmin, GradSalin_U, GradSalin_L, D_GradSalin)



        MinDepth = -1.0*(SZZ(i,j,KUB) + ZCellCenter(i,j,kspl))
        MaxDepth = -1.0*(SZZ(i,j,KUB) + ZCellCenter(i,j,kmin))

        do k=kmin, kspl

            Z_n = (-1.0*(SZZ(i,j,KUB) + ZCellCenter(i,j,k))) + DeltaThermocline
 
            if (Z_n .le. MinDepth) then
                Temperature_n = ColumnTemperature(kspl)
                Salinity_n = ColumnSalinity(kspl)
            else
                if (Z_n .ge. MaxDepth) then
                    Temperature_n = ColumnTemperature(kmin)
                    Salinity_n = ColumnSalinity(kmin)
                else

                    Call Splint (Depth, ColumnTemperature, D_GradTemp, kmin, -Z_n,           &
                                 Temperature_n)
                    Call Splint (Depth, ColumnSalinity, D_GradSalin, kmin, -Z_n, Salinity_n)

                endif
            endif

          DeltaTemperature(k) = Temperature_n - ColumnTemperature(k)
          DeltaSalinity(k) = Salinity_n - ColumnSalinity(k)

          DepthLevel = -1.0*(SZZ(i,j,KUB) + ZCellCenter(i,j,k))
          ColumnDensity_n = Density (UNESCOState_, .true., Temperature_n, Salinity_n, DepthLevel)

          DeltaDensity(k) = ColumnDensity_n - ColumnDensity(k)

        enddo

        if (kspl .eq. 1) then
          DeltaPressure(kspl) = ColumnPressure +                                            &
                                Gravity * DeltaDensity(kspl) * DWZ(i,j,kspl)
        else
          DeltaPressure(kspl) = ColumnPressure +                                            &
                                0.5 * Gravity * DeltaDensity(kspl) * DWZ(i,j,kspl)
        endif

    !   Makes the integral of Gravity*DeltaDensity to obtain DeltaPressureBottom

        do k=kspl-1,kmin,-1
          DeltaPressure(k) = DeltaPressure(k+1) + 0.5 * Gravity *                           &
                             (DeltaDensity(k) + DeltaDensity(k+1)) * DWZ(i,j,k)
        enddo


    End Subroutine ChangeTemperatureSalinity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    Subroutine Spline_Alt (Depth, ColumnProperty, kmin, GradProp_U, GradProp_L, D_GradProp)

        !Arguments-------------------------------------------------------------
        real                                    :: GradProp_U, GradProp_L
        real, dimension(:), pointer             :: Depth, ColumnProperty,  D_GradProp
        integer                                 :: kmin

        !Local-----------------------------------------------------------------
        Real  sig, p, D_Grad
        Integer n, KUB
        real, dimension(:), pointer             :: u

        KUB  = Me%WorkSize%KUB 
        u => Me%Altimetric_Assim%AuxSpline 
        
        Depth=-1.*Depth       

         if (GradProp_L.gt..99e30) then
            D_GradProp(kmin) = 0.
            u(kmin) = 0.
          else
            D_GradProp(kmin) = -0.5
            u(kmin) = (3./(Depth(kmin+1) - Depth(kmin))) * ((ColumnProperty(kmin+1) -           &
                      ColumnProperty(kmin))/(Depth(kmin+1) - Depth(kmin)) - GradProp_U)
          endif

          do n=kmin+1,KUB-1
            sig = (Depth(n) - Depth(n-1))/(Depth(n+1) - Depth(n-1))
            p = sig * D_GradProp(n-1) + 2.
            D_GradProp(n) = (sig - 1.)/p
            u(n) = (6. * ((ColumnProperty(n+1) - ColumnProperty(n))/(Depth(n+1) -           &
                   Depth(n)) - (ColumnProperty(n) - ColumnProperty(n-1))/                   &
                   (Depth(n) - Depth(n-1)))/(Depth(n+1) - Depth(n-1)) - sig * u(n-1))/p
          enddo

          if (GradProp_U.gt..99e30) then
            D_Grad = 0.
            u(KUB) = 0.
          else
            D_Grad = -0.5
            u(KUB) = (3./(Depth(KUB-1) - Depth(KUB))) * ((ColumnProperty(KUB-1) -       &
                      ColumnProperty(KUB))/(Depth(KUB-1) - Depth(KUB)) - GradProp_U)
          endif

          D_GradProp(KUB)=(u(KUB)-D_Grad*u(KUB-1))/(D_Grad*D_GradProp(KUB-1)+1.)

          do n=KUB-1,kmin,-1

             D_GradProp(n) = D_GradProp(n) * D_GradProp(n+1) + u(n)

          enddo


        Depth=-1.*Depth       

    End Subroutine Spline_Alt

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    Subroutine Splint (Depth, ColumnProperty, D_GradProp, kmin, Z_n, Property_U)

        !Arguments-------------------------------------------------------------
        real                                    :: Z_n, Property_U
        real, dimension(:), pointer             :: Depth, ColumnProperty, D_GradProp
        integer                                 :: kmin

        !Local------------------------------------------------------------------
        Real  h, a, b
        integer KUB

        Integer klo, khi, k

        Depth=-1.*Depth       

        KUB = Me%WorkSize%KUB

          klo=kmin
          khi=KUB
          do while (khi-klo.gt.1) 
            k=(khi+klo)/2
            if(Depth(k).gt.Z_n)then
              khi=k
            else
              klo=k
            endif
          enddo
          h=Depth(khi)-Depth(klo)
          if (h.eq.0.) then 
            write(*,*) 'Error: bad Depth input'
            stop 'Splint - Module Assimilation - ERR33.'
          end if
          a=(Depth(khi)-Z_n)/h
          b=(Z_n-Depth(klo))/h
          Property_U=a*ColumnProperty(klo)+b*ColumnProperty(khi)+                           &
            ((a**3-a)*D_GradProp(klo)+(b**3-b)*D_GradProp(khi))*(h**2)/6.

        Depth=-1.*Depth       


    End Subroutine Splint

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AssimilationFromFile(PropertyX)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer               :: PropertyX
        
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, i, j, k
        integer, dimension(:,:  ), pointer      :: WaterPoints2D, WaterFaces2D_U, WaterFaces2D_V, PointsToFill2D
        integer, dimension(:,:,:), pointer      :: WaterPoints3D, WaterFaces3D_U, WaterFaces3D_V, PointsToFill3D 
        real,    dimension(:,:  ), pointer      :: Matrix2D
        real,    dimension(:,:,:), pointer      :: Matrix3D
        integer                                    :: CHUNK
         
        !Begin------------------------------------------------------------------------

        call GetComputeCurrentTime(Me%ObjTime, Me%ActualTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'
        
        call GetWaterPoints3D(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

        call GetWaterFaces2D(Me%ObjHorizontalMap,                                               &
                               WaterFaces2DU = WaterFaces2D_U,      &
                               WaterFaces2DV = WaterFaces2D_V,      &
                               STAT            = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR90'

        call GetWaterFaces3D(Me%ObjMap,                                               &
                               WaterFacesU3D = WaterFaces3D_U,                      &
                               WaterFacesV3D = WaterFaces3D_V,                      &
                               STAT            = STAT_CALL)

        if (MonitorPerformance) then
            call StartWatch ("ModuleAssimilation", "AssimilationFromFile")
        endif

        !Verifies if its necessary to update the property
cd1:    if (Me%ActualTime > PropertyX%LastActualization) then


            if (PropertyX%ID%SolutionFromFile) then


                if (PropertyX%Dim == Dim_2D) then 

                    if (PropertyX%Field%TypeZUV == TypeZ_) then

                        PointsToFill2D => WaterPoints2D

                        if (GetPropertyIDNumber(PropertyX%ID%Name) == BarotropicVelocityU_ .or. &
                            GetPropertyIDNumber(PropertyX%ID%Name) == BarotropicVelocityV_) then

                            allocate (Matrix2D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                            Matrix2D(:,:) = FillValueReal
                        
                        else

                            Matrix2D => PropertyX%Field%R2D

                        endif

                    else if (PropertyX%Field%TypeZUV == TypeU_) then

                        PointsToFill2D => WaterFaces2D_U
                        Matrix2D       => PropertyX%Field%R2D

                    else if (PropertyX%Field%TypeZUV == TypeV_) then

                        PointsToFill2D => WaterFaces2D_V
                        Matrix2D       => PropertyX%Field%R2D

                    else

                        stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR230'

                    endif




                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,         &
                                           Matrix2D       = Matrix2D,                            &
                                           PointsToFill2D = PointsToFill2D,                      &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

                    if (PropertyX%Field%TypeZUV == TypeZ_) then
                        if (GetPropertyIDNumber(PropertyX%ID%Name) == BarotropicVelocityU_) then
                            
                            CHUNK = CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
                            !$OMP PARALLEL PRIVATE(i,j)
                            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                            do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                    
                                PropertyX%Field%R2D(i,j) = 0.

                                if (PointsToFill2D(i,j-1) == OpenPoint .and.          &
                                    PointsToFill2D(i,j  ) == OpenPoint) then
                                    PropertyX%Field%R2D(i,j) = (Matrix2D(i,j-1)+Matrix2D(i,j))/2.
                                endif

                            enddo
                            enddo
                            !$OMP END DO
                            !$OMP END PARALLEL

                            deallocate(Matrix2D)

                        endif

                        if (GetPropertyIDNumber(PropertyX%ID%Name) == BarotropicVelocityV_) then

                            CHUNK = CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
                            !$OMP PARALLEL PRIVATE(i,j)
                            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)                        
                            do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                    
                                PropertyX%Field%R2D(i,j) = 0.

                                if (PointsToFill2D(i-1,j) == OpenPoint .and.          &
                                    PointsToFill2D(i  ,j) == OpenPoint) then
                                    PropertyX%Field%R2D(i,j) = (Matrix2D(i-1,j)+Matrix2D(i,j))/2.
                                endif

                            enddo
                            enddo
                            !$OMP END DO
                            !$OMP END PARALLEL

                            deallocate(Matrix2D)

                        endif

                    endif

                    nullify   (PointsToFill2D)
                    nullify   (Matrix2D)
                
                else if (PropertyX%Dim == Dim_3D) then

                    if (PropertyX%Field%TypeZUV == TypeZ_) then

                        PointsToFill3D => WaterPoints3D

                        if (GetPropertyIDNumber(PropertyX%ID%Name) == VelocityU_ .or. &
                            GetPropertyIDNumber(PropertyX%ID%Name) == VelocityV_) then

                            allocate (Matrix3D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
                            Matrix3D(:,:,:) = FillValueReal

                        else

                            Matrix3D => PropertyX%Field%R3D

                        endif

                    else if (PropertyX%Field%TypeZUV == TypeU_) then

                        PointsToFill3D => WaterFaces3D_U
                        Matrix3D       => PropertyX%Field%R3D

                    else if (PropertyX%Field%TypeZUV == TypeV_) then

                        PointsToFill3D => WaterFaces3D_V
                        Matrix3D       => PropertyX%Field%R3D

                    else

                        stop 'ReadInitialImposedSolution  - ModuleHydrodynamic - ERR230'

                    endif


                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,         &
                                           Matrix3D       = Matrix3D,                &
                                           PointsToFill3D = PointsToFill3D,                      &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

                    if (PropertyX%Field%TypeZUV == TypeZ_) then
                        if (GetPropertyIDNumber(PropertyX%ID%Name) == VelocityU_) then

                            CHUNK = CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
                            !$OMP PARALLEL PRIVATE(i,j,k)
                            do k = Me%WorkSize%KLB,Me%WorkSize%KUB
                            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                            do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                    
                                PropertyX%Field%R3D(i,j,k) = 0.

                                if (PointsToFill3D(i,j-1,k) == OpenPoint .and.          &
                                    PointsToFill3D(i,j  ,k) == OpenPoint) then
                                    PropertyX%Field%R3D(i,j,k) = (Matrix3D(i,j-1,k)+Matrix3D(i,j,k))/2.
                                endif

                            enddo
                            enddo
                            !$OMP END DO
                            enddo
                            !$OMP END PARALLEL

                            deallocate(Matrix3D)

                        endif

                        if (GetPropertyIDNumber(PropertyX%ID%Name) == VelocityV_) then

                            CHUNK = CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
                            !$OMP PARALLEL PRIVATE(i,j,k)
                            do k = Me%WorkSize%KLB,Me%WorkSize%KUB
                            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                            do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                    
                                PropertyX%Field%R3D(i,j,k) = 0.

                                if (PointsToFill3D(i-1,j,k) == OpenPoint .and.          &
                                    PointsToFill3D(i  ,j,k) == OpenPoint) then
                                    PropertyX%Field%R3D(i,j,k) = (Matrix3D(i-1,j,k)+Matrix3D(i,j,k))/2.
                                endif

                            enddo
                            enddo
                            !$OMP END DO
                            enddo
                            !$OMP END PARALLEL

                            deallocate(Matrix3D)

                        endif
                    endif

                    nullify   (PointsToFill3D)
                    nullify   (Matrix3D)

                endif

            endif
            
            PropertyX%LastActualization = Me%ActualTime


        endif cd1

        if (MonitorPerformance) then
            call StopWatch ("ModuleAssimilation", "AssimilationFromFile")
        endif

        !Do if necessary the output of the Assimilation propertyX
        if (Me%OutPut%ON) call OutPutResultsHDF   
        
        !Output TimeSerie
        call OutPut_TimeSeries  (PropertyX)


        call UngetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

        call UngetHorizontalMap(Me%ObjHorizontalMap, WaterFaces2D_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

        call UngetHorizontalMap(Me%ObjHorizontalMap, WaterFaces2D_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

        call UngetMap(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

        call UngetMap(Me%ObjMap, WaterFaces3D_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

        call UngetMap(Me%ObjMap, WaterFaces3D_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AssimilationFromFile - ModuleAssimilation - ERR01'

    end subroutine AssimilationFromFile


    !--------------------------------------------------------------------------

    
    subroutine OutPutResultsHDF

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !External--------------------------------------------------------------
        real,       dimension(:,:,:), pointer       :: SZZ
        integer,    dimension(:,:,:), pointer       :: OpenPoints3D
        integer                                     :: STAT_CALL
        integer,    dimension(6)                    :: TimeAux
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second
                                                    
        !Local-----------------------------------------------------------------
        type(T_Time)                                :: Actual, EndTime
        integer                                     :: OutPutNumber
        integer                                     :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        real,       dimension(6),     target        :: AuxTime
        real,       dimension(:),     pointer       :: TimePtr

        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 
     
        Actual   = Me%ActualTime
        EndTime  = Me%EndTime

        !SZZ
        call GetGeometryDistances (Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR01'

        !OpenPoints3D
        call GetOpenPoints3D(Me%ObjMap, OpenPoints3D = OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR02'
          
        OutPutNumber = Me%OutPut%Next

        if (OutPutNumber <= Me%OutPut%Number) then
TOut:   if (Actual >= Me%OutPut%OutTime(OutPutNumber)) then 

            Me%OutPut%Next = Me%OutPut%Next + 1
            
            call ExtractDate(Actual, Year = Year, Month  = Month,  Day    = Day, &
                                     Hour = Hour, Minute = Minute, Second = Second)

            TimeAux(1) = int(Year  )
            TimeAux(2) = int(Month )
            TimeAux(3) = int(Day   )
            TimeAux(4) = int(Hour  )
            TimeAux(5) = int(Minute)
            TimeAux(6) = int(Second)

            !Writes current time
            call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),                 &
                                AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR03'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                                 Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR04'


            !Writes SZZ
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                                 WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR05'

            call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",                 &
                                 "m", Array3D = SZZ,                                        &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR06'

            !Writes OpenPoints
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB, WorkKLB, WorkKUB,                        &
                                 STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR07'

            call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",              &
                                 "-", Array3D = OpenPoints3D, OutputNumber = OutPutNumber,  &
                                 STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR08'

            PropertyX => Me%FirstAssimilationProp
  
            do while (associated(PropertyX)) 

                if (PropertyX%OutputHDF) then
                if (PropertyX%Dim == Dim_2D) then

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/"//PropertyX%ID%Name, PropertyX%ID%Name, &
                                         PropertyX%ID%Units, Array2D = PropertyX%Field%R2D,             &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR09'


                    call HDF5WriteData  (Me%ObjHDF5, "/DecayCoefs/"//PropertyX%ID%Name, PropertyX%ID%Name, &
                                         PropertyX%ID%Units, Array2D = PropertyX%CoefField%R2D,            &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR14'

                elseif (PropertyX%Dim == Dim_3D) then

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/"//PropertyX%ID%Name, PropertyX%ID%Name, &
                                         PropertyX%ID%Units, Array3D = PropertyX%Field%R3D,             &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR10'


                    call HDF5WriteData  (Me%ObjHDF5, "/DecayCoefs/"//PropertyX%ID%Name, PropertyX%ID%Name, &
                                         PropertyX%ID%Units, Array3D = PropertyX%CoefField%R3D,            &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR10'

                else 

                    stop 'OutPutResultsHDF - ModuleAssimilation - ERR11'

                endif                
                endif



                PropertyX => PropertyX%Next

            enddo

            nullify (PropertyX)

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR11'
        
        endif  TOut
        endif

        !SZZ
        call UnGetGeometry (Me%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR12'

        !OpenPoints3D
        call UnGetMap(Me%ObjMap, OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OutPutResultsHDF - ModuleAssimilation - ERR13'
      
            
       !----------------------------------------------------------------------

    end subroutine OutPutResultsHDF

    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries(PropertyX)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

cd1:    if (Me%ActualTime > PropertyX%LastTimeSerieOutput) then
            
cd2:        if (PropertyX%TimeSerie) then

cd3:            if (PropertyX%Dim == Dim_2D) then 

                    call WriteTimeSerie(Me%ObjTimeSerie,                        &
                                        Data2D  = PropertyX%Field%R2D,          &
                                        STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'OutPut_TimeSeries - ModuleAssimilation - ERR01'

                else if (PropertyX%Dim == Dim_3D) then cd3

                    call WriteTimeSerie(Me%ObjTimeSerie,                        &
                                        Data3D  = PropertyX%Field%R3D,          &
                                        STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'OutPut_TimeSeries - ModuleAssimilation - ERR02'


                endif cd3

            endif cd2

            PropertyX%LastTimeSerieOutput = Me%ActualTime

        endif cd1
    
    end subroutine OutPut_TimeSeries

    !-------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------

    subroutine ReadLockExternalVar              ! J. Nogueira Assim
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Now
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAssimilation - ERR01'

        !OpenPoints3D
        call GetOpenPoints3D(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAssimilation - ERR02'

        !SZZ
        call GetGeometryDistances (Me%ObjGeometry, SZZ = Me%ExternalVar%SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAssimilation - ERR03'

        !DWZ
        call GetGeometryDistances (Me%ObjGeometry, DWZ = Me%ExternalVar%DWZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAssimilation - ERR04'
        
        !ZCellCenter
        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  ZCellCenter   = Me%ExternalVar%ZCellCenter,           &
                                  STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAssimilation - ERR05'

        !KFloor_Z
        call GetGeometryKFloor(Me%ObjGeometry, Z = Me%ExternalVar%KFloor_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAssimilation - ERR06'

        !GridCellArea
        call GetGridCellArea (Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea , STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAssimilation - ERR07'
    
    end subroutine ReadLockExternalVar


    !--------------------------------------------------------------------------


    subroutine ReadUnlockExternalVar            ! J. Nogueira Assim
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------


        call UnGetMap(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAssimilation - ERR01'
        
        call UnGetGeometry(Me%ObjGeometry,Me%ExternalVar%DWZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAssimilation - ERR02'
        
        call UnGetGeometry(Me%ObjGeometry,Me%ExternalVar%SZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAssimilation - ERR03'

        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%ZCellCenter, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAssimilation - ERR04'

        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%KFloor_Z, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAssimilation - ERR05'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAssimilation - ERR06'
                               
        call null_time(Me%ExternalVar%Now)

    end subroutine ReadUnlockExternalVar

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillAssimilation(AssimilationID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: AssimilationID             
        integer, optional, intent(OUT)          :: STAT

        !External--------------------------------------------------------------
        integer                                 :: ready_             
        integer                                 :: STAT_CALL, nUsers

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_


        call Ready(AssimilationID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mASSIMILATION_,  Me%InstanceID)

            if (nUsers == 0) then

                call WriteFinalHDFOutPut
                
                ! guillaume : aqui chama-se a rotina de dealocacao de memoria de variaveis
                if(Me%AltimetricAssimilation)                                &
                        call KillAltimetricAssimilation

                if (Me%ObjTimeSerie /= 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillAssimilation - ModuleAssimilation - ERR01'
                endif

                call DeallocateVariables

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillAssimilation - ModuleAssimilation - ERR02'

                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjGridData)
                if (nUsers == 0) stop 'KillAssimilation - ModuleAssimilation - ERR03'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillAssimilation - ModuleAssimilation - ERR04'

                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillAssimilation - ModuleAssimilation - ERR05'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjGeometry)
                if (nUsers == 0) stop 'KillAssimilation - ModuleAssimilation - ERR06'

                nUsers = DeassociateInstance(mMAP_,             Me%ObjMap)
                if (nUsers == 0) stop 'KillAssimilation - ModuleAssimilation - ERR07'

                !Deallocates Instance
                call DeallocateInstance ()

                AssimilationID  = 0
                STAT_           = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           
        !----------------------------------------------------------------------

    end subroutine KillAssimilation

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! nogueira e guillaume : dealocacao de variaveis alocadas para Cooper-Haines
    subroutine KillAltimetricAssimilation

 
       !Deallocation of internal variables--------------------------------------------

        Deallocate (Me%Altimetric_Assim%ColumnTemperature,                                               &
                    Me%Altimetric_Assim%ColumnSalinity,                                                  &
                    Me%Altimetric_Assim%ColumnDensity,                                                   &
                    Me%Altimetric_Assim%DeltaTemperature,                                                &
                    Me%Altimetric_Assim%DeltaSalinity,                                                   &
                    Me%Altimetric_Assim%DeltaDensity,                                                    &
                    Me%Altimetric_Assim%DeltaPressure)

        Deallocate (Me%Altimetric_Assim%Delta_WaterLevel,                                       &
                    Me%Altimetric_Assim%Delta_WaterLevelToAssimilate,                           &
                    Me%Altimetric_Assim%PressureAnomaly)

        Deallocate (Me%Altimetric_Assim%Observation_Error,                                      &
                    Me%Altimetric_Assim%Model_Error,                                            &
                    Me%Altimetric_Assim%Gain)
          

        Deallocate (Me%Altimetric_Assim%AuxDepth,                                                &
                    Me%Altimetric_Assim%D_GradTemp,                                              &
                    Me%Altimetric_Assim%AuxSpline,                                               &
                    Me%Altimetric_Assim%D_GradSalin)

        Deallocate (Me%Altimetric_Assim%SigmaDensAnalyzed)

        !----------------------------------------------------------------------

    end subroutine KillAltimetricAssimilation

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Assimilation), pointer          :: AuxObjAssimilation
        type (T_Assimilation), pointer          :: PreviousObjAssimilation

        !Updates pointers
        if (Me%InstanceID == FirstObjAssimilation%InstanceID) then
            FirstObjAssimilation    => FirstObjAssimilation%Next
        else
            PreviousObjAssimilation => FirstObjAssimilation
            AuxObjAssimilation      => FirstObjAssimilation%Next
            do while (AuxObjAssimilation%InstanceID /= Me%InstanceID)
                PreviousObjAssimilation => AuxObjAssimilation
                AuxObjAssimilation      => AuxObjAssimilation%Next
            enddo

            !Now update linked list
            PreviousObjAssimilation%Next => AuxObjAssimilation%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance


    !--------------------------------------------------------------------------

    subroutine DeallocateVariables

        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------

        ! Deallocates all properties 
        PropertyX => Me%FirstAssimilationProp

do1 :   do while(associated(PropertyX))  

cd1:        if (PropertyX%Dim == Dim_2D) then        

                deallocate(PropertyX%Field%R2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    call SetError(FATAL_, INTERNAL_, 'DeAllocateVariables - ModuleAssimilation - ERR01') 

                nullify   (PropertyX%Field%R2D)

                deallocate(PropertyX%CoefField%R2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    call SetError(FATAL_, INTERNAL_, 'DeAllocateVariables - ModuleAssimilation - ERR02') 

                nullify   (PropertyX%CoefField%R2D)


            else if (PropertyX%Dim == Dim_3D) then  cd1                      

                deallocate(PropertyX%Field%R3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    call SetError(FATAL_, INTERNAL_, 'DeAllocateVariables - ModuleAssimilation - ERR03') 

                nullify   (PropertyX%Field%R3D)


                deallocate(PropertyX%CoefField%R3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    call SetError(FATAL_, INTERNAL_, 'DeAllocateVariables - ModuleAssimilation - ERR04') 

                nullify   (PropertyX%CoefField%R3D)


            endif cd1

            if(PropertyX%ID%SolutionFromFile)then

                call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    call SetError(FATAL_, INTERNAL_, 'DeAllocateVariables - ModuleAssimilation - ERR05') 

            end if

            if(PropertyX%CoefID%SolutionFromFile)then

                call KillFillMatrix(PropertyX%CoefID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    call SetError(FATAL_, INTERNAL_, 'DeAllocateVariables - ModuleAssimilation - ERR06') 

            end if


            PropertyX => PropertyX%Next

        end do do1

        !Sets the number of properties equal to the FillValueInt
        Me%PropertiesNumber = FillValueInt

        Nullify   (Me%FirstAssimilationProp,Me%LastAssimilationProp)

        !----------------------------------------------------------------------

    end subroutine DeallocateVariables   

    !--------------------------------------------------------------------------

    subroutine WriteFinalHDFOutPut

        !Local----------------------------------------------------------------

        !----------------------------------------------------------------------

        if (Me%OutPut%ON) then
        
            call OutPutResultsHDF


            ! Close the output transient HDF output file
            call KillHDF5(Me%ObjHDF5)

        endif


    end subroutine WriteFinalHDFOutPut


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Ready (ObjAssimilation_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjAssimilation_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjAssimilation_ID > 0) then
            call LocateObjAssimilation (ObjAssimilation_ID)
            ready_ = VerifyReadLock (mASSIMILATION_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjAssimilation (ObjAssimilationID)

        !Arguments-------------------------------------------------------------
        integer                    :: ObjAssimilationID

        !Local-----------------------------------------------------------------

        Me => FirstObjAssimilation
        do while (associated (Me))
            if (Me%InstanceID == ObjAssimilationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleAssimilation - LocateObjAssimilation - ERR01'

    end subroutine LocateObjAssimilation

    !--------------------------------------------------------------------------

end module ModuleAssimilation

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

