!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Hydrodynamic File 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module which writes/reads the hydrodynamic solution to/from a binary file
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
! Module Hydrodynamic File
!
!   - Output File  structure
!
!       - Header (written by WriteHeader)
!           - Time limits of the run
!           - DT of the output
!           - Initial (Instant before first output):
!               - WaterLevel
!               - Velocity
!               - ComputeFaces
!           
!
!       - n Blocks (written by WriteBlock)
!           - Current Time
!           - WaterLevel
!           - Horizontal WaterFluxes
!           - Discharges
!           - HorizontalComputeFaces
!
!
!   - Input File type
!
!           - Begin - End
!               - all the Hydrodynamic is stored in the file
!               - One call to ReadHeader
!               - for every instant one call to read ReadsBlock
!
!           - Mohid3D M2
!               - one cycle of M2 stored in the file (12h25min30sec)
!               - searches for the right instant and reads respectiv block
!
Module ModuleHydrodynamicFile

    use ModuleGlobalData
    use ModuleTime               
    use ModuleEnterData
    use ModuleGridData,         only : GetGridData, WriteGridData, UnGetGridData
    use ModuleHorizontalMap
    use ModuleHorizontalGrid     
    use ModuleGeometry,         only : GetGeometrySize, GetGeometryVolumes, UnGetGeometry
    use ModuleMap,              only : GetOpenPoints3D, GetWaterPoints3D, UngetMap
    use ModuleHydroIntegration      

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartHydrodynamicFile      
    private ::      HydrodynamicFileOptions
    private ::      VerifyHydrodynamicFileOptions
    private ::      StartHydrodynamicFileInput
    private ::          AllocateVariablesInput
    private ::      StartHydrodynamicFileOutput
    private ::          AllocateVariablesOutput
    private ::          IntegrateBathymetry


    !Selector
    public  :: GetHydrodynamicFileIOState
    public  :: GetFileWaterLevel
    public  :: GetFileFluxes
    public  :: GetFileMapping

    public  :: UngetHydrodynamicFile
                     
    
    !Modifier
    public  :: ModifyHydrodynamicFile            
    private ::      ModifyMohid3d_M2
    private ::          M2_Iteration                !Function
    private ::      IntegrateInSpace

    private :: ReadHeader                           !Input
    private :: ReadsBlock                           !Input
            
    private :: WriteHeader                          !Output
    private :: WriteBlock                           !Output    

    !Destructor
    public  ::  KillHydrodynamicFile                                                     
    private ::      DeallocateVariablesInput

 
    !Management
    private ::      Ready
    private ::          LocateObjHydrodynamicFile


    !Interfaces----------------------------------------------------------------

    private :: UngetHydrodynamicFile2Dreal4
    private :: UngetHydrodynamicFile2Dreal8
    private :: UngetHydrodynamicFile3Dreal4
    private :: UngetHydrodynamicFile3Dreal8
    private :: UnGetHydrodynamicFile3Dinteger
    interface  UngetHydrodynamicFile
        module procedure UngetHydrodynamicFile2Dreal4
        module procedure UngetHydrodynamicFile2Dreal8
        module procedure UngetHydrodynamicFile3Dreal4
        module procedure UngetHydrodynamicFile3Dreal8
        module procedure UngetHydrodynamicFile3Dinteger
    end interface UngetHydrodynamicFile


    private :: ModifyHydrodynamicFileInput
    private :: ModifyHydrodynamicFileOutput
    interface  ModifyHydrodynamicFile
        module procedure ModifyHydrodynamicFileInput
        module procedure ModifyHydrodynamicFileOutput
    end interface ModifyHydrodynamicFile

    !Parameter-----------------------------------------------------------------

    !Hydrodynamic
    real,    parameter :: M2_Period_Hours = 12.425 ! hours of the M2 tide wave 
    
    ! File types 
    integer, parameter :: M2_Tide_type    = 1
    integer, parameter :: BeginEnd_type   = 2

    ! InicJuan
    ! Bathymetry integration type
    integer, parameter :: MaxVal_Type  = 1
    integer, parameter :: MeanVal_Type = 2
    
    character(LEN = StringLength), parameter :: Char_MaxVal_Type  = trim(adjustl('MaxVal_Type' ))
    character(LEN = StringLength), parameter :: Char_MeanVal_Type = trim(adjustl('MeanVal_Type'))
    ! FimJuan

    ! File type
    character(LEN = StringLength), parameter :: Char_M2_Tide_type    = trim(adjustl('M2_Tide_type'   ))
    character(LEN = StringLength), parameter :: Char_BeginEnd_type   = trim(adjustl('BeginEnd_type'  ))

    !Input/Output State
    integer, parameter :: InputState  = 1
    integer, parameter :: OutputState = 2

    !Output control
    integer, parameter :: BeginBlock = 456
    integer, parameter :: EndBlock   = 654

    !Types---------------------------------------------------------------------

    type       T_State
        logical :: INPUT      = OFF
        logical :: OUTPUT     = OFF
        logical :: TimeIntegration = OFF
        logical :: SpaceIntegration = OFF
    end type T_State

    type       T_Files
         character(len=StringLength) :: InPutFields             = trim(adjustl('********************'))
         integer                     :: InPutFieldsUnit         = null_int
         integer                     :: InputVersionNumber      = 2
         character(len=StringLength) :: OutPutFields            = trim(adjustl('********************'))
         integer                     :: OutPutFieldsUnit        = null_int
         integer                     :: OutputVersionNumber     = 2
         character(len=StringLength) :: IntegratedBathymetry    = trim(adjustl('********************'))
    end type T_Files

    type       T_External
        !ObjTime
        type(T_Time)                       :: Now
        real                               :: DT = null_real
    end type T_External

    type T_BlockInMemory
        real,    pointer, dimension(:)              :: WaterLevel
        real(8), pointer, dimension(:)              :: WaterFluxX
        real(8), pointer, dimension(:)              :: WaterFluxY
        real(8), pointer, dimension(:)              :: Discharges
        integer, pointer, dimension(:)              :: ComputeFacesU3D
        integer, pointer, dimension(:)              :: ComputeFacesV3D
    end type T_BlockInMemory


    !Input Data
    type T_Input
        type(T_Time)                        :: StartTime
        type(T_Time)                        :: EndTime
        type(T_Time)                        :: NextInput
        integer                             :: File_Type
        logical                             :: LoadAllToMemory
        integer                             :: NextBlockInMemory
        real                                :: DT_HYDROFILE
        !Initial
        real,    pointer, dimension(:,:  )  :: InitialWaterLevel
        real(8), pointer, dimension(:,:,:)  :: InitialWaterFluxX
        real(8), pointer, dimension(:,:,:)  :: InitialWaterFluxY
        real(8), pointer, dimension(:,:,:)  :: InitialDischarges
        integer, pointer, dimension(:,:,:)  :: InitialComputeFacesU3D
        integer, pointer, dimension(:,:,:)  :: InitialComputeFacesV3D
        integer, pointer, dimension(:,:  )  :: WaterPoints2D

        !Transient
        real,    pointer, dimension(:,:  )  :: WaterLevel_New
        real(8), pointer, dimension(:,:,:)  :: WaterFluxX
        real(8), pointer, dimension(:,:,:)  :: WaterFluxY
        real(8), pointer, dimension(:,:,:)  :: Discharges
        integer, pointer, dimension(:,:,:)  :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:)  :: ComputeFacesV3D
        type(T_Time)                        :: TransientTime

        !AllInMemory
        type (T_BlockInMemory), dimension(:), pointer :: BlockInMemory
    end type T_Input

       
    !Output information from 3D type
    type       T_STIntegration
        real,    pointer, dimension(:,:  ) :: WaterLevel
        real(8), pointer, dimension(:,:,:) :: WaterFluxX
        real(8), pointer, dimension(:,:,:) :: WaterFluxY
        real(8), pointer, dimension(:,:,:) :: Discharges
        integer, dimension(:,:,:), pointer :: ComputeFacesU3D
        integer, dimension(:,:,:), pointer :: ComputeFacesV3D
        integer, dimension(:,:  ), pointer :: WaterPoints2D
        logical                            :: WindowOn 
        type(T_Size3D)                     :: WindowSize
        type(T_Size3D)                     :: WindowWorkSize
        real                               :: IntegrationStep = null_real
    end type T_STIntegration

    type       T_TimeIntegration
        type(T_STIntegration) :: Integrate3DTime
        type(T_Time)          :: NextOutput
    end type T_TimeIntegration

    type       T_SpaceIntegration
        type(T_STIntegration)            :: Integrate3DSpace
        integer                          :: BatIntegrationType
        integer, pointer, dimension(:)   :: VectorIntegrationX
        integer, pointer, dimension(:)   :: VectorIntegrationY
        real,    pointer, dimension(:,:) :: BatMin
    end type T_SpaceIntegration


    type      T_HydrodynamicFile
        integer                                     :: InstanceID   
        type(T_Size3D  )                            :: Size   
        type(T_Size3D  )                            :: WorkSize   
        type(T_State   )                            :: State
        type(T_Files   )                            :: Files
        type(T_External)                            :: ExternalVar
        type(T_Input    )                           :: Input
        type(T_TimeIntegration)                     :: TimeIntegration
        type(T_SpaceIntegration)                    :: SpaceIntegration

        !Instance of other modules
        integer                                     :: ObjTopography        = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjGeometry          = 0
        integer                                     :: ObjMap               = 0
        integer                                     :: ObjHorizontalMap     = 0
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjHydroIntegration  = 0
        type (T_HydrodynamicFile), pointer          :: Next

    end type T_HydrodynamicFile

    !Global Module Variables
    type (T_HydrodynamicFile), pointer              :: FirstHydrodynamicFile
    type (T_HydrodynamicFile), pointer              :: Me

    !--------------------------------------------------------------------------
   
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartHydrodynamicFile(HydrodynamicFileID,                        &
                                     TopographyID,                              &
                                     HorizontalGridID,                          &
                                     GeometryID,                                &
                                     MapID,                                     & 
                                     HorizontalMapID,                           & 
                                     TimeID,                                    &
                                     InPutOutPutState,                          &
                                     InitialWaterLevel,                         &
                                     InitialWaterFluxX,                         &
                                     InitialWaterFluxY,                         &
                                     InitialDischarges,                         &
                                     InitialComputeFacesU3D,                    &
                                     InitialComputeFacesV3D,                    &
                                     STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        integer                                     :: TopographyID
        integer                                     :: HorizontalGridID
        integer                                     :: GeometryID
        integer                                     :: MapID                              
        integer                                     :: HorizontalMapID                              
        integer                                     :: TimeID                              
        integer                                     :: InPutOutPutState
        real,    pointer, dimension(:,:  )          :: InitialWaterLevel
        real(8), pointer, dimension(:,:,:)          :: InitialWaterFluxX
        real(8), pointer, dimension(:,:,:)          :: InitialWaterFluxY
        real(8), pointer, dimension(:,:,:)          :: InitialDischarges
        integer, pointer, dimension(:,:,:)          :: InitialComputeFacesU3D
        integer, pointer, dimension(:,:,:)          :: InitialComputeFacesV3D
        integer, optional, intent(OUT)              :: STAT     

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        character(PathLength)                       :: DataFile

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mHydrodynamicFile_)) then
            nullify (FirstHydrodynamicFile)
            call RegisterModule (mHydrodynamicFile_) 
        endif

        call Ready(HydrodynamicFileID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance

            !Associates External Instances
            Me%ObjTopography     = AssociateInstance (mGRIDDATA_,       TopographyID    )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMAP_,            MapID           )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )


            !Current Time
            call GetComputeCurrentTime  (Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFile - ModuleHydrodynamicFile - ERR01'

            !Compute Step
            call GetComputeTimeStep     (Me%ObjTime, Me%ExternalVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFile - ModuleHydrodynamicFile - ERR02'

            call GetGeometrySize(Me%ObjGeometry,                        &
                                 Size     = Me%Size,                    &
                                 WorkSize = Me%WorkSize,                &
                                 STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFile - ModuleHydrodynamicFile - ERR03' 


            call ReadFileName('IN_HYDRO_FILE', DataFile,                                 &
                               Message = 'HydrodynamicFile Data File', STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFile - ModuleHydrodynamicFile - ERR04' 

            call ConstructEnterData(Me%ObjEnterData, DataFile, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFile - ModuleHydrodynamicFile - ERR03' 
              
            !Sets the Input Output State
            if     (InPutOutPutState == InputState) then
                Me%State%Input  = .true.
                Me%State%Output = .false.
            elseif (InPutOutPutState == OutputState) then
                Me%State%Input  = .false.
                Me%State%Output = .true.
            else
                stop 'StartHydrodynamicFile - ModuleHydrodynamicFile - ERR05'
            end if

            !Gets the options from the file
            call HydrodynamicFileOptions

            !Starts Input
            if (Me%State%INPUT) then

                call StartHydrodynamicFileInput

                InitialWaterLevel     (:,:  ) =  Me%Input%InitialWaterLevel     (:,:  )    
                InitialWaterFluxX     (:,:,:) =  Me%Input%InitialWaterFluxX     (:,:,:)
                InitialWaterFluxY     (:,:,:) =  Me%Input%InitialWaterFluxY     (:,:,:)
                InitialDischarges     (:,:,:) =  Me%Input%InitialDischarges     (:,:,:)
                !This module can not change the mapping information
                InitialComputeFacesU3D        => Me%Input%InitialComputeFacesU3D 
                InitialComputeFacesV3D        => Me%Input%InitialComputeFacesV3D 

            endif

            !Starts Output
            if (Me%State%OUTPUT) then

                Me%Input%InitialWaterLevel       => InitialWaterLevel        
                Me%Input%InitialWaterFluxX       => InitialWaterFluxX     
                Me%Input%InitialWaterFluxY       => InitialWaterFluxY     
                Me%Input%InitialDischarges       => InitialDischarges     
                Me%Input%InitialComputeFacesU3D  => InitialComputeFacesU3D
                Me%Input%InitialComputeFacesV3D  => InitialComputeFacesV3D

                call StartHydrodynamicFileOutput

                nullify(Me%Input%InitialWaterLevel     )
                nullify(Me%Input%InitialWaterFluxX     )
                nullify(Me%Input%InitialWaterFluxY     )
                nullify(Me%Input%InitialDischarges     )
                nullify(Me%Input%InitialComputeFacesU3D)
                nullify(Me%Input%InitialComputeFacesV3D)

            endif

            !Verifies numerical options
            call VerifyHydrodynamicFileOptions

            !Kills EnterData
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartHydrodynamicFile - ModuleHydrodynamicFile - ERR07'

            STAT_ = SUCCESS_

            HydrodynamicFileID = Me%InstanceID

        else cd0
            

            stop 'StartHydrodynamicFile - HydrodynamicFile - ERR99' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine StartHydrodynamicFile

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        type (T_HydrodynamicFile), pointer          :: NewObjHydrodynamicFile
        type (T_HydrodynamicFile), pointer          :: PrevObjHydrodynamicFile


        !Allocates new instance
        allocate (NewObjHydrodynamicFile)
        nullify  (NewObjHydrodynamicFile%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstHydrodynamicFile)) then
            FirstHydrodynamicFile       => NewObjHydrodynamicFile
            Me                          => NewObjHydrodynamicFile
        else
            PrevObjHydrodynamicFile     => FirstHydrodynamicFile
            Me                          => FirstHydrodynamicFile%Next
            do while (associated(Me))
                PrevObjHydrodynamicFile  => Me
                Me                      => Me%Next
            enddo
            Me                          => NewObjHydrodynamicFile
            PrevObjHydrodynamicFile%Next=> NewObjHydrodynamicFile
        endif

        Me%InstanceID = RegisterNewInstance (mHYDRODYNAMICFILE_)

    end subroutine AllocateInstance

    !----------------------------------------------------------------------------

    subroutine HydrodynamicFileOptions

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: FromFile
        integer                                     :: flag
        integer, dimension(4)                       :: aux
        character(LEN = StringLength)               :: String

        !------------------------------------------------------------------------

        call GetExtractType    (FromFile = FromFile)

        !Input Variables
ifin:   if (Me%State%Input) then
            
            !Input file
            call GetData(Me%Files%InPutFields,                                           &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='IN_FIELD',                                       &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR01'

            !Input File type
            call GetData(String,                                                         &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='IN_FILE_TYPE',                                   &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         default      = Char_M2_Tide_type,                               &
                         STAT         = STAT_CALL)             
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR02'

            String = trim(adjustl(String))
            select case(String)
                case(Char_M2_Tide_type)

                    Me%Input%File_Type = M2_Tide_type

                case(Char_BeginEnd_type)

                    Me%Input%File_Type = BeginEnd_type

                case default

                    stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR03'
    
            end select

            !Load all to memory
            call GetData(Me%Input%LoadAllToMemory,                                       &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='LOAD_TO_MEMORY',                                 &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         default      = .false.,                                         &
                         STAT         = STAT_CALL)             
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR04'

            !Input File Version
            call GetData(Me%Files%InputVersionNumber,                                    &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='IN_FILE_VERSION',                                &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         default      = 2,                                               &
                         STAT         = STAT_CALL)             
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR05'


        endif ifin

        !Output Variables
ifout:  if (Me%State%Output) then
            
            !OutputFile
            call GetData(Me%Files%OutPutFields,                                          &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='OUT_FIELD',                                      &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR06'

            !<BeginKeyword>
               !Keyword          : TIME_INTEGRATION
               !<BeginDescription>       
               !Verifies if the integration of fluxes in time is to be done
               !<EndDescription>
               !Type             : Logical
               !Default          : False
               !File keyword     : IN_HYDRO_FILE
               !Multiple Options : do not have
               !Search Type      : From File
            !<EndKeyword>
            
            !Time Integration on
            call GetData(Me%State%TimeIntegration,                                       &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='TIME_INTEGRATION',                               &
                         Default      = OFF,                                             &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR07'
            
            !<BeginKeyword>
               !Keyword          : OUT_FILE_VERSION
               !<BeginDescription>       
               !Controls the version of the output file
               !<EndDescription>
               !Type             : Integer
               !Default          : 2
               !File keyword     : IN_HYDRO_FILE
               !Multiple Options : do not have
               !Search Type      : From File
            !<EndKeyword>

            !Output File Version
            call GetData(Me%Files%OutputVersionNumber,                                   &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='OUT_FILE_VERSION',                               &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         default      = 2,                                               &
                         STAT         = STAT_CALL)             
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR08'

            

            !<BeginKeyword>
               !Keyword          : SPACE_INTEGRATION
               !<BeginDescription>       
               !Verifies if the integration of fluxes in space is to be done
               !<EndDescription>
               !Type             : Logical
               !Default          : False
               !File keyword     : IN_HYDRO_FILE
               !Multiple Options : do not have
               !Search Type      : From File
            !<EndKeyword>

            call GetData(Me%State%SpaceIntegration,                                      &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='SPACE_INTEGRATION',                              &
                         Default      = OFF,                                             &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR09'
            
                       

            !DT to record file
            call GetData(Me%TimeIntegration%Integrate3DTime%IntegrationStep,             &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='DT_HYDROFILE',                                   &
                         ClientModule ='ModuleHydrodynamicFile',                         &
                         STAT         = STAT_CALL)             
            if (STAT_CALL /= SUCCESS_ .or. flag == 0)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR10'
            

            !<BeginKeyword>
               !Keyword          : N_ITEGRATION_CELLS
               !<BeginDescription>       
               !Number of cells that will be integrated (the integration space step)
               !<EndDescription>
               !Type             : Integer
               !Default          : False
               !File keyword     : IN_HYDRO_FILE
               !Multiple Options : do not have
               !Search Type      : From File
            !<EndKeyword>

            if (Me%State%SpaceIntegration) then

                call GetData(Me%SpaceIntegration%Integrate3DSpace%IntegrationStep,       &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromFile,                                    &
                             keyword      ='N_ITEGRATION_CELLS',                         &
                             ClientModule ='ModuleHydrodynamicFile',                     &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_ .or. flag == 0)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR11'

            !<BeginKeyword>
               !Keyword          : BAT_INTEGRATION_TYPE
               !<BeginDescription>       
               ! It is posible to calculate the new bathymetry using two options. 
               ! "MaxVal_Type"  -> Each new integrated cell has the maximum value of the cells 
               !                   used to do the integration of that cell.  
               ! "MeanVal_Type" -> The depth of the integrated cell is obtained by the average 
               !                   of the cells used to do the integration of that cell.
               !<EndDescription>
               !Type             : String
               !Default          : MaxVal_Type
               !File keyword     : IN_HYDRO_FILE
               !Multiple Options : MaxVal_Type, MeanVal_Type
               !Search Type      : From File
            !<EndKeyword>

                call GetData(String,                                                     &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromFile,                                    &
                             keyword      ='BAT_INTEGRATION_TYPE',                       &
                             ClientModule ='ModuleHydrodynamicFile',                     &
                             default      = Char_MeanVal_Type,                           &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR13'

                select case(trim(adjustl(String)))
                    case(Char_MeanVal_Type)

                        Me%SpaceIntegration%BatIntegrationType = MeanVal_Type

                    case(Char_MaxVal_Type)

                        Me%SpaceIntegration%BatIntegrationType = MaxVal_Type

                    case default

                        call SetError(FATAL_, INTERNAL_, "HydrodynamicFileOptions - ModuleHydrodynamicFile -ERR14")
    
                end select

            endif

            

            !Window to Record
            Me%TimeIntegration%Integrate3DTime%WindowOn = .false.
            call GetData(aux,                                                            &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='WINDOW',                                         &
                         ClientModule ='ModuleHydrodynamicFile')

            if (flag .NE. 0) then
                if (flag .EQ. 4) then
                    Me%TimeIntegration%Integrate3DTime%WindowOn   = .true.

                    Me%TimeIntegration%Integrate3DTime%WindowWorkSize%ILB = aux(1)
                    Me%TimeIntegration%Integrate3DTime%WindowWorkSize%IUB = aux(2)

                    Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JLB = aux(3)
                    Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JUB = aux(4)

                    Me%TimeIntegration%Integrate3DTime%WindowSize%ILB = aux(1) - 1
                    Me%TimeIntegration%Integrate3DTime%WindowSize%IUB = aux(2) + 1

                    Me%TimeIntegration%Integrate3DTime%WindowSize%JLB = aux(3) - 1
                    Me%TimeIntegration%Integrate3DTime%WindowSize%JUB = aux(4) + 1
                else
                    call SetError(FATAL_, INTERNAL_, "HydrodynamicFileOptions - ModuleHydrodynamicFile -ERR15")
                end if
            end if


            !Gets name of the new bathymetry
            if (Me%TimeIntegration%Integrate3DTime%WindowOn .or. Me%State%SpaceIntegration) then
        
                call GetData(Me%Files%IntegratedBathymetry,                              &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromFile,                                    &
                             keyword      ='NEW_BATIM',                                  &
                             ClientModule ='ModuleHydrodynamicFile',                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptions - ModuleHydrodynamicFile - ERR14'
            endif

        end if ifout

        !------------------------------------------------------------------------

    end subroutine HydrodynamicFileOptions         

    !--------------------------------------------------------------------------

    subroutine VerifyHydrodynamicFileOptions

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: exists
        integer                                     :: STAT_CALL
        type(T_Time)                                :: BeginTimeModel, EndTimeModel


        call GetComputeTimeStep(Me%ObjTime, DT = Me%ExternalVar%DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR00'


        !Verifies if input or output is on
        if ((.NOT. Me%State%INPUT) .AND. (.NOT. Me%State%OUTPUT))                        &
            stop 'VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR01'


        
        if (Me%State%INPUT) then

            !Input File exists
            inquire(FILE = Me%Files%InPutFields, EXIST = exists)
            if (.NOT. exists) then
                write(*,*)trim(Me%Files%InPutFields)//' does not exists'
                stop 'VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR02'
            end if

            !Time step
            if (Me%Input%DT_HYDROFILE /= Me%ExternalVar%DT) then
                write(*,*)
                write(*,*)'DT_HYDROFILE = ', Me%Input%DT_HYDROFILE
                write(*,*)'DT           = ', Me%ExternalVar%DT
                write(*,*)'DT_HYDROFILE and DT are different.' 
                stop 'VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR03'
            end if

            !Get Time Limits
            call GetComputeTimeLimits(Me%ObjTime,                                    &
                                      BeginTime = BeginTimeModel,                    &
                                      EndTime = EndTimeModel,                        &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR04'


            !Condition for Begin End Type
            !EndTime   (File) >= EndTime Model
            if (Me%Input%File_Type == BeginEnd_type) then
                if (EndTimeModel   - Me%Input%EndTime   >  0) then
                    write(*, *)'EndTime   (File) < EndTime Model'
                    stop 'VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR05'
                endif
            endif

            if (Me%Input%File_Type == M2_Tide_type) then

                !Verifies if Input (File) >= Start Time (Model)
                if (Me%Input%StartTime .gt. BeginTimeModel) then
                    write(*,*)'Start Time Model before StartTime HydrodynamicFile'
                    stop 'VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR06'
                endif

                !Verifies if StartTime consist with M2_iteration
                if (mod((BeginTimeModel - Me%Input%StartTime),          &
                         dble(Me%Input%DT_HYDROFILE)) /= 0.0) then
                    write(*,*)'StartTime inconsistent with M2_iteration'
                    stop 'VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR07'
                endif
            endif

        end if


        if (Me%State%OUTPUT) then
            if ((mod(Me%TimeIntegration%Integrate3DTime%IntegrationStep, Me%ExternalVar%DT) .EQ. 0.0) .AND.   &
                    (Me%TimeIntegration%Integrate3DTime%IntegrationStep .GE. Me%ExternalVar%DT)) then
                Me%TimeIntegration%NextOutput = Me%ExternalVar%Now + Me%TimeIntegration%Integrate3DTime%IntegrationStep

            else if (mod(Me%TimeIntegration%Integrate3DTime%IntegrationStep, Me%ExternalVar%DT) .NE. 0.0) then
                
                call SetError(ErrorMagnitude = FATAL_,                           &
                              ErrorType      = INTERNAL_,                        &
                              SmallMessage   ="mod(DT_HYDROFILE, DT) /= 0.") 

            else if (.NOT. (Me%TimeIntegration%Integrate3DTime%IntegrationStep .GE. Me%ExternalVar%DT)) then

                call SetError(ErrorMagnitude = FATAL_,                           &
                              ErrorType      = INTERNAL_,                        &
                              SmallMessage   ="DT_HYDROFILE < DT")

            else

                call SetError(FATAL_, INTERNAL_, "VerifyHydrodynamicFileOptions - ModuleHydrodynamicFile - ERR05") 

            end if 
        end if       

    end subroutine VerifyHydrodynamicFileOptions

    !--------------------------------------------------------------------------

    subroutine StartHydrodynamicFileInput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: exists, FILEOPEN, StartInputFound
        integer                                     :: STAT_CALL
        integer                                     :: M2_It, i, j
        integer, dimension(:, :), pointer           :: WaterPoints2D
        type(T_Time)                                :: BeginTimeModel


        !Allocates input variables    
        call AllocateVariablesInput


        !Verifies if the input file exists
        inquire(FILE = Me%Files%InPutFields, EXIST = exists)
        if (.NOT. exists) stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR01'

        !Verifies if the file is already open
        inquire(FILE = Me%Files%InPutFields, OPENED = FILEOPEN)
        if (.NOT. exists) stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR02'

        !Opens file
        call UnitsManager(Me%Files%InPutFieldsUnit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR03'

        open(UNIT   = Me%Files%InPutFieldsUnit,                                          &
             FILE   = Me%Files%InPutFields,                                              &
             FORM   = 'UNFORMATTED',                                                     &
             STATUS = 'OLD',                                                             &
             ACTION = 'READ',                                                            &
             IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR04'

        !Get Time Limits
        call GetComputeTimeLimits(Me%ObjTime,BeginTime = BeginTimeModel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR05'

        !Reads File Header
        call ReadHeader

BE:     if (Me%Input%File_Type == BeginEnd_type) then

            !Time step
            if (Me%Input%DT_HYDROFILE /= Me%ExternalVar%DT) then
                write(*,*)
                write(*,*)'DT_HYDROFILE = ', Me%Input%DT_HYDROFILE
                write(*,*)'DT           = ', Me%ExternalVar%DT
                write(*,*)'DT_HYDROFILE and DT are different.' 
                stop      'Subroutine StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR06'
            end if

            !Condition for Begin End Type
            !StartTime (File) <= StartTime Model
            if (BeginTimeModel < Me%Input%StartTime) then
                write(*, *)'StartTime (File) > StartTime Model'
                stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR07'
            endif

ST:         if (BeginTimeModel > Me%Input%StartTime) then

                StartInputFound = .false.
        
DW:             do while (.not. StartInputFound) 
                
                    call ReadsBlock

TT:                 if (Me%Input%TransientTime == BeginTimeModel) then

                        StartInputFound = .true.

                    else if (Me%Input%TransientTime > BeginTimeModel) then  TT

                        write(*, *)'The initial condition was not found'
                        stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR08'
                   
                    endif  TT
                
                enddo   DW


SF:             if (StartInputFound) then

                    !InitialWater Level is the last one readed
                    Me%Input%InitialWaterLevel(:,:)         = Me%Input%WaterLevel_New(:,:)
                    Me%Input%InitialWaterFluxX(:,:,:)       = Me%Input%WaterFluxX(:,:,:)
                    Me%Input%InitialWaterFluxY(:,:,:)       = Me%Input%WaterFluxY(:,:,:)
                    Me%Input%InitialDischarges(:,:,:)       = Me%Input%Discharges(:,:,:)
                    Me%Input%InitialComputeFacesU3D(:,:,:)  = Me%Input%ComputeFacesU3D(:,:,:)
                    Me%Input%InitialComputeFacesV3D(:,:,:)  = Me%Input%ComputeFacesV3D(:,:,:)

                endif  SF

            endif  ST

        endif  BE


        !If every thing is too load to memory, start encoding + allocate Variables
        if (Me%Input%LoadAllToMemory) then
            if (Me%Files%InputVersionNumber == 1) then
                write(*,*)'Cant encode Version 1 of Hydrodynamic File'
                stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR09'
            endif
            call StartEncoding 
            call UnitsManager  (Me%Files%InPutFieldsUnit, CLOSE_FILE, STAT = STAT_CALL)
        endif

        if (Me%Input%File_Type == M2_Tide_type .and.                    &
            (.not.Me%Input%LoadAllToMemory)) then

            !Reads the blocks until the initial time of the model
            if (Me%ExternalVar%Now .GT.                                 &
                Me%Input%StartTime) then
                M2_It = M2_Iteration(Me%ExternalVar%Now,                &
                                     Me%Input%StartTime,                &
                                     Me%Input%DT_HYDROFILE)
                do i = 1, M2_It
                    call ReadsBlock
                enddo

                !InitialWater Level is the last one readed
                Me%Input%InitialWaterLevel(:,:)         = Me%Input%WaterLevel_New(:,:)
                Me%Input%InitialWaterFluxX(:,:,:)       = Me%Input%WaterFluxX(:,:,:)
                Me%Input%InitialWaterFluxY(:,:,:)       = Me%Input%WaterFluxY(:,:,:)
                Me%Input%InitialDischarges(:,:,:)       = Me%Input%Discharges(:,:,:)
                Me%Input%InitialComputeFacesU3D(:,:,:)  = Me%Input%ComputeFacesU3D(:,:,:)
                Me%Input%InitialComputeFacesV3D(:,:,:)  = Me%Input%ComputeFacesV3D(:,:,:)

            endif

        endif

        if (Me%Files%InputVersionNumber > 1) then

            call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR10'

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (WaterPoints2D(i, j) /= Me%Input%WaterPoints2D(i, j)) then
                    write(*,*)"Input WaterPoint From file does not match Bathymetry"
                    write(*,*)"Point i = ", i
                    write(*,*)"Point j = ", j
                    write(*,*)'StartHydrodynamicFileInput - ModuleHydrodynamicFile - WRN05'
                endif
            enddo
            enddo

            call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicFileInput - ModuleHydrodynamicFile - ERR11'

        endif

        !Actualizes next output
        Me%Input%NextInput = Me%ExternalVar%Now + Me%Input%DT_HYDROFILE

    end subroutine StartHydrodynamicFileInput

    !----------------------------------------------------------------------------
    
    subroutine StartEncoding

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: NumberOfInPuts, iInput
        integer                                     :: NumberOfWaterPoints, nPoint
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k, Layers
        integer                                     :: STAT_CALL
        type (T_Time)                               :: BeginTime, EndTime
        

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB
        Layers = Me%Size%KUB - Me%Size%KLB + 1

        !NumberOfWaterPoints
        NumberOfWaterPoints = 0
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%Input%WaterPoints2D(i, j) == WaterPoint) then
                NumberOfWaterPoints = NumberOfWaterPoints + Layers
            endif
        enddo
        enddo

        if (Me%Input%File_Type == M2_Tide_type) then
            NumberOfInPuts = M2_Period_Hours*3600./Me%Input%DT_HYDROFILE + 1
        else

            !Gets Time Limits
            call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, &
                                      EndTime = EndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartEncoding - ModuleHydrodynamicFile - ERR01'
            NumberOfInputs = (EndTime-BeginTime)/Me%Input%DT_HYDROFILE + 1
        endif

        !Allocates Data
        allocate (Me%Input%BlockInMemory(NumberOfInPuts))

        do i = 1, NumberOfInputs
            allocate (Me%Input%BlockInMemory(i)%WaterLevel      (NumberOfWaterPoints))
            allocate (Me%Input%BlockInMemory(i)%WaterFluxX      (NumberOfWaterPoints))
            allocate (Me%Input%BlockInMemory(i)%WaterFluxY      (NumberOfWaterPoints))
            allocate (Me%Input%BlockInMemory(i)%Discharges      (NumberOfWaterPoints))
            allocate (Me%Input%BlockInMemory(i)%ComputeFacesU3D (NumberOfWaterPoints))
            allocate (Me%Input%BlockInMemory(i)%ComputeFacesV3D (NumberOfWaterPoints))
        enddo 

        !Encode Header
        nPoint = 0
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
        if (Me%Input%WaterPoints2D(i, j) == WaterPoint) then
            nPoint = nPoint + 1
            Me%Input%BlockInMemory(1)%WaterLevel      (nPoint) = Me%Input%InitialWaterLevel(i, j)
            Me%Input%BlockInMemory(1)%WaterFluxX      (nPoint) = Me%Input%InitialWaterFluxX(i, j, k)
            Me%Input%BlockInMemory(1)%WaterFluxY      (nPoint) = Me%Input%InitialWaterFluxY(i, j, k)
            Me%Input%BlockInMemory(1)%Discharges      (nPoint) = Me%Input%InitialDischarges(i, j, k)
            Me%Input%BlockInMemory(1)%ComputeFacesU3D (nPoint) = Me%Input%InitialComputeFacesU3D(i, j, k)
            Me%Input%BlockInMemory(1)%ComputeFacesV3D (nPoint) = Me%Input%InitialComputeFacesV3D(i, j, k)
        endif
        enddo
        enddo
        enddo

        do iInput = 2, NumberOfInputs
            call ReadsBlock 
            nPoint = 0
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
            if (Me%Input%WaterPoints2D(i, j) == WaterPoint) then
                nPoint = nPoint + 1
                Me%Input%BlockInMemory(iInput)%WaterLevel      (nPoint) = Me%Input%WaterLevel_New  (i, j)
                Me%Input%BlockInMemory(iInput)%WaterFluxX      (nPoint) = Me%Input%WaterFluxX      (i, j, k)
                Me%Input%BlockInMemory(iInput)%WaterFluxY      (nPoint) = Me%Input%WaterFluxY      (i, j, k)
                Me%Input%BlockInMemory(iInput)%Discharges      (nPoint) = Me%Input%Discharges      (i, j, k)
                Me%Input%BlockInMemory(iInput)%ComputeFacesU3D (nPoint) = Me%Input%ComputeFacesU3D (i, j, k)
                Me%Input%BlockInMemory(iInput)%ComputeFacesV3D (nPoint) = Me%Input%ComputeFacesV3D (i, j, k)
            endif
            enddo
            enddo
            enddo
        enddo

        Me%Input%NextBlockInMemory = 2

    end subroutine StartEncoding

    !----------------------------------------------------------------------------

    subroutine DecodeBlock 

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: iInput
        integer                                     :: nPoint
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        iInput = Me%Input%NextBlockInMemory
        nPoint = 0
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
        if (Me%Input%WaterPoints2D(i, j) == WaterPoint) then
            nPoint = nPoint + 1
            Me%Input%WaterLevel_New  (i, j)      = Me%Input%BlockInMemory(iInput)%WaterLevel      (nPoint)
            Me%Input%WaterFluxX      (i, j, k)   = Me%Input%BlockInMemory(iInput)%WaterFluxX      (nPoint)
            Me%Input%WaterFluxY      (i, j, k)   = Me%Input%BlockInMemory(iInput)%WaterFluxY      (nPoint)
            Me%Input%Discharges      (i, j, k)   = Me%Input%BlockInMemory(iInput)%Discharges      (nPoint)
            Me%Input%ComputeFacesU3D (i, j, k)   = Me%Input%BlockInMemory(iInput)%ComputeFacesU3D (nPoint)
            Me%Input%ComputeFacesV3D (i, j, k)   = Me%Input%BlockInMemory(iInput)%ComputeFacesV3D (nPoint)
        endif
        enddo
        enddo
        enddo


        Me%Input%NextBlockInMemory =                                    &
            Me%Input%NextBlockInMemory + 1

    end subroutine DecodeBlock

    !----------------------------------------------------------------------------

    subroutine StartHydrodynamicFileOutput

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer, dimension(:, :), pointer           :: BoundaryPoints2D
        integer                                     :: STAT_CALL

        !Sets Size

        Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KLB = Me%WorkSize%KLB
        Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KUB = Me%WorkSize%KUB

        Me%TimeIntegration%Integrate3DTime%WindowSize%KLB     = Me%Size%KLB
        Me%TimeIntegration%Integrate3DTime%WindowSize%KUB     = Me%Size%KUB

        if (.not.Me%TimeIntegration%Integrate3DTime%WindowOn) then

            Me%TimeIntegration%Integrate3DTime%WindowWorkSize%ILB = Me%WorkSize%ILB
            Me%TimeIntegration%Integrate3DTime%WindowWorkSize%IUB = Me%WorkSize%IUB

            Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JLB = Me%WorkSize%JLB
            Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JUB = Me%WorkSize%JUB

            Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KLB = Me%WorkSize%KLB
            Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KUB = Me%WorkSize%KUB

            Me%TimeIntegration%Integrate3DTime%WindowSize%ILB     = Me%Size%ILB
            Me%TimeIntegration%Integrate3DTime%WindowSize%IUB     = Me%Size%IUB

            Me%TimeIntegration%Integrate3DTime%WindowSize%JLB     = Me%Size%JLB
            Me%TimeIntegration%Integrate3DTime%WindowSize%JUB     = Me%Size%JUB

            Me%TimeIntegration%Integrate3DTime%WindowSize%KLB     = Me%Size%KLB
            Me%TimeIntegration%Integrate3DTime%WindowSize%KUB     = Me%Size%KUB
            
        end if

        if (Me%State%SpaceIntegration) then

            Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%ILB = 1
            Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%IUB = int((Me%TimeIntegration%Integrate3DTime%WindowWorkSize%IUB -          &
                                                                                            Me%TimeIntegration%Integrate3DTime%WindowWorkSize%ILB + 1)       &
                                                                                          / Me%SpaceIntegration%Integrate3DSpace%IntegrationStep)

            Me%SpaceIntegration%Integrate3DSpace%WindowSize%ILB     = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%ILB - 1
            Me%SpaceIntegration%Integrate3DSpace%WindowSize%IUB     = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%IUB + 1

        !-----------------------------
            Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%JLB = 1
            Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%JUB = int((Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JUB -          &
                                                                                            Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JLB + 1)       &
                                                                                          / Me%SpaceIntegration%Integrate3DSpace%IntegrationStep)

            Me%SpaceIntegration%Integrate3DSpace%WindowSize%JLB     = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%JLB - 1
            Me%SpaceIntegration%Integrate3DSpace%WindowSize%JUB     = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%JUB + 1

        !-----------------------------
            Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%KLB = Me%WorkSize%KLB
            Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%KUB = Me%WorkSize%KUB

            Me%SpaceIntegration%Integrate3DSpace%WindowSize%KLB     = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%KLB - 1
            Me%SpaceIntegration%Integrate3DSpace%WindowSize%KUB     = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%KUB + 1

        end if


        !Allocates input variables    
        call AllocateVariablesOutput

        !Starts Integration
        if (Me%State%TimeIntegration) then
            call StartHydroIntegration(Me%ObjHydroIntegration, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptionsOutput - ModuleHydrodynamicFile -ERR06'


            !Gets BoundaryPoints2D
            call GetBoundaries(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptionsOutput - ModuleHydrodynamicFile -ERR07'

            call StartHydroIntegrationList(ObjHydroIntegrationID  = Me%ObjHydroIntegration,         &
                   Size                   = Me%TimeIntegration%Integrate3DTime%WindowWorkSize,      &
                   DT                     = Me%TimeIntegration%Integrate3DTime%IntegrationStep,     &
                   CurrentTime            = Me%ExternalVar%Now,                                     &
                   DT_ComputeStep         = Me%ExternalVar%DT,                                      &
                   BoundaryPoints2D       = BoundaryPoints2D,                                       &
                   STAT                   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptionsOutput - ModuleHydrodynamicFile -ERR08'

            !Gets BoundaryPoints2D
            call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'HydrodynamicFileOptionsOutput - ModuleHydrodynamicFile -ERR09'

        end if

        !Writes new Bathymetry
        if (Me%TimeIntegration%Integrate3DTime%WindowOn .or. Me%State%SpaceIntegration) &
            call IntegrateBathymetry
                

        if (Me%State%SpaceIntegration) then
            call IntegrateInSpace (InitialFields = Me%Input)
            !Writes Header 
            call WriteHeader(IntegratedInitialFields = Me%SpaceIntegration%Integrate3DSpace)
        else
            !Writes Header 
            call WriteHeader(InitialFields = Me%Input)
        endif

        !Updates output time
        Me%TimeIntegration%NextOutput = Me%ExternalVar%Now + Me%TimeIntegration%Integrate3DTime%IntegrationStep


    end subroutine StartHydrodynamicFileOutput            

    !----------------------------------------------------------------------------

    subroutine AllocateVariablesInput

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: ILB, IUB
        integer                                     :: JLB, JUB
        integer                                     :: KLB, KUB

        !------------------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB


        !Initial
        nullify(Me%Input%InitialWaterLevel      )
        nullify(Me%Input%InitialWaterFluxX      )
        nullify(Me%Input%InitialWaterFluxY      )
        nullify(Me%Input%InitialDischarges      )
        nullify(Me%Input%InitialComputeFacesU3D )
        nullify(Me%Input%InitialComputeFacesV3D )
        nullify(Me%Input%WaterPoints2D          )


        !Transient
        nullify(Me%Input%WaterLevel_New    )
        nullify(Me%Input%WaterFluxX        )
        nullify(Me%Input%WaterFluxY        )
        nullify(Me%Input%Discharges        )
        nullify(Me%Input%ComputeFacesU3D   )
        nullify(Me%Input%ComputeFacesV3D   )


        !WaterLevel_New
        allocate(Me%Input%WaterLevel_New        (ILB:IUB, JLB:JUB         )) 
        allocate(Me%Input%InitialWaterLevel     (ILB:IUB, JLB:JUB         )) 
        allocate(Me%Input%WaterFluxX            (ILB:IUB, JLB:JUB, KLB:KUB)) 
        allocate(Me%Input%WaterFluxY            (ILB:IUB, JLB:JUB, KLB:KUB)) 
        allocate(Me%Input%Discharges            (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Input%ComputeFacesU3D       (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Input%ComputeFacesV3D       (ILB:IUB, JLB:JUB, KLB:KUB)) 
        allocate(Me%Input%InitialWaterFluxX     (ILB:IUB, JLB:JUB, KLB:KUB)) 
        allocate(Me%Input%InitialWaterFluxY     (ILB:IUB, JLB:JUB, KLB:KUB)) 
        allocate(Me%Input%InitialDischarges     (ILB:IUB, JLB:JUB, KLB:KUB)) 
        allocate(Me%Input%InitialComputeFacesU3D(ILB:IUB, JLB:JUB, KLB:KUB)) 
        allocate(Me%Input%InitialComputeFacesV3D(ILB:IUB, JLB:JUB, KLB:KUB)) 
        allocate(Me%Input%WaterPoints2D         (ILB:IUB, JLB:JUB         )) 

        Me%Input%WaterLevel_New(:,:)            = null_real
        Me%Input%InitialWaterLevel(:,:)         = null_real
        Me%Input%WaterFluxX(:,:,:)              = null_real
        Me%Input%WaterFluxY(:,:,:)              = null_real
        Me%Input%Discharges(:,:,:)              = null_real
        Me%Input%ComputeFacesU3D(:,:,:)         = 0
        Me%Input%ComputeFacesV3D(:,:,:)         = 0
        Me%Input%InitialWaterFluxX(:,:,:)       = null_real
        Me%Input%InitialWaterFluxY(:,:,:)       = null_real
        Me%Input%InitialDischarges(:,:,:)       = null_real
        Me%Input%InitialComputeFacesU3D(:,:,:)  = 0
        Me%Input%InitialComputeFacesV3D(:,:,:)  = 0
        Me%Input%WaterPoints2D(:,:)             = 0

        !----------------------------------------------------------------------

    end subroutine AllocateVariablesInput

    !--------------------------------------------------------------------------    

    subroutine AllocateVariablesOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: NewILB, NewIUB
        integer                                     :: NewJLB, NewJUB
        integer                                     :: KLB, KUB

        !----------------------------------------------------------------------

        if (Me%State%SpaceIntegration) then

            NewILB = Me%SpaceIntegration%Integrate3DSpace%WindowSize%ILB
            NewIUB = Me%SpaceIntegration%Integrate3DSpace%WindowSize%IUB

            NewJLB = Me%SpaceIntegration%Integrate3DSpace%WindowSize%JLB
            NewJUB = Me%SpaceIntegration%Integrate3DSpace%WindowSize%JUB

            KLB = Me%SpaceIntegration%Integrate3DSpace%WindowSize%KLB
            KUB = Me%SpaceIntegration%Integrate3DSpace%WindowSize%KUB

            nullify(Me%SpaceIntegration%Integrate3DSpace%WaterLevel     )
            nullify(Me%SpaceIntegration%Integrate3DSpace%WaterFluxX     )
            nullify(Me%SpaceIntegration%Integrate3DSpace%WaterFluxY     )
            nullify(Me%SpaceIntegration%Integrate3DSpace%Discharges     )
            nullify(Me%SpaceIntegration%Integrate3DSpace%ComputeFacesU3D)
            nullify(Me%SpaceIntegration%Integrate3DSpace%ComputeFacesV3D)
            nullify(Me%SpaceIntegration%Integrate3DSpace%WaterPoints2D  )

            !The next three variables are not output variables, but they are needed to make the integration of the fluxes in space

            allocate(Me%SpaceIntegration%VectorIntegrationY                 (NewILB:NewIUB))
            allocate(Me%SpaceIntegration%VectorIntegrationX                 (NewJLB:NewJUB))
            allocate(Me%SpaceIntegration%BatMin                             (NewILB:NewIUB, NewJLB:NewJUB))
            allocate(Me%SpaceIntegration%Integrate3DSpace%WaterLevel        (NewILB:NewIUB, NewJLB:NewJUB))
            allocate(Me%SpaceIntegration%Integrate3DSpace%WaterFluxX        (NewILB:NewIUB, NewJLB:NewJUB, KLB:KUB))
            allocate(Me%SpaceIntegration%Integrate3DSpace%WaterFluxY        (NewILB:NewIUB, NewJLB:NewJUB, KLB:KUB))
            allocate(Me%SpaceIntegration%Integrate3DSpace%Discharges        (NewILB:NewIUB, NewJLB:NewJUB, KLB:KUB))
            allocate(Me%SpaceIntegration%Integrate3DSpace%ComputeFacesU3D   (NewILB:NewIUB, NewJLB:NewJUB, KLB:KUB))
            allocate(Me%SpaceIntegration%Integrate3DSpace%ComputeFacesV3D   (NewILB:NewIUB, NewJLB:NewJUB, KLB:KUB))
            allocate(Me%SpaceIntegration%Integrate3DSpace%WaterPoints2D     (NewILB:NewIUB, NewJLB:NewJUB))

            Me%SpaceIntegration%VectorIntegrationY                  = null_int
            Me%SpaceIntegration%VectorIntegrationX                  = null_int
            Me%SpaceIntegration%VectorIntegrationY                  = null_int
            Me%SpaceIntegration%Integrate3DSpace%WaterLevel         = null_real
            Me%SpaceIntegration%Integrate3DSpace%WaterFluxX         = null_real
            Me%SpaceIntegration%Integrate3DSpace%WaterFluxY         = null_real
            Me%SpaceIntegration%Integrate3DSpace%Discharges         = null_real
            Me%SpaceIntegration%Integrate3DSpace%ComputeFacesU3D    = 0
            Me%SpaceIntegration%Integrate3DSpace%ComputeFacesV3D    = 0
            Me%SpaceIntegration%Integrate3DSpace%WaterPoints2D      = 0

        endif

    end subroutine AllocateVariablesOutput

    !----------------------------------------------------------------------------

    subroutine ReadHeader

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        real, dimension(6)                          :: TimeAux
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        integer                                     :: Unit
        integer                                     :: VersionNumber

        !------------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
        
        Unit            = Me%Files%InPutFieldsUnit
        VersionNumber   = Me%Files%InputVersionNumber

        !Rewinds unit
        rewind(Unit, IOSTAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR01.'


        !Reads DT_HYDROFILE
        read (Unit, IOSTAT = STAT_CALL) Me%Input%DT_HYDROFILE
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR03.'

         !StartTime
        read (Unit, IOSTAT = STAT_CALL) TimeAux
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR04.'
        call SetDate(Me%Input%StartTime,                                &
                     Year   = TimeAux(1),                                                &
                     Month  = TimeAux(2),                                                &
                     Day    = TimeAux(3),                                                &
                     Hour   = TimeAux(4),                                                &
                     Minute = TimeAux(5),                                                &
                     Second = TimeAux(6))


        !EndTime
        read (Unit, IOSTAT = STAT_CALL) TimeAux
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR05.'
        call SetDate(Me%Input%EndTime,                                  &
                     Year   = TimeAux(1),                                                &
                     Month  = TimeAux(2),                                                &
                     Day    = TimeAux(3),                                                &
                     Hour   = TimeAux(4),                                                &
                     Minute = TimeAux(5),                                                &
                     Second = TimeAux(6))


        !Reads WaterPoints2D
        if (VersionNumber > 1) then
            do j = JLB, JUB
            do i = ILB, IUB
                read(Unit, IOSTAT = STAT_CALL) Me%Input%WaterPoints2D(i, j)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR05a'
            enddo
            enddo
        endif

        !InitialWaterLevel
        do j=JLB, JUB
        do i=ILB, IUB
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(i, j) == 1) then
            read(Unit, IOSTAT = STAT_CALL) Me%Input%InitialWaterLevel(i, j)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR06.'
        endif
        end do
        end do

        !Initial WaterFluxX
        do k=KLB, KUB
        do j=JLB, JUB+1
        do i=ILB, IUB
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(i, min(j, JUB)) == 1) then
            read(Unit, IOSTAT = STAT_CALL) Me%Input%InitialWaterFluxX(i, j, k)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR07.'
        endif
        end do
        end do
        end do

        !Initial WaterFluxY
        do k=KLB, KUB
        do j=JLB, JUB
        do i=ILB, IUB+1
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(min(i, IUB), j) == 1) then
            read(Unit, IOSTAT = STAT_CALL) Me%Input%InitialWaterFluxY(i, j, k)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR08.'
        endif
        end do
        end do
        end do

        !Initial Discharges
        do k=KLB, KUB
        do j=JLB, JUB
        do i=ILB, IUB
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(i, j) == 1) then
            read(Unit, IOSTAT = STAT_CALL) Me%Input%InitialDischarges(i, j, k)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR09.'
        endif
        end do
        end do
        end do

        !Initial ComputeFacesU3D
        do k=KLB, KUB
        do j=JLB, JUB+1
        do i=ILB, IUB
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(i, min(j, JUB)) == 1) then
            read(Unit, IOSTAT = STAT_CALL) Me%Input%InitialComputeFacesU3D(i, j, k)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR10.'
        endif
        end do
        end do
        end do

        !Initial ComputeFacesV3D
        do k=KLB, KUB
        do j=JLB, JUB
        do i=ILB, IUB+1
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(min(i, IUB), j) == 1) then
            read(Unit, IOSTAT = STAT_CALL) Me%Input%InitialComputeFacesV3D(i, j, k)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeader - ModuleHydrodynamicFile - ERR11.'
        endif
        end do
        end do
        end do

        !------------------------------------------------------------------------

    end subroutine ReadHeader

    !--------------------------------------------------------------------------

    subroutine IntegrateBathymetry

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: NewILB, NewIUB
        integer                                     :: NewJLB, NewJUB
        integer                                     :: Zone
        integer                                     :: CoordType
        real, pointer, dimension(:,:)               :: NewBathymetry
        real, pointer, dimension(:  )               :: NewXX
        real, pointer, dimension(:  )               :: NewYY
        real                                        :: Xorig 
        real                                        :: Yorig 
        real                                        :: Latitude 
        real                                        :: Longitude 
        real                                        :: GRID_ANGLE
        integer                                     :: ILB, IUB, I
        integer                                     :: JLB, JUB, J
        integer                                     :: Iaux, Jaux
        integer                                     :: NIntegrationCells
        integer, pointer, dimension(:)              :: Vector_X, Vector_Y
        real,    pointer, dimension(:,:)            :: DUX, DVY
        integer                                     :: L, M
        real                                        :: CellArea, IntegratedCellArea
        real                                        :: haux, hmax, hmin
        real                                        :: RealAux
        integer                                     :: IntegerAux
        real,    pointer, dimension(:,:  )          :: Bathymetry
        real,    pointer, dimension(:    )          :: XX
        real,    pointer, dimension(:    )          :: YY
        type (T_Size2D)                             :: NewWorkSize
        
        !------------------------------------------------------------------------                         

        if(Me%State%SpaceIntegration) then
            NIntegrationCells = int(Me%SpaceIntegration%Integrate3DSpace%IntegrationStep)
            call GetHorizontalGrid(Me%ObjHorizontalGrid, DUX = DUX, DVY = DVY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR00'
        else
            NIntegrationCells = 1
        endif
        
        !Size of the window
        ILB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%ILB
        IUB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%IUB

        JLB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JLB
        JUB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JUB
        
        !Gets old Grid/Bathymetry
        call GetLatitudeLongitude(Me%ObjHorizontalGrid,                                  &
                                  Latitude  = Latitude,                                  &
                                  Longitude = Longitude,                                 &
                                  STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR01'

        call GetGridZone(Me%ObjHorizontalGrid, Zone = Zone, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR02'

        call GetGridAngle(Me%ObjHorizontalGrid, Angle = GRID_ANGLE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR03'

        call GetGridCoordType(Me%ObjHorizontalGrid, CoordType = CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR04'

        call GetGridOrigin(Me%ObjHorizontalGrid, Xorig = Xorig, Yorig = Yorig, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR05'

        call GetHorizontalGrid (Me%ObjHorizontalGrid, XX = XX, YY = YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR06'

        call GetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR07'


        !Limits of the new bathymetry
        NewILB    = 1
        NewIUB    = (IUB - ILB + 1) / NIntegrationCells

        NewJLB    = 1
        NewJUB    = (JUB - JLB + 1) / NIntegrationCells

        allocate(Vector_Y(NewILB - 1:NewIUB + 1))
        allocate(Vector_X(NewJLB - 1:NewJUB + 1))
       
        Vector_Y(NewILB) = ILB
        do I = NewILB, NewIUB
            Vector_Y(I + 1      ) = NIntegrationCells * I + ILB
        enddo
            Vector_Y(NewIUB + 1 ) = IUB + 1

            Vector_X(NewJLB) = JLB
        do J = NewJLB, NewJUB
            Vector_X(J + 1      ) = NIntegrationCells * J + JLB
        enddo
            Vector_X(NewJUB + 1 ) = JUB + 1

        if(Me%State%SpaceIntegration)then
            Me%SpaceIntegration%VectorIntegrationY = Vector_Y
            Me%SpaceIntegration%VectorIntegrationX = Vector_X
        endif
        


        !Allocates new Bathymetry
        allocate(NewBathymetry(NewILB:NewIUB+1, NewJLB:NewJUB+1))
        allocate(NewXX        (                 NewJLB:NewJUB+2))
        allocate(NewYY        (NewILB:NewIUB+2                 ))

        !Computes new XX, YY
        
        Iaux        = 1
        NewYY(Iaux) = 0.
do1 :   do I = NewILB, NewIUB
            Iaux        = Iaux + 1
            NewYY(Iaux) = YY(Vector_Y(I+1)) - YY(Vector_Y(NewILB))
        end do do1


        Jaux        = 1
        NewXX(Jaux) = 0.
do4 :   do J = NewJLB, NewJUB
            Jaux        = Jaux + 1
            NewXX(Jaux) = XX(Vector_X(J+1)) - XX(Vector_X(NewJLB))
        end do do4
              
        
cd1 :   if (Me%State%SpaceIntegration) then

cd2 :       if(Me%SpaceIntegration%BatIntegrationType /= MeanVal_Type)then
                !Fills new integrated bathymetry using the maximum value
                Iaux = 1
                Jaux = 1
                do J = NewJLB, NewJUB
                    do I = NewILB, NewIUB

                        IntegerAux = 0
                        do M = Vector_X(J), Vector_X(J + 1) - 1
                            do L = Vector_Y(I), Vector_Y(I + 1) - 1
                                if(Bathymetry(L, M) > -55.0) IntegerAux = 1
                            end do
                        end do

                        if(IntegerAux == 1)then

                            do M = Vector_X(J), Vector_X(J + 1) - 1
                                do L = Vector_Y(I), Vector_Y(I + 1) - 1
                                    if(Bathymetry(L, M) > -55.0)then
                                        haux = Bathymetry(L, M)
                                        if(M == Vector_X(J))then
                                            hmax = haux
                                            hmin = haux
                                        else
                                            if(haux > hmax) hmax = haux
                                            if(haux < hmin) hmin = haux
                                        endif
                                    endif
                                enddo
                            enddo
                            
                            IntegratedCellArea = 0.0
                            RealAux = 0.0
                            do M = Vector_X(J), Vector_X(J + 1) - 1
                                do L = Vector_Y(I), Vector_Y(I + 1) - 1
                                    CellArea = DUX(L, M) * DVY(L, M)
                                    if(Bathymetry(L, M) > -55.0)then
                                        RealAux = RealAux + CellArea * hmax
                                    else
                                        RealAux = RealAux + CellArea * hmin
                                    endif
                                    IntegratedCellArea = IntegratedCellArea + CellArea
                                end do
                            end do

                            NewBathymetry(iaux, jaux) = RealAux / IntegratedCellArea
                            Me%SpaceIntegration%BatMin(iaux, jaux) = hmin

                        else

                            NewBathymetry(iaux, jaux) = -99.00000

                        endif

                        Iaux = Iaux + 1
                    end do
                    Iaux = 1
                    Jaux = Jaux + 1
                end do

            else cd2
                !Fills new integrated bathymetry using the mean value
                Iaux = 1
                Jaux = 1
                do J = NewJLB, NewJUB
                    do I = NewILB, NewIUB

                        IntegerAux = 0
                        do M = Vector_X(J), Vector_X(J + 1) - 1
                            do L = Vector_Y(I), Vector_Y(I + 1) - 1
                                if(Bathymetry(L, M) > -55.0) IntegerAux = 1
                            end do
                        end do

                        if(IntegerAux == 1)then

                            do M = Vector_X(J), Vector_X(J + 1) - 1
                                do L = Vector_Y(I), Vector_Y(I + 1) - 1
                                    if(Bathymetry(L, M) > -55.0)then
                                        haux = Bathymetry(L, M)
                                        if(M == Vector_X(J))then
                                            hmin = haux
                                        else
                                            if(haux < hmin) hmin = haux
                                        endif
                                    endif
                                enddo
                            enddo
                            
                            IntegratedCellArea = 0.0
                            RealAux = 0.0
                            do M = Vector_X(J), Vector_X(J + 1) - 1
                                do L = Vector_Y(I), Vector_Y(I + 1) - 1
                                    CellArea = DUX(L, M) * DVY(L, M)
                                    if(Bathymetry(L, M) > -55.0)then
                                        RealAux = RealAux + CellArea * Bathymetry(L, M)
                                    else
                                        RealAux = RealAux + CellArea * hmin
                                    endif
                                    IntegratedCellArea = IntegratedCellArea + CellArea
                                enddo
                            enddo

                            NewBathymetry(iaux, jaux) = RealAux / IntegratedCellArea
                            Me%SpaceIntegration%BatMin(iaux, jaux) = hmin

                        else

                            NewBathymetry(iaux, jaux) = -99.00000

                        endif

                        Iaux = Iaux + 1
                    end do
                    Iaux = 1
                    Jaux = Jaux + 1
                end do
            
            endif cd2

        else cd1        

        !Fills new Bathymetry
        Iaux = 1
        Jaux = 1
        do J = JLB, JUB
            do I = ILB, IUB
                NewBathymetry(Iaux, Jaux) = Bathymetry(I, J)
                Iaux = Iaux + 1
            end do
            Iaux = 1
            Jaux = Jaux + 1
        end do

        endif cd1

        
        !Moves Origin
        Xorig = Xorig + XX(JLB)
        Yorig = Yorig + YY(ILB)

        NewWorkSize = T_Size2D(NewILB, NewIUB, NewJLB, NewJUB)
        call WriteGridData  (Me%Files%IntegratedBathymetry,                         &
                             XX = NewXX, YY = NewYY,                                &
                             COMENT1  = 'New Integrated Bathymetry',                &
                             COMENT2  = 'Generated by HydrodynamicFile',            &
                             WorkSize = NewWorkSize,                                &
                             CoordType= CoordType,                                  &
                             Xorig    = Xorig, Yorig = Yorig,                       &
                             Zone     = Zone, GRID_ANGLE = GRID_ANGLE,              &
                             Latitude = Latitude, Longitude = Longitude,            &
                             FillValue  = -99.,                                     &
                             GridData2D_Real = NewBathymetry,                       &
                             Overwrite  = OFF,                                      &
                             STAT = STAT_CALL)                                                                           
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR08'

        

        deallocate(Vector_X     )                    
        deallocate(Vector_Y     )                    
        deallocate(NewBathymetry)                    
        deallocate(NewXX        )                    
        deallocate(NewYY        )                    

        if((Me%State%SpaceIntegration) .and.                             &
           (Me%SpaceIntegration%BatIntegrationType == MeanVal_Type))then

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DUX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR09'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DVY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR10'

        endif

        call UngetGridData (Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR11'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, XX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR12'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'IntegrateBathymetry - ModuleHydrodynamicFile - ERR13'

        !------------------------------------------------------------------------

    end subroutine IntegrateBathymetry

    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------

    subroutine WriteHeader(InitialFields, IntegratedInitialFields)

        !Arguments---------------------------------------------------------------
        type(T_STIntegration), optional, intent(IN)  :: IntegratedInitialFields
        type(T_Input),         optional, intent(IN)  :: InitialFields

        !External----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_Time)                               :: BeginTime, EndTime
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        integer                                     :: Unit
        integer                                     :: VersionNumber
        real,    pointer, dimension(:,:  )          :: InitialWaterLevel
        real(8), pointer, dimension(:,:,:)          :: InitialWaterFluxX
        real(8), pointer, dimension(:,:,:)          :: InitialWaterFluxY
        real(8), pointer, dimension(:,:,:)          :: InitialDischarges
        integer, dimension(:,:,:), pointer          :: InitialComputeFacesU3D
        integer, dimension(:,:,:), pointer          :: InitialComputeFacesV3D
        integer, dimension(:,:  ), pointer          :: WaterPoints2D

        !------------------------------------------------------------------------                         

        VersionNumber = Me%Files%OutputVersionNumber

        if (present(IntegratedInitialFields)) then

            ILB = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%ILB
            IUB = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%IUB
        
            JLB = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%JLB
            JUB = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%JUB
        
            KLB = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%KLB
            KUB = Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%KUB


            InitialWaterLevel      => IntegratedInitialFields%WaterLevel
            InitialWaterFluxX      => IntegratedInitialFields%WaterFluxX
            InitialWaterFluxY      => IntegratedInitialFields%WaterFluxY
            InitialDischarges      => IntegratedInitialFields%Discharges
            InitialComputeFacesU3D => IntegratedInitialFields%ComputeFacesU3D
            InitialComputeFacesV3D => IntegratedInitialFields%ComputeFacesV3D
            WaterPoints2D          => IntegratedInitialFields%WaterPoints2D

        else if (present(InitialFields)) then

            ILB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%ILB
            IUB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%IUB
        
            JLB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JLB
            JUB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JUB
        
            KLB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KLB
            KUB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KUB


            InitialWaterLevel      => InitialFields%InitialWaterLevel
            InitialWaterFluxX      => InitialFields%InitialWaterFluxX
            InitialWaterFluxY      => InitialFields%InitialWaterFluxY
            InitialDischarges      => InitialFields%InitialDischarges
            InitialComputeFacesU3D => InitialFields%InitialComputeFacesU3D
            InitialComputeFacesV3D => InitialFields%InitialComputeFacesV3D

            call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR00'

        endif

        !Gets Unit for the output file
        call UnitsManager(Me%Files%OutPutFieldsUnit, OPEN_FILE, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR01'

        Unit = Me%Files%OutPutFieldsUnit

        !Opens the file
        open(UNIT   = Unit,                                                              &
             FILE   = Me%Files%OutPutFields,                                             &
             FORM   ='UNFORMATTED',                                                      &
             STATUS ='UNKNOWN',                                                          &
             ACTION ='WRITE',                                                            &
             IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR02'


        !Gets Time Limits
        call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, &
                                  EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR03'

        !DT_HYDROFILE
        write(Unit, IOSTAT = STAT_CALL) Me%TimeIntegration%Integrate3DTime%IntegrationStep
        if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR04'


        !StartTime
        call ExtractDate(BeginTime, Year, Month, Day, Hour, Minute, Second)
        write(Unit, IOSTAT = STAT_CALL) Year, Month, Day, Hour, Minute, Second
        if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR05'


        !EndTime
        call ExtractDate(EndTime,   Year, Month, Day, Hour, Minute, Second)
        write(Unit, IOSTAT = STAT_CALL) Year, Month, Day, Hour, Minute, Second
        if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR06'


        !Writes WaterPoints2D
        if (VersionNumber > 1) then
            do j = JLB, JUB
            do i = ILB, IUB
                write(Unit, IOSTAT = STAT_CALL) WaterPoints2D(i, j)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR07'
            enddo
            enddo
        endif


        !Writes initial Waterlevel
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints2D(i, j) == WaterPoint .or. VersionNumber == 1) then
                write(Unit, IOSTAT = STAT_CALL) InitialWaterLevel(i, j)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR08'
            endif
        enddo
        enddo

        !Initial WaterFluxX
        do k=KLB, KUB
        do j=JLB, JUB+1
        do i=ILB, IUB
            if (WaterPoints2D(i, min(j, JUB)) == WaterPoint .or. VersionNumber == 1) then
                write(Unit, IOSTAT = STAT_CALL) InitialWaterFluxX(i, j, k)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR09'
            end if
        end do
        end do
        end do

        !Initial WaterFluxY
        do k=KLB, KUB
        do j=JLB, JUB
        do i=ILB, IUB+1
            if (WaterPoints2D(min(i, IUB), j) == WaterPoint .or. VersionNumber == 1) then
                write(Unit, IOSTAT = STAT_CALL) InitialWaterFluxY(i, j, k)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR10'
            end if
        end do
        end do
        end do

        !Initial Discharges
        do k=KLB, KUB
        do j=JLB, JUB
        do i=ILB, IUB
            if (WaterPoints2D(i, j) == WaterPoint .or. VersionNumber == 1) then
                write(Unit, IOSTAT = STAT_CALL) InitialDischarges(i, j, k)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR11'
            end if
        end do
        end do
        end do

        !Initial ComputeFacesU3D
        do k=KLB, KUB
        do j=JLB, JUB+1
        do i=ILB, IUB
            if (WaterPoints2D(i, min(j, JUB)) == WaterPoint .or. VersionNumber == 1) then
                write(Unit, IOSTAT = STAT_CALL) InitialComputeFacesU3D(i, j, k)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR12'
            end if
        end do
        end do
        end do

        !Initial ComputeFacesV3D
        do k=KLB, KUB
        do j=JLB, JUB
        do i=ILB, IUB+1
            if (WaterPoints2D(min(i, IUB), j) == WaterPoint .or. VersionNumber == 1) then
                write(Unit, IOSTAT = STAT_CALL) InitialComputeFacesV3D(i, j, k)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR13'
            end if
        end do
        end do
        end do

        if (present(InitialFields)) then
            call UnGetHorizontalMap (Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHeader - ModuleHydrodynamicFile - ERR14'
        endif
        !------------------------------------------------------------------------

    end subroutine WriteHeader


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetHydrodynamicFileIOState(InputStateEx, OutputStateEx)

        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: InputStateEx
        integer, optional, intent(OUT) :: OutputStateEx

        !----------------------------------------------------------------------

        if (present(InputStateEx )) InputStateEx  = InputState
        if (present(OutputStateEx)) OutputStateEx = OutputState

        !----------------------------------------------------------------------

    end subroutine GetHydrodynamicFileIOState

    !--------------------------------------------------------------------------

    subroutine GetFileWaterLevel(HydrodynamicFileID, WaterLevel, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        real, dimension(:,:),   pointer             :: WaterLevel
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHYDRODYNAMICFILE_, Me%InstanceID)
            WaterLevel => Me%Input%WaterLevel_New

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetFileWaterLevel

    !--------------------------------------------------------------------------

    subroutine GetFileFluxes(HydrodynamicFileID, WaterFluxX, WaterFluxY,                &
                             Discharges, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        real(8), dimension(:,:,:),   pointer        :: WaterFluxX, WaterFluxY, Discharges
        integer, optional, intent(OUT)              :: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHYDRODYNAMICFILE_, Me%InstanceID)
            WaterFluxX  => Me%Input%WaterFluxX

            call Read_Lock(mHYDRODYNAMICFILE_, Me%InstanceID)
            WaterFluxY  => Me%Input%WaterFluxY

            call Read_Lock(mHYDRODYNAMICFILE_, Me%InstanceID)
            Discharges  => Me%Input%Discharges

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetFileFluxes

    !--------------------------------------------------------------------------

    subroutine GetFileMapping(HydrodynamicFileID, ComputeFacesU3D, ComputeFacesV3D, &
                              STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        integer, dimension(:, :, :), pointer        :: ComputeFacesU3D
        integer, dimension(:, :, :), pointer        :: ComputeFacesV3D
        integer, optional, intent (OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_


        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR.  (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHYDRODYNAMICFILE_, Me%InstanceID)
            ComputeFacesU3D => Me%Input%ComputeFacesU3D
            
            call Read_Lock(mHYDRODYNAMICFILE_, Me%InstanceID)
            ComputeFacesV3D => Me%Input%ComputeFacesV3D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine GetFileMapping

    !--------------------------------------------------------------------------

    subroutine UngetHydrodynamicFile3Dreal4(HydrodynamicFileID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        real(4), pointer, dimension(:,:,:)          :: Array
        integer, optional, intent (OUT)             :: STAT
   
        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mHYDRODYNAMICFILE_, Me%InstanceID, "UngetHydrodynamicFile3Dreal4")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHydrodynamicFile3Dreal4

    !--------------------------------------------------------------------------

    subroutine UngetHydrodynamicFile2Dreal4(HydrodynamicFileID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        real(4), pointer, dimension(:,:)            :: Array
        integer, optional, intent (OUT)             :: STAT
   
        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mHYDRODYNAMICFILE_, Me%InstanceID, "UngetHydrodynamicFile2Dreal4")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHydrodynamicFile2Dreal4

    !--------------------------------------------------------------------------

    subroutine UngetHydrodynamicFile2Dreal8(HydrodynamicFileID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        real(8), pointer, dimension(:,:)            :: Array
        integer, optional, intent (OUT)             :: STAT
   
        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mHYDRODYNAMICFILE_, Me%InstanceID, "UngetHydrodynamicFile2Dreal8")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHydrodynamicFile2Dreal8

    !--------------------------------------------------------------------------

    subroutine UngetHydrodynamicFile3Dreal8(HydrodynamicFileID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        real(8), pointer, dimension(:,:,:)          :: Array
        integer, optional, intent (OUT)             :: STAT
   
        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mHYDRODYNAMICFILE_, Me%InstanceID, "UngetHydrodynamicFile3Dreal8")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHydrodynamicFile3Dreal8

    !--------------------------------------------------------------------------

    subroutine UngetHydrodynamicFile3Dinteger(HydrodynamicFileID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        integer, pointer, dimension(:,:,:)          :: Array
        integer, optional, intent (OUT)             :: STAT
   
        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mHYDRODYNAMICFILE_, Me%InstanceID, "UngetHydrodynamicFile3Dinteger")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHydrodynamicFile3Dinteger


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine ModifyHydrodynamicFileInput(HydrodynamicFileID, HydrodynamicTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        type (T_Time), optional                     :: HydrodynamicTime
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_CALL
        integer                                     :: STAT_
        integer                                     :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

cd2 :       if (present(HydrodynamicTime)) then

                Me%ExternalVar%Now = HydrodynamicTime

            else cd2

                ! Actualized the time
                call GetComputeCurrentTime(Me%ObjTime,         &
                                           Me%ExternalVar%Now, &
                                           STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'Subroutine ModifyHydrodynamicFileInput; Module ModuleHydrodynamicFile. ERR01' 
            end if cd2



            if (Me%Input%File_Type == M2_Tide_type   ) then

                call ModifyMohid3d_M2

            else if (Me%Input%File_Type == BeginEnd_type  ) then

                if (Me%ExternalVar%Now .GE. Me%Input%NextInput) then
                    
                    if (Me%Input%LoadAllToMemory) then
                        
                        call DecodeBlock 

                    else                    

                        !Reads a new block of information
                        call ReadsBlock  

                    endif

                    !Actualizes next output
                    Me%Input%NextInput = Me%ExternalVar%Now   &
                                                         + Me%Input%DT_HYDROFILE
                endif
            endif



            call null_time   (Me%ExternalVar%Now)
           
            STAT_ = SUCCESS_

        else cd1
                      
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ModifyHydrodynamicFileInput

    !----------------------------------------------------------------------------

    subroutine ModifyMohid3d_M2

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------
        type (T_Time)                               :: ActualTime
        real                                        :: DT, DT_HYDROFILE

        !Local-------------------------------------------------------------------
        integer                                     :: M2_It
        integer                                     :: ILB, IUB
        integer                                     :: JLB, JUB

        !------------------------------------------------------------------------

        !Shorten the names variables 
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


        ActualTime   = Me%ExternalVar%Now        
        DT           = Me%ExternalVar%DT
        DT_HYDROFILE = Me%Input%DT_HYDROFILE
        M2_It        = null_int


cd1 :   if ((ActualTime + DT) .GT.  Me%Input%StartTime) then

            M2_It = M2_Iteration(ActualTime, Me%Input%StartTime, DT_HYDROFILE)

        else cd1

            M2_It = -100

        end if cd1


        !If reached a point where M2_It is equal to one, its necessary to rewind the file 
        !and reread the header.
        if (M2_It .EQ. 1) then
            if (Me%Input%LoadAllToMemory) then
                Me%Input%NextBlockInMemory = 1
            else
                call ReadHeader
            endif
        endif

        if (Me%Input%LoadAllToMemory) then
            call DecodeBlock
        else
            !Read one block of information
            call ReadsBlock
        endif

        !------------------------------------------------------------------------      

    end subroutine ModifyMohid3d_M2

    !----------------------------------------------------------------------------

    integer function M2_Iteration(Time,TimeBegin,DT_Seconds)

        !Arguments---------------------------------------------------------------

        type (T_Time)       :: Time,TimeBegin
        real, intent(in)    :: DT_Seconds

        !Local-------------------------------------------------------------------

        type (T_Time) :: Time1,Time2
        real          :: DT_Interval,NumberPeriods,NumberOfOutPuts

        !------------------------------------------------------------------------

        Time1 = Time

        Time2 = TimeBegin

        if (Time1.EQ.Time2) stop 'error in - function M2_Iteration - M2_Hydrodynamic'
    
        DT_Interval=Time1-Time2

        NumberPeriods = DT_Interval/(M2_Period_Hours*3600.)

        NumberOfOutPuts = M2_Period_Hours*3600./DT_Seconds

        M2_Iteration  = nint((NumberPeriods-int(NumberPeriods)) * NumberOfOutPuts)

        if (M2_Iteration==0) M2_Iteration = NumberOfOutPuts

        !------------------------------------------------------------------------

    End function M2_Iteration

    !----------------------------------------------------------------------------

    subroutine ReadsBlock

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: ILB, IUB, I
        integer                                     :: JLB, JUB, J
        integer                                     :: KLB, KUB, K
        integer                                     :: Unit, VersionNumber
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second

        integer                                     :: AuxBlock, ios

        !------------------------------------------------------------------------

        ILB       = Me%WorkSize%ILB
        IUB       = Me%WorkSize%IUB

        JLB       = Me%WorkSize%JLB
        JUB       = Me%WorkSize%JUB
         
        KLB       = Me%WorkSize%KLB
        KUB       = Me%WorkSize%KUB

        Unit            = Me%Files%InPutFieldsUnit
        VersionNumber   = Me%Files%InputVersionNumber


        !Reads begin block
        read(Unit, err = 50, iostat = ios) AuxBlock
 50     if (ios /= 0 .or. AuxBlock /= BeginBlock) then
            stop 'ReadsBlock - ModuleHydrodynamicFile - ERR01' 
        end if

        

        !Reads the current time
        read(Unit) Year, Month, Day, Hour, Minute, Second
        call SetDate(Me%Input%TransientTime,                                             &
                     Year   = Year,                                                      &
                     Month  = Month,                                                     &
                     Day    = Day,                                                       &
                     Hour   = Hour,                                                      &
                     Minute = Minute,                                                    &
                     Second = Second)



        !Reads the Waterlevel
        do j = JLB, JUB
        do i = ILB, IUB
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(i, j) == 1) then
            read(Unit) Me%Input%WaterLevel_New(i, j)
        endif
        enddo
        enddo

        !Read WaterFlux X
        do k = KLB, KUB
        do j = JLB, JUB + 1
        do i = ILB, IUB
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(i, min(j, JUB)) == 1) then
            read(Unit) Me%Input%WaterFluxX(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Reads WaterFlux Y
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB + 1
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(min(i, IUB), j) == 1) then
            read(Unit) Me%Input%WaterFluxY(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Reads Discharges
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(i, j) == 1) then
            read(Unit) Me%Input%Discharges(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Writes ComputeFacesU3D
        do k = KLB, KUB
        do j = JLB, JUB + 1
        do i = ILB, IUB
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(i, min(j, JUB)) == 1) then
            read(Unit) Me%Input%ComputeFacesU3D(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Writes ComputeFacesV3D
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB + 1
        if (VersionNumber == 1 .or. Me%Input%WaterPoints2D(min(i, IUB), j) == 1) then
            read(Unit) Me%Input%ComputeFacesV3D(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Reads end block
        read(Unit)AuxBlock

        if (AuxBlock /= EndBlock) then
            stop 'ReadsBlock - ModuleHydrodynamicFile - ERR02' 
        endif


        !------------------------------------------------------------------------

    end subroutine ReadsBlock

    !----------------------------------------------------------------------------

    subroutine ModifyHydrodynamicFileOutput(HydrodynamicFileID,                         &
                                            WaterLevel,                                 &
                                            WaterFluxX,                                 &
                                            WaterFluxY,                                 &
                                            Discharges,                                 &
                                            ComputeFacesU3D,                            &
                                            ComputeFacesV3D,                            &
                                            HydrodynamicTime,                           &
                                            STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HydrodynamicFileID
        real,    dimension(:,:),   pointer          :: WaterLevel
        real(8), dimension(:,:,:), pointer          :: WaterFluxX
        real(8), dimension(:,:,:), pointer          :: WaterFluxY
        real(8), dimension(:,:,:), pointer          :: Discharges
        integer, dimension(:,:,:), pointer          :: ComputeFacesU3D
        integer, dimension(:,:,:), pointer          :: ComputeFacesV3D
        type (T_Time), optional                     :: HydrodynamicTime
        integer, optional, intent(OUT)              :: STAT

        !External----------------------------------------------------------------

        integer                                     :: ready_ 
        integer                                     :: STAT_CALL             

        !Local-------------------------------------------------------------------

        integer                                     :: STAT_
        
        integer                                     :: ILB, IUB
        integer                                     :: JLB, JUB
        integer                                     :: KLB, KUB
        real(8), dimension(:, :, :), pointer        :: VolumeZ, VolumeZOld
        integer, dimension(:, :, :), pointer        :: WaterPoints3D

        !------------------------------------------------------------------------                

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%ILB
            IUB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%IUB

            JLB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JLB
            JUB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%JUB
         
            KLB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KLB
            KUB = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KUB

cd2 :       if (present(HydrodynamicTime)) then

                Me%ExternalVar%Now = HydrodynamicTime

            else cd2

                ! Actualized the time
                call GetComputeCurrentTime(Me%ObjTime,                                  &
                                           Me%ExternalVar%Now,                          &
                                           STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR01' 
            end if cd2

            
            nullify(Me%TimeIntegration%Integrate3DTime%WaterLevel)
            nullify(Me%TimeIntegration%Integrate3DTime%WaterFluxX)
            nullify(Me%TimeIntegration%Integrate3DTime%WaterFluxY)
            nullify(Me%TimeIntegration%Integrate3DTime%Discharges)
            nullify(Me%TimeIntegration%Integrate3DTime%ComputeFacesU3D)
            nullify(Me%TimeIntegration%Integrate3DTime%ComputeFacesV3D)

            Me%TimeIntegration%Integrate3DTime%WaterLevel      => WaterLevel
            Me%TimeIntegration%Integrate3DTime%WaterFluxX      => WaterFluxX
            Me%TimeIntegration%Integrate3DTime%WaterFluxY      => WaterFluxY
            Me%TimeIntegration%Integrate3DTime%Discharges      => Discharges
            Me%TimeIntegration%Integrate3DTime%ComputeFacesU3D => ComputeFacesU3D
            Me%TimeIntegration%Integrate3DTime%ComputeFacesV3D => ComputeFacesV3D

cd4  :      if (Me%State%TimeIntegration) then

                !Get VolumeZ, VolumeZOld
                call GetGeometryVolumes(Me%ObjGeometry,                                 &
                                        VolumeZ    = VolumeZ,                           &
                                        VolumeZOld = VolumeZOld,                        &
                                        STAT       = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR10' 

                call GetWaterPoints3D(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
                
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR20' 


                !Do one integration step
                call ModifyHydroIntegration(Me%ObjHydroIntegration,                     &
                                            WaterFluxX, WaterFluxY, Discharges,         &
                                            ComputeFacesU3D, ComputeFacesV3D,           &
                                            WaterPoints3D, VolumeZ, VolumeZOld,         &
                                            Me%ExternalVar%Now,                         &
                                            STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR30' 

                !Unget VolumeZ
                call UnGetGeometry(Me%ObjGeometry, VolumeZ,                             &
                                   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR40' 

                !Unget VolumeZOld
                call UnGetGeometry(Me%ObjGeometry, VolumeZOld,                          &
                                   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR50' 

                call UnGetMap     (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
                
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR55' 



            end if cd4


cd3 :       if (Me%ExternalVar%Now .GE. Me%TimeIntegration%NextOutput) then

                if (Me%State%TimeIntegration) then

                    !Get integrated WaterFluxes
                    call GetHydroIntegrationWaterFluxes(                                &
                            Me%ObjHydroIntegration,                                     &
                            Me%TimeIntegration%Integrate3DTime%IntegrationStep,         &
                            Me%TimeIntegration%Integrate3DTime%WindowWorkSize,          &
                            Me%TimeIntegration%Integrate3DTime%WaterFluxX,              &
                            Me%TimeIntegration%Integrate3DTime%WaterFluxY,              &
                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR60' 


                    !Get integrated Mapping
                    call GetHydroIntegrationComputeFaces(                               &
                            Me%ObjHydroIntegration,                                     &
                            Me%TimeIntegration%Integrate3DTime%IntegrationStep,         &
                            Me%TimeIntegration%Integrate3DTime%WindowWorkSize,          &
                            Me%TimeIntegration%Integrate3DTime%ComputeFacesU3D,         &
                            Me%TimeIntegration%Integrate3DTime%ComputeFacesV3D,         &

                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR70' 


                    !Get integrated Discharges
                    call GetHydroIntegrationDischarges(                                 &
                            Me%ObjHydroIntegration,                                     &
                            Me%TimeIntegration%Integrate3DTime%IntegrationStep,         &
                            Me%TimeIntegration%Integrate3DTime%WindowWorkSize,          &
                            Me%TimeIntegration%Integrate3DTime%Discharges,              &
                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR80' 


                endif

                

                if (Me%State%SpaceIntegration) then

                    call IntegrateInSpace

                endif

                
                
                !Writes the data file
                if (Me%State%SpaceIntegration) then

                    call WriteBlock(Me%SpaceIntegration%Integrate3DSpace)
                else

                    call WriteBlock(Me%TimeIntegration%Integrate3DTime)
                endif

                !Actualizes next output    
                Me%TimeIntegration%NextOutput = Me%TimeIntegration%NextOutput  +        &
                                                Me%TimeIntegration%Integrate3DTime%IntegrationStep


                !UnGets pointers                
                if (Me%State%TimeIntegration) then

                    !WaterFluxes
                    call UnGetHydroIntegration(                                         &
                            Me%ObjHydroIntegration,                                     &
                            Me%TimeIntegration%Integrate3DTime%WaterFluxX,              &
                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR90' 
                    
                    call UnGetHydroIntegration(                                         &
                            Me%ObjHydroIntegration,                                     &
                            Me%TimeIntegration%Integrate3DTime%WaterFluxY,              &
                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR100' 

                    !Mapping
                    call UnGetHydroIntegration(                                         &
                            Me%ObjHydroIntegration,                                     &
                            Me%TimeIntegration%Integrate3DTime%ComputeFacesU3D,         &
                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR110' 

                    call UnGetHydroIntegration(                                         &
                            Me%ObjHydroIntegration,                                     &
                            Me%TimeIntegration%Integrate3DTime%ComputeFacesV3D,         &
                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR120' 

                    !Discharges
                    call UnGetHydroIntegration(                                         &
                            Me%ObjHydroIntegration,                                     &
                            Me%TimeIntegration%Integrate3DTime%Discharges,              &
                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Subroutine ModifyHydrodynamicFileOutput; Module ModuleHydrodynamicFile. ERR130' 


                endif


            end if cd3
     

            call null_time   (Me%ExternalVar%Now)

            STAT_ = SUCCESS_

        else cd1
                      
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                              &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine ModifyHydrodynamicFileOutput

    !----------------------------------------------------------------------------

    subroutine IntegrateInSpace(InitialFields)

        !Arguments---------------------------------------------------------------
        type(T_Input), optional, intent(IN)         :: InitialFields

        !Local-------------------------------------------------------------------
        integer                              :: STAT_CALL

        real,    pointer, dimension(:,:  )   :: SIWaterLevel
        real(8), pointer, dimension(:,:,:)   :: SIWaterFluxX
        real(8), pointer, dimension(:,:,:)   :: SIWaterFluxY
        real(8), pointer, dimension(:,:,:)   :: SIDischarges
        integer, dimension(:,:,:), pointer   :: SIComputeFacesU3D
        integer, dimension(:,:,:), pointer   :: SIComputeFacesV3D
        integer, dimension(:,:  ), pointer   :: SIWaterPoints2D

        real,    pointer, dimension(:,:  )   :: TIWaterLevel
        real(8), pointer, dimension(:,:,:)   :: TIWaterFluxX
        real(8), pointer, dimension(:,:,:)   :: TIWaterFluxY
        real(8), pointer, dimension(:,:,:)   :: TIDischarges
        integer, dimension(:,:,:), pointer   :: TIComputeFacesU3D
        integer, dimension(:,:,:), pointer   :: TIComputeFacesV3D

        real,    pointer, dimension(:,:)     :: DUX, DVY

        integer, pointer, dimension(:)       :: Vector_X, Vector_Y

        integer, pointer, dimension(:,:,:)   :: OpenPoints
        integer, pointer, dimension(:,:,:)   :: WaterPoints3D

        real                                 :: CellArea, IntegratedCellArea

        integer                              :: i, j, l, m
        integer                              :: NewILB, NewIUB, NewJLB, NewJUB, KUB

        integer                              :: SumComputeFacesU3D, SumComputeFacesV3D

        integer                              :: IntegerAux
        real(8)                              :: Real8Aux
        real                                 :: RealAux


        !------------------------------------------------------------------------
        SIWaterLevel      => Me%SpaceIntegration%Integrate3DSpace%WaterLevel
        SIWaterFluxX      => Me%SpaceIntegration%Integrate3DSpace%WaterFluxX
        SIWaterFluxY      => Me%SpaceIntegration%Integrate3DSpace%WaterFluxY
        SIDischarges      => Me%SpaceIntegration%Integrate3DSpace%Discharges
        SIComputeFacesU3D => Me%SpaceIntegration%Integrate3DSpace%ComputeFacesU3D
        SIComputeFacesV3D => Me%SpaceIntegration%Integrate3DSpace%ComputeFacesV3D
        SIWaterPoints2D   => Me%SpaceIntegration%Integrate3DSpace%WaterPoints2D

        if (present(InitialFields)) then

            TIWaterLevel      => InitialFields%InitialWaterLevel
            TIWaterFluxX      => InitialFields%InitialWaterFluxX
            TIWaterFluxY      => InitialFields%InitialWaterFluxY
            TIDischarges      => InitialFields%InitialDischarges
            TIComputeFacesU3D => InitialFields%InitialComputeFacesU3D
            TIComputeFacesV3D => InitialFields%InitialComputeFacesV3D

        else

            TIWaterLevel      => Me%TimeIntegration%Integrate3DTime%WaterLevel
            TIWaterFluxX      => Me%TimeIntegration%Integrate3DTime%WaterFluxX
            TIWaterFluxY      => Me%TimeIntegration%Integrate3DTime%WaterFluxY
            TIDischarges      => Me%TimeIntegration%Integrate3DTime%Discharges
            TIComputeFacesU3D => Me%TimeIntegration%Integrate3DTime%ComputeFacesU3D
            TIComputeFacesV3D => Me%TimeIntegration%Integrate3DTime%ComputeFacesV3D

        endif

        Vector_X          => Me%SpaceIntegration%VectorIntegrationX
        Vector_Y          => Me%SpaceIntegration%VectorIntegrationY

        NewILB            =  Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%ILB
        NewIUB            =  Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%IUB

        NewJLB            =  Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%JLB
        NewJUB            =  Me%SpaceIntegration%Integrate3DSpace%WindowWorkSize%JUB


        KUB               = Me%TimeIntegration%Integrate3DTime%WindowWorkSize%KUB

        !-------------------------------------------------------------------------
        SIWaterLevel      = 0.0
        SIWaterFluxX      = 0.0
        SIWaterFluxY      = 0.0
        SIDischarges      = 0.0
        SIComputeFacesU3D = 0
        SIComputeFacesV3D = 0
        SIWaterPoints2D   = 0

            
        !Get Horizontal Areas-------------------------------------------------------
        call GetHorizontalGrid(Me%ObjHorizontalGrid, DUX = DUX, DVY = DVY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'Subroutine IntegrateInSpace; module ModuleHydrodynamicFile. ERR01.'
                 
        !Get OpenPoints-----------------------------------------------------------
        if (Me%State%TimeIntegration .and. (.not. present(InitialFields))) then
            
            call GetHydroIntegrationComputeFaces(                                                      &
                               Me%ObjHydroIntegration,                                &
                               Me%TimeIntegration%Integrate3DTime%IntegrationStep,    &
                               Me%TimeIntegration%Integrate3DTime%WindowWorkSize,     &
                               OpenPoints3D = OpenPoints,                                              &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                           &
            stop 'Subroutine IntegrateInSpace; Module ModuleHydrodynamicFile. ERR02'
        
        else
        
            call GetOpenPoints3D(Me%ObjMap,                                           &
                                 OpenPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                                 &
            stop 'Subroutine IntegrateInSpace; Module ModuleHydrodynamicFile. ERR03'

        endif

        !Get WaterPoints----------------------------------------------------------
        call GetWaterPoints3D(Me%ObjMap, WaterPoints3D,          &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                            &
            Stop 'Subroutine IntegrateInSpace; module ModuleHydrodynamicFile. ERR00a.'

        !WaterFluxX---------------------------------------------------------------
        do J = NewJLB, NewJUB
            do I = NewILB, NewIUB
                Real8Aux = 0.0
                do M = Vector_X(J), Vector_X(J)
                    do L = Vector_Y(I), Vector_Y(I + 1) - 1
                        
                        if(TIComputeFacesU3D(L, M, KUB) == 1) &
                        Real8Aux =  Real8Aux + TIWaterFluxX(L, M, KUB)

                    enddo
                end do
                SIWaterFluxX(I, J, KUB) = Real8Aux
            end do
        end do

        !Me%SpaceIntegration%Integrate3DSpace%WaterFluxX = SIWaterFluxX

        !WaterFluxY---------------------------------------------------------------
        do J = NewJLB, NewJUB
            do I = NewILB, NewIUB
                Real8Aux = 0.0
                do M = Vector_X(J), Vector_X(J + 1) - 1
                    do L = Vector_Y(I), Vector_Y(I)
                        
                        if(TIComputeFacesV3D(L, M, KUB) == 1) &
                        Real8Aux = Real8Aux + TIWaterFluxY(L, M, KUB)

                    enddo
                end do
                SIWaterFluxY(I, J, KUB) = Real8Aux
            end do
        end do

        !Me%SpaceIntegration%Integrate3DSpace%WaterFluxY = SIWaterFluxY

        !ComputeFacesU3D----------------------------------------------------------
        do J = NewJLB, NewJUB
            do I = NewILB, NewIUB
                SumComputeFacesU3D = sum(TIComputeFacesU3D(Vector_Y(I):Vector_Y(I + 1) - 1,    &
                                                           Vector_X(J):Vector_X(J    ),        &
                                                           KUB)                                )

                if(SumComputeFacesU3D > 0)then
                    SIComputeFacesU3D(I, J, KUB) = 1
                else
                    SIComputeFacesU3D(I, J, KUB) = 0
                endif

            end do
        end do

        !Me%SpaceIntegration%Integrate3DSpace%ComputeFacesU3D = SIComputeFacesU3D
        
        !ComputeFacesV3D----------------------------------------------------------
        do J = NewJLB, NewJUB
            do I = NewILB, NewIUB
                SumComputeFacesV3D = sum(TIComputeFacesV3D(Vector_Y(I):Vector_Y(I    ),        &
                                                           Vector_X(J):Vector_X(J + 1) - 1,    &
                                                           KUB)                                )

                if(SumComputeFacesV3D > 0)then
                    SIComputeFacesV3D(I, J, KUB) = 1
                else
                    SIComputeFacesV3D(I, J, KUB) = 0
                endif

            end do
        end do

        !Me%SpaceIntegration%Integrate3DSpace%ComputeFacesV3D = SIComputeFacesV3D

        !Discharges---------------------------------------------------------------
        do J = NewJLB, NewJUB
            do I = NewILB, NewIUB
                Real8Aux = 0.0
                do M = Vector_X(J), Vector_X(J + 1) - 1
                    do L = Vector_Y(I), Vector_Y(I + 1) - 1
                        
                        if(OpenPoints(L, M, KUB) == 1) &
                        Real8Aux = Real8Aux + TIDischarges(L, M, KUB)

                    enddo
                end do
                SIDischarges(I, J, KUB) = Real8Aux
            end do
        end do

        !Me%SpaceIntegration%Integrate3DSpace%Discharges = SIDischarges

        !WaterLevel---------------------------------------------------------------
        do J = NewJLB, NewJUB
            do I = NewILB, NewIUB

                IntegerAux = 0
                do M = Vector_X(J), Vector_X(J + 1) - 1
                    do L = Vector_Y(I), Vector_Y(I + 1) - 1
                        if(WaterPoints3D(L, M, KUB) == 1) IntegerAux = 1
                    end do
                end do
                
                if(IntegerAux == 1)then

                IntegratedCellArea = 0.0
                RealAux    = 0.0
                do M = Vector_X(J), Vector_X(J + 1) - 1
                    do L = Vector_Y(I), Vector_Y(I + 1) - 1
                        
                        CellArea = DUX(L, M) * DVY(L, M)
                        if(WaterPoints3D(L, M, KUB) == 1)then
                            RealAux = RealAux +                         &
                                                TIWaterLevel(L, M   ) * &
                                                CellArea
                        else
                            RealAux = RealAux -                                                           &
                                                Me%SpaceIntegration%BatMin(I, J) * &
                                                CellArea

                        endif
                        IntegratedCellArea = IntegratedCellArea + CellArea

                    enddo
                end do
                
                SIWaterLevel(I, J) = RealAux / IntegratedCellArea
                
                endif

            end do
        end do

        !Me%SpaceIntegration%Integrate3DSpace%WaterLevel = SIWaterLevel

        !WaterPoints2D---------------------------------------------------------
        do J = NewJLB, NewJUB
            do I = NewILB, NewIUB

                IntegerAux = 0
                do M = Vector_X(J), Vector_X(J + 1) - 1
                    do L = Vector_Y(I), Vector_Y(I + 1) - 1
                        if(WaterPoints3D(L, M, KUB) == 1) IntegerAux = 1
                    end do
                end do
                if(IntegerAux == 1)then
                    SIWaterPoints2D(i, j) = 1
                endif
                    
            end do
        end do


        !Unget Horizontal Areas----------------------------------------------------
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,               &
                                 DUX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                           &
                stop 'Subroutine IntegrateInSpace; Module ModuleHydrodynamicFile. ERR04' 

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,               &
                                 DVY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                           &
                stop 'Subroutine IntegrateInSpace; Module ModuleHydrodynamicFile. ERR05' 

        !Unget OpenPoints---------------------------------------------------------
        if (Me%State%TimeIntegration .and. (.not. present(InitialFields))) then

            call UnGetHydroIntegration(                                              &
                        Me%ObjHydroIntegration,                     &
                        OpenPoints,                                                  &
                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                           &
                    stop 'Subroutine IntegrateInSpace; Module ModuleHydrodynamicFile. ERR06'

        else

            call UnGetMap(Me%ObjMap,                                 &
                          OpenPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                &
                stop 'Subroutine IntegrateInSpace; Module ModuleHydrodynamicFile. ERR07'

        endif

        !Unget WaterPoints---------------------------------------------------------
        call UnGetMap(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            Stop 'Subroutine IntegrateInSpace; Module ModuleHydrodynamicFile. ERR08'


        !Nullify auxiliar pointers------------------------------------------------
        nullify(SIWaterLevel)
        nullify(SIWaterFluxX)
        nullify(SIWaterFluxY)
        nullify(SIDischarges)
        nullify(SIComputeFacesU3D)
        nullify(SIComputeFacesV3D)

        nullify(TIWaterLevel)
        nullify(TIWaterFluxX)
        nullify(TIWaterFluxY)
        nullify(TIDischarges)
        nullify(TIComputeFacesU3D)
        nullify(TIComputeFacesV3D)

        nullify(Vector_X)
        nullify(Vector_Y)
        
        !-------------------------------------------------------------------------

    end subroutine IntegrateInSpace

    !----------------------------------------------------------------------------

    subroutine WriteBlock(IntegratedFields)

        !Arguments---------------------------------------------------------------
        type(T_STIntegration),    intent(IN)        :: IntegratedFields

        !Local-------------------------------------------------------------------
        real                                        :: Year
        real                                        :: Month
        real                                        :: Day
        real                                        :: Hour
        real                                        :: Minute
        real                                        :: Second
        integer                                     :: ILB, IUB
        integer                                     :: JLB, JUB
        integer                                     :: KLB, KUB
        integer                                     :: i, j, k, STAT_CALL
        integer                                     :: Unit
        integer                                     :: VersionNumber
        integer, dimension(:, :), pointer           :: WaterPoints2D

        !------------------------------------------------------------------------                         

        ILB = IntegratedFields%WindowWorkSize%ILB
        IUB = IntegratedFields%WindowWorkSize%IUB
        
        JLB = IntegratedFields%WindowWorkSize%JLB
        JUB = IntegratedFields%WindowWorkSize%JUB
        
        KLB = IntegratedFields%WindowWorkSize%KLB
        KUB = IntegratedFields%WindowWorkSize%KUB

        Unit            = Me%Files%OutPutFieldsUnit
        VersionNumber   = Me%Files%OutputVersionNumber

        if (Me%State%SpaceIntegration) then
            WaterPoints2D => Me%SpaceIntegration%Integrate3DSpace%WaterPoints2D
        else
            call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D,   &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteBlock - ModuleHydrodynamicFile - ERR01'
        endif
            

        !Gets the current time
        call ExtractDate(Me%ExternalVar%Now,                                    &
                         Year   = Year,                                         &
                         Month  = Month,                                        &
                         Day    = Day,                                          &
                         Hour   = Hour,                                         &
                         Minute = Minute,                                       &
                         Second = Second)

        !Writes Begin block
        write(Unit)BeginBlock

        !Writes the current time
        write(Unit) Year, Month, Day, Hour, Minute, Second

        !Writes the Waterlevel
        do j = JLB, JUB
        do i = ILB, IUB
        if (WaterPoints2D(i, j) == WaterPoint .or. VersionNumber == 1) then
            write(Unit)IntegratedFields%WaterLevel(i, j)
        endif
        enddo
        enddo

        !Writes WaterFlux X
        do k = KLB, KUB
        do j = JLB, JUB + 1
        do i = ILB, IUB
        if (WaterPoints2D(i, min(j, JUB)) == WaterPoint .or. VersionNumber == 1) then
            write(Unit)IntegratedFields%WaterFluxX(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Writes WaterFlux Y
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB + 1
        if (WaterPoints2D(min(i, IUB), j) == WaterPoint .or. VersionNumber == 1) then
            write(Unit)IntegratedFields%WaterFluxY(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Writes Discharges
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
        if (WaterPoints2D(i, j) == WaterPoint .or. VersionNumber == 1) then
            write(Unit)IntegratedFields%Discharges(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Writes ComputeFacesU3D
        do k = KLB, KUB
        do j = JLB, JUB + 1
        do i = ILB, IUB
        if (WaterPoints2D(i, min(j, JUB)) == WaterPoint .or. VersionNumber == 1) then
            write(Unit)IntegratedFields%ComputeFacesU3D(i, j, k)
        endif
        enddo
        enddo
        enddo

        !Writes ComputeFacesV3D
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB + 1
        if (WaterPoints2D(min(i, IUB), j) == WaterPoint .or. VersionNumber == 1) then
            write(Unit)IntegratedFields%ComputeFacesV3D(i, j, k)
        endif
        enddo
        enddo
        enddo


        write(Unit)EndBlock


        if (Me%State%SpaceIntegration) then
            nullify (WaterPoints2D)
        else
            call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D,   &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteBlock - ModuleHydrodynamicFile - ERR02'
        endif



        !------------------------------------------------------------------------

    end subroutine WriteBlock


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillHydrodynamicFile(HydrodynamicFileID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HydrodynamicFileID 
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_
        integer                                     :: nUsers, STAT_CALL
        logical                                     :: opened

        !------------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(HydrodynamicFileID, ready_)    

        if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mHYDRODYNAMICFILE_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%ObjHydroIntegration /= 0) then
                    call KillHydroIntegration(Me%ObjHydroIntegration, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR01'
                endif


                !Deassociates External Instances
                nUsers = DeassociateInstance (mGRIDDATA_,     Me%ObjTopography    )
                if (nUsers == 0) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR02'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap )
                if (nUsers == 0) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR03'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR04'

                nUsers = DeassociateInstance (mGEOMETRY_,       Me%ObjGeometry      )
                if (nUsers == 0) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR05'

                nUsers = DeassociateInstance (mMAP_,            Me%ObjMap           )
                if (nUsers == 0) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR06'

                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime          )
                if (nUsers == 0) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR07'

        
                if (Me%State%INPUT) call DeallocateVariablesInput 

                if (Me%State%OUTPUT) then

                    inquire(UNIT   = Me%Files%OutPutFieldsUnit, OPENED = opened,         &
                            IOSTAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR08'

cd3 :               if (opened) then
                        call UnitsManager(Me%Files%OutPutFieldsUnit, CLOSE_FILE, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillHydrodynamicFile - ModuleHydrodynamicFile - ERR09'
                    end if cd3
                    
                endif

                call DeallocateInstance

                HydrodynamicFileID = 0
                STAT_              = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if



        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillHydrodynamicFile

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HydrodynamicFile), pointer          :: AuxObjHydrodynamicFile
        type (T_HydrodynamicFile), pointer          :: PreviousObjHydrodynamicFile

        !Updates pointers
        if (Me%InstanceID == FirstHydrodynamicFile%InstanceID) then
            FirstHydrodynamicFile => FirstHydrodynamicFile%Next
        else
            PreviousObjHydrodynamicFile => FirstHydrodynamicFile
            AuxObjHydrodynamicFile      => FirstHydrodynamicFile%Next
            do while (AuxObjHydrodynamicFile%InstanceID /= Me%InstanceID)
                PreviousObjHydrodynamicFile => AuxObjHydrodynamicFile
                AuxObjHydrodynamicFile      => AuxObjHydrodynamicFile%Next
            enddo

            !Now update linked list
            PreviousObjHydrodynamicFile%Next => AuxObjHydrodynamicFile%Next

        endif
            
        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

    end subroutine DeallocateInstance

    !----------------------------------------------------------------------------

    subroutine DeallocateVariablesInput

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        deallocate(Me%Input%WaterLevel_New) 
        deallocate(Me%Input%InitialWaterLevel) 
        deallocate(Me%Input%WaterFluxX) 
        deallocate(Me%Input%WaterFluxY) 
        deallocate(Me%Input%Discharges) 
        deallocate(Me%Input%ComputeFacesU3D) 
        deallocate(Me%Input%ComputeFacesV3D) 
        deallocate(Me%Input%InitialWaterFluxX) 
        deallocate(Me%Input%InitialWaterFluxY) 
        deallocate(Me%Input%InitialDischarges) 
        deallocate(Me%Input%InitialComputeFacesU3D) 
        deallocate(Me%Input%InitialComputeFacesV3D) 
        deallocate(Me%Input%WaterPoints2D)

        nullify(Me%Input%InitialWaterLevel        )
        nullify(Me%Input%InitialWaterFluxX        )
        nullify(Me%Input%InitialWaterFluxY        )
        nullify(Me%Input%InitialDischarges        )
        nullify(Me%Input%InitialComputeFacesU3D   )
        nullify(Me%Input%InitialComputeFacesV3D   )

        nullify(Me%Input%WaterLevel_New    )
        nullify(Me%Input%WaterFluxX        )
        nullify(Me%Input%WaterFluxY        )
        nullify(Me%Input%Discharges        )
        nullify(Me%Input%ComputeFacesU3D   )
        nullify(Me%Input%ComputeFacesV3D   )

        !----------------------------------------------------------------------

    end subroutine DeallocateVariablesInput

    !--------------------------------------------------------------------------    

    subroutine DeallocateVariablesOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------


        deallocate(Me%TimeIntegration%Integrate3DTime%WaterLevel)
        deallocate(Me%TimeIntegration%Integrate3DTime%WaterFluxX)
        deallocate(Me%TimeIntegration%Integrate3DTime%WaterFluxY)
        deallocate(Me%TimeIntegration%Integrate3DTime%Discharges)
        deallocate(Me%TimeIntegration%Integrate3DTime%ComputeFacesU3D)
        deallocate(Me%TimeIntegration%Integrate3DTime%ComputeFacesV3D)


        nullify(Me%TimeIntegration%Integrate3DTime%WaterLevel     )
        nullify(Me%TimeIntegration%Integrate3DTime%WaterFluxX     )
        nullify(Me%TimeIntegration%Integrate3DTime%WaterFluxY     )
        nullify(Me%TimeIntegration%Integrate3DTime%Discharges     )
        nullify(Me%TimeIntegration%Integrate3DTime%ComputeFacesU3D)
        nullify(Me%TimeIntegration%Integrate3DTime%ComputeFacesV3D)

        if (Me%State%SpaceIntegration) then

            deallocate(Me%SpaceIntegration%VectorIntegrationY)
            deallocate(Me%SpaceIntegration%VectorIntegrationX)
            deallocate(Me%SpaceIntegration%BatMin)
            deallocate(Me%SpaceIntegration%Integrate3DSpace%WaterLevel)
            deallocate(Me%SpaceIntegration%Integrate3DSpace%WaterFluxX)
            deallocate(Me%SpaceIntegration%Integrate3DSpace%WaterFluxY)
            deallocate(Me%SpaceIntegration%Integrate3DSpace%Discharges)
            deallocate(Me%SpaceIntegration%Integrate3DSpace%ComputeFacesU3D)
            deallocate(Me%SpaceIntegration%Integrate3DSpace%ComputeFacesV3D)
            deallocate(Me%SpaceIntegration%Integrate3DSpace%WaterPoints2D  )

            nullify(Me%SpaceIntegration%VectorIntegrationY              )
            nullify(Me%SpaceIntegration%VectorIntegrationX              )
            nullify(Me%SpaceIntegration%BatMin                          )
            nullify(Me%SpaceIntegration%Integrate3DSpace%WaterLevel     )
            nullify(Me%SpaceIntegration%Integrate3DSpace%WaterFluxX     )
            nullify(Me%SpaceIntegration%Integrate3DSpace%WaterFluxY     )
            nullify(Me%SpaceIntegration%Integrate3DSpace%Discharges     )
            nullify(Me%SpaceIntegration%Integrate3DSpace%ComputeFacesU3D)
            nullify(Me%SpaceIntegration%Integrate3DSpace%ComputeFacesV3D)
            nullify(Me%SpaceIntegration%Integrate3DSpace%WaterPoints2D  )

        endif


    end subroutine DeallocateVariablesOutput


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (ObjHydrodynamicFile_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHydrodynamicFile_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjHydrodynamicFile_ID > 0) then
            call LocateObjHydrodynamicFile(ObjHydrodynamicFile_ID)
            ready_ = VerifyReadLock (mHYDRODYNAMICFILE_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjHydrodynamicFile (ObjHydrodynamicFileID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHydrodynamicFileID

        !Local-----------------------------------------------------------------

        Me => FirstHydrodynamicFile
        do while (associated (Me))
            if (Me%InstanceID == ObjHydrodynamicFileID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'HydrodynamicFile - LocateObjHydrodynamicFile - ERR01'

    end subroutine LocateObjHydrodynamicFile

    !--------------------------------------------------------------------------

end Module ModuleHydrodynamicFile 

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

