!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Turbulence 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2003
! REVISION      : Pedro Galvao / Frank Braunschweig - v4.0
! DESCRIPTION   : Module which calculates the horizontal and vertical eddy 
!                 viscosities and diffusivities from different turbulent closures 
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
!MODULE TURBULENCE
!
! Computes horizontal and vertical eddy viscosities and diffusivities from different turbulent closures 
!0
! EXISTING GETS:
!    public  :: GetHorizontalViscosity
!    public  :: GetVerticalViscosity
!    public  :: GetVerticalDiffusivity
!
! !HISTORY:
!
! 1999: MOHID2000 Module Turbulence by Ricardo Miranda and Paulo Chambel
! Jan 2000:MRV -Call to new TurbGOTM module added (subroutine MellorYamadaTL removed). 
!               The new keyword to call TurbGOTM model is "turbulence_equation"
! Feb 2000:MRV -Instead of calculating verticalviscosity and verticalprandtlnumber, 
!               verticalviscosity and verticaldiffusivity are calculated. 
!               Now there is a possibility of adding a constant background viscosity
!               to those calculated. By default it is the molecular. It will be possible to add
!               a constant diffusivity in waterproperties Diffusivity=a*TurbulentDiff+b, where a and
!               b could be chosen by the user for every property... 
!               
!              -Brunt-Vaisalla and PRandtl frequency are calculated here and transported by argument to TurbGOTM     
!  
! 2001 : MRV modifications to Richardson, Leendertsee, call to TurbGOTM, output of turbulence magnitudes
! 
!
! !LAST MODIFICATION
! 
! March 2001: Manuel Ruiz Villarreal
!
! !BUGS:
!

Module ModuleTurbulence

    use ModuleEnterData  
    use ModuleFunctions,        only : SigmaUNESCO, SigmaUNESCOPressureCorrection,       &
                                       SigmaMel96PressureCorrection, ConvertTemperature, &
                                       SigmaJMD95PressureCorrection, SigmaLeendertse,    &
                                       SigmaWang, SetMatrixValue, CHUNK_J, CHUNK_K
    use ModuleGlobalData
    use ModuleGridData,         only : GetGridData, UngetGridData   
    use ModuleGeometry,         only : GetGeometrySize, GetGeometryWaterColumn,          &
                                       GetGeometryDistances, GetGeometryKFloor,          &
                                       GetLayer4Level, UnGetGeometry
    use ModuleHDF5
    use ModuleProfile,          only : StartProfile, WriteProfile, KillProfile,          &
                                       GetProfileNextOutputTime
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, UngetHorizontalGrid,           &
                                       WriteHorizontalGrid, GetXYCellZ
    use ModuleHorizontalMap,    only : GetWaterPoints2D, UnGetHorizontalMap
    use ModuleMap,              only : GetWaterPoints3D, GetOpenPoints3D,                &
                                       GetComputeFaces3D, GetImposedTangentialFaces,     &
                                       GetImposedNormalFaces, UnGetMap             
    use ModuleStatistic,        only : ConstructStatistic, GetStatisticMethod,           &
                                       GetStatisticParameters, ModifyStatistic,          &
                                       KillStatistic
    use ModuleStopWatch,        only : StartWatch, StopWatch         
    use ModuleTime                 
    use ModuleTimeSerie,        only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,    &
                                       GetTimeSerieLocation, CorrectsCellsTimeSerie,     &
                                       GetNumberOfTimeSeries, TryIgnoreTimeSerie

    use ModuleTurbGOTM         
    use ModuleFillMatrix,       only : ConstructFillMatrix, GetIfMatrixRemainsConstant,  &
                                       GetDefaultValue, KillFillMatrix

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructTurbulence
    private  ::     AllocateInstance
    private ::      AllocateVariables
    private ::          ConstructMLDStatistics
    private ::      TurbulenceOptions
    private ::          CheckVerticalTurbulenceOption
    private ::          CheckHorizontalTurbulenceOption
    private ::      InicHorizontalModels
    private ::          InicConstantHorizontalModel
    private ::          InicEstuarySmagorinskyModel
    private ::          InicSmagorinskyModel
    private ::          InicEstuaryModel
    private ::          InicFileHorizontalModel
    private ::      InicVerticalModels
    private ::          InicConstantVerticalModel
    private ::          InicNihoulLeendertseeModel
    private ::      TurbulentViscosity_CellCorner
    private ::      InicMixingLengthHorizontal
    private ::      Construct_Time_Serie
    private ::      Construct_Output_Profile
    private ::          Open_HDF5_OutPut_File

    !Selector
    public  :: GetHorizontalViscosity
    public  :: GetVerticalViscosity
    public  :: GetVerticalDiffusivity
    public  :: GetMixingLengthVertical              !rcm
    public  :: GetMixingLengthHorizontal            !rcm
    public  :: GetMLD_Surf
    public  :: GetContinuousGOTM
    public  :: GetTurbulenceOptions

    public  :: SetTurbulenceBottomRugosity

    public  :: UngetTurbulence
                     
    
    !Modifier
    public  :: Turbulence
    private ::      LeendertseeModel                !Vertical
    private ::      BackhausModel                   !Vertical
    private ::      PacanowskiModel                 !Vertical
    private ::      NihoulModel                     !Vertical
    private ::      TurbulenceEquationModel         !Vertical  
    private ::      Richardson                      !Vertical
    private ::      ComputeMixedLayerDepth          !Vertical
    private ::          Output_Statistics
    private ::      EstuaryModel                    !Horizontal
    private ::      SmagorinskyModel                !Horizontal
    private ::      OutPut_TimeSeries
    private ::      Output_Profile
    private ::      OutPut_Results_HDF
    private ::          Write_HDF5_Format

    !Destructor
    public  ::  KillTurbulence

    !Management
    private ::      Ready
    private ::          LocateObjTurbulence

    !Interfaces----------------------------------------------------------------

    private :: UngetTurbulence2D
    private :: UngetTurbulence3D
    interface  UngetTurbulence
        module procedure UngetTurbulence2D
        module procedure UngetTurbulence3D
    end interface UngetTurbulence

    !Parameter-----------------------------------------------------------------

    !Constants                          
    real,    parameter :: CNH           = -0.8           !NIHOUL constant
    real,    parameter :: CVK           =  0.41          !VON KARMAN constant
    real,    parameter :: DELTA         =  0.5           !Used by WALLEFFE
    real,    parameter :: G             =  9.81          !Gravitic constant
    real,    parameter :: RichardsonMAX =  10.0
    real,    parameter :: CBK1          =  7.0           !Backhaus constant 1
    real,    parameter :: CBK2          = -0.25          !Backhaus constant 2
    real,    parameter :: CBK3          =  0.0025        !Backhaus constant 3 !NOT USED IT'S A BACKGROUND VISCOSITY
    real,    parameter :: CBK4          =  0.0002669     !Backhaus constant 4
    real,    parameter :: CBK5          =  1.75          !Backhaus constant 5
    real,    parameter :: CPW1          =  5.0           !Pacanowski constant 1
    real,    parameter :: ICPW2         = -2.0           !Pacanowski constant 2
    
    !It's molecular viscosity    real,    parameter :: CPW3          = 0.00001       !Pacanowski constant 3
    real,    parameter :: CPW4          =  0.00001         !Pacanowski constant 4
    real,    parameter :: GAMA          =  1.4          ! Vertical Diffusivity 
    real,    parameter :: GAMA2         =  1.2          ! Vertical Diffusivity  
    
  
    !Horizontal turbulence options
    integer, parameter :: Constant_        = 1
    integer, parameter :: Estuary_         = 2
    integer, parameter :: Smagorinsky_     = 3
   
    !InicJuan
    integer, parameter :: File2D_          = 10
    !FimJuan

    !Vertical turbulence options
    integer, parameter :: leendertsee_        = 4
    integer, parameter :: backhaus_           = 5
    integer, parameter :: pacanowski_         = 6
    integer, parameter :: nihoul_             = 7
    integer, parameter :: TurbulenceEquation_ = 8
    integer, parameter :: Layers_             = 9

    !3D turbulence options
    integer, parameter :: smagorinsky3d_      = 1
    integer, parameter :: lomax3d_            = 2
    integer, parameter :: smagorinsky3d_vh_   = 3

    !Mixed layer options
    integer, parameter :: tke_mld_        = 1
    integer, parameter :: rich_mld_       = 2
    integer, parameter :: brunt_mld_      = 3

    character(LEN = StringLength), parameter        :: viscosity_h_begin            = '<begin_viscosity_h>'
    character(LEN = StringLength), parameter        :: viscosity_h_end              = '<end_viscosity_h>'
    character(LEN = StringLength), parameter        :: viscosity_v_begin            = '<begin_viscosity_v>'
    character(LEN = StringLength), parameter        :: viscosity_v_end              = '<end_viscosity_v>'

    !Types---------------------------------------------------------------------
    type       T_Files
         character(len=PathLength)                  :: ConstructData
         character(len=PathLength)                  :: OutPutFields
    end type T_Files


    type       T_TurbOptions
         integer                                    :: MODTURB = null_int
         integer                                    :: MODVISH = null_int
         integer                                    :: MLD_Method
         logical                                    :: Continuous_Compute
         logical                                    :: MLD_Calc
         logical                                    :: MLD_Calc_Bot
         integer                                    :: DensityMethod
         logical                                    :: PressureCorrec
    end type T_TurbOptions

    type       T_Statistics
         integer                                    :: ID, FrBV_ID
         character(LEN=StringLength)                :: File 
         logical                                    :: ON
    end type T_Statistics

    type       T_Viscosity
        real, dimension(:,:,:), pointer             :: Vertical
        real, dimension(:,:,:), pointer             :: HorizontalCenter
        real, dimension(:,:,:), pointer             :: HorizontalCorner
        real                                        :: Background
    end type T_Viscosity

    type       T_Diffusivity
        real, dimension(:,:,:), pointer             :: Vertical
    end type T_Diffusivity

    type       T_TurbVar
        real, dimension(:,:,:), pointer             :: Richardson
        real, dimension(:,:,:), pointer             :: FPRANDTL           !Prandtl freq.
        real, dimension(:,:,:), pointer             :: FBRUNTV            !Brunt-Vaisalla freq.
        real, dimension(:,:,:), pointer             :: MixingLengthX      !MixingLength
        real, dimension(:,:,:), pointer             :: MixingLengthY      !MixingLength
        real, dimension(:,:,:), pointer             :: MixingLengthZ      !MixingLength
        real, dimension(:,:,:), pointer             :: Ldownward          !H. Coelho
        real, dimension(:,:,:), pointer             :: VMOD
        real, dimension(:,:,:), pointer             :: VertPrandtlNumber  !Vertical Prandtl Number 
        real, dimension(:,:)  , pointer             :: MLD_Surf
        real, dimension(:,:)  , pointer             :: MLD_Bot
        !Variables from TurbGOTM. Needed for output.
        real, dimension(:,:,:), pointer             :: TKE, L, eps, P, B          
        real                                        :: MAXMixingLength        = null_real     !Nihoul & Leendertsee
        real                                        :: MINHorizontalViscosity = null_real     !Estuary & Smagorinsky
        real                                        :: ReferenceDepth         = null_real     !Estuary
        real                                        :: ReferenceVelocity      = null_real     !Estuary
        real                                        :: HORCON                 = null_real     !Smagorinsky
        real                                        :: TKE_MLD,RICH_MLD
    end type T_TurbVar

    type       T_External
        !Time
        type(T_Time         )                       :: Now
        type(T_Time         )                       :: BeginTime  
        type(T_Time         )                       :: EndTime
        
        !Hydrodynamic
        real,    pointer, dimension(:,:,:)          :: VelocityX
        real,    pointer, dimension(:,:,:)          :: VelocityY
        real,    pointer, dimension(:,:,:)          :: VelocityZ
        real,    pointer, dimension(:,:  )          :: Chezy
        !Geometry
        real,    pointer, dimension(:,:,:)          :: ZCellCenter
        real,    pointer, dimension(:,:,:)          :: SZZ
        real,    pointer, dimension(:,:,:)          :: DWZ
        real,    pointer, dimension(:,:,:)          :: DZZ 
        real,    pointer, dimension(:,:  )          :: HT
        integer, pointer, dimension(:,:  )          :: KFloorZ

        !Bathymetry    
        real,    pointer, dimension(:,:  )          :: Bathymetry        

        !HorizontalGrid
        real,    pointer, dimension(:,:  )          :: DUX
        real,    pointer, dimension(:,:  )          :: DVY
        real,    pointer, dimension(:,:  )          :: DYY
        real,    pointer, dimension(:,:  )          :: DZX
        real,    pointer, dimension(:,:  )          :: DZY

        !WaterProperties
        real,    dimension(:,:,:), pointer          :: Salinity
        real,    dimension(:,:,:), pointer          :: Temperature

        !Bottom
        real                                        :: BottomRugosity

        !Map
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D   
        integer, pointer, dimension(:,:,:)          :: OpenPoints3D
        integer, pointer, dimension(:,:,:)          :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:)          :: ComputeFacesV3D
        integer, pointer, dimension(:,:,:)          :: ComputeFacesW3D
        integer, pointer, dimension(:,:,:)          :: ImposedTangentialFacesU
        integer, pointer, dimension(:,:,:)          :: ImposedTangentialFacesV
        integer, pointer, dimension(:,:,:)          :: ImposedNormalFacesU
        integer, pointer, dimension(:,:,:)          :: ImposedNormalFacesV

        integer, pointer, dimension(:,:  )          :: WaterPoints2D   
    end type T_External

    type       T_OutPut         
         integer                                    :: NextOutPut, Number, NextRestartOutput
         logical                                    :: ON                       = .false.
         logical                                    :: WriteRestartFile         = .false.
         logical                                    :: RestartOverwrite         = .false.
         logical                                    :: Run_End
         type (T_Time), dimension(:), pointer       :: OutTime, RestartOutTime
         real, dimension(:, :, :), pointer          :: Aux3D
         logical                                    :: TimeSerie                = .false.
         logical                                    :: ProfileON                = .false.
    end type T_OutPut

    type      T_Turbulence
        integer                                     :: InstanceID
        character(PathLength)                       :: ModelName
        type(T_Size3D     )                         :: Size   
        type(T_Size3D     )                         :: WorkSize   
        logical                                     :: MixingLengthH = .false.
        type(T_Files      )                         :: Files
        type(T_OutPut     )                         :: OutPut
        type(T_Statistics )                         :: StatMLD
        type(T_TurbOptions)                         :: TurbOptions
        type(T_Viscosity  )                         :: Viscosity
        type(T_Diffusivity)                         :: Diffusivity
        type(T_External   )                         :: ExternalVar
        type(T_TurbVar    )                         :: TurbVar
    

        !Instance of Module_Time
        integer                                     :: ObjTime              = 0
                                                    
        !Instance of Module_EnterData               
        integer                                     :: ObjEnterData         = 0
                                                    
        !Instance of ModuleHorizontalGrid           
        integer                                     :: ObjHorizontalGrid    = 0
                                                    
        !Instance of ModuleGeometry                 
        integer                                     :: ObjGeometry          = 0
                                                    
        !Instance of ModuleMap                      
        integer                                     :: ObjMap               = 0
                                                    
        !Instance of ModuleHorizontalMap            
        integer                                     :: ObjHorizontalMap     = 0    
                                                    
        !Instance of ModuleBathymetry               
        integer                                     :: ObjGridData          = 0
                                                    
        !Instance of ModuleTurbGOTM                 
        integer                                     :: ObjTurbGOTM          = 0
                                                    
        !Instance of ModuleTimeSerie                
        integer                                     :: ObjTimeSerie         = 0
                                                    
        !Instance of ModuleHDF5                     
        integer                                     :: ObjHDF5              = 0

        !Instance of ModuleProfile
        integer                                     :: ObjProfile           = 0

                                                    
        type(T_Turbulence), pointer                 :: Next

    end type T_Turbulence

   
    !Global Module Variables
   
    type (T_Turbulence), pointer                    :: FirstObjTurbulence
    type (T_Turbulence), pointer                    :: Me        

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructTurbulence(ModelName,               &
                                   TurbulenceID,            &
                                   TurbGOTMID,              &
                                   HorizontalGridID,        &
                                   GeometryID,              &
                                   MapID,                   &
                                   GridDataID,              &
                                   HorizontalMapID,         &
                                   TimeID,                  &
                                   STAT)

        !Arguments-------------------------------------------------------------
        character(Len=*)                :: ModelName
        integer                         :: TurbulenceID 
        integer                         :: TurbGOTMID    
        integer                         :: HorizontalGridID 
        integer                         :: GeometryID       
        integer                         :: MapID            
        integer                         :: GridDataID     
        integer                         :: HorizontalMapID  
        integer                         :: TimeID                                    
        integer, optional, intent(OUT)  :: STAT     

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: ready_         

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
 
        !----------------------------------------------------------------------

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTurbulence_)) then
            nullify (FirstObjTurbulence)
            call RegisterModule (mTurbulence_) 
        endif

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ModelName = ModelName

            !Associates External Instances

            Me%ObjTime              = AssociateInstance (mTIME_,           TimeID           )
            Me%ObjHorizontalMap     = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID  )
            Me%ObjHorizontalGrid    = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID )
            Me%ObjGeometry          = AssociateInstance (mGEOMETRY_,       GeometryID       )
            Me%ObjGridData          = AssociateInstance (mGRIDDATA_,       GridDataID       )
            Me%ObjMap               = AssociateInstance (mMAP_,            MapID            )

            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR01'
   
            call GetComputeTimeLimits(  Me%ObjTime,                             & 
                                        BeginTime = Me%ExternalVar%BeginTime,   &
                                        EndTime   = Me%ExternalVar%EndTime,     & 
                                        STAT = STAT_CALL )

            !Reads File Names From Nomfich
            call ReadFileName('IN_TURB', Me%Files%ConstructData, "Turbulence Data File", &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR02'

            call ReadFileName('TURB_HDF', Me%Files%OutPutFields, "Turbulence HDF File",  &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR03'
        

            !Gets the size from the Geometry   
            call GetGeometrySize(Me%ObjGeometry, Size = Me%Size, WorkSize = Me%WorkSize, &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR04'

            !Gets Horizontal Grid
            call GetHorizontalGrid (Me%ObjHorizontalGrid, DUX = Me%ExternalVar%DUX,      &        
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR05'

            call GetHorizontalGrid (Me%ObjHorizontalGrid, DYY = Me%ExternalVar%DYY,      &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR06'

            !Opens Turbulence data file
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR07'

            call GetWaterPoints3D(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR08'


            call TurbulenceOptions
            call AllocateVariables 
            call InicHorizontalModels
            call InicVerticalModels
            call InicMixingLengthHorizontal
            
                        
            !Construct the Time Serie Obj
            if (Me%Output%TimeSerie)   call Construct_Time_Serie

            if (Me%Output%ProfileON)   call Construct_Output_Profile

            if (Me%OutPut%ON) call Open_HDF5_OutPut_File

            if(Me%TurbOptions%MLD_Calc) call ConstructMLDStatistics

            call OutPut_Results_HDF

            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%DUX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR09'


            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%DYY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR10'

            call UngetMap(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR11'

            !Kills Instance of EnterData
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbulence - ModuleTurbulence - ERR12'

            !Returns ID
            TurbGOTMID      = Me%ObjTurbGOTM
            TurbulenceID    = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            
             stop 'ModuleTurbulence - ConstructTurbulence - ERR13'

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructTurbulence

    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Turbulence), pointer                :: NewObjTurbulence
        type (T_Turbulence), pointer                :: PreviousObjTurbulence


        !Allocates new instance
        allocate (NewObjTurbulence)
        nullify  (NewObjTurbulence%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjTurbulence)) then
        
            FirstObjTurbulence      => NewObjTurbulence
            Me                      => NewObjTurbulence
        
        else
        
            PreviousObjTurbulence   => FirstObjTurbulence
            Me                      => FirstObjTurbulence%Next

            do while (associated(Me))
                PreviousObjTurbulence   => Me
                Me                      => Me%Next
            enddo
            
            Me                          => NewObjTurbulence
            PreviousObjTurbulence%Next  => NewObjTurbulence
        
        endif

        Me%InstanceID = RegisterNewInstance (mTURBULENCE_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine Construct_Time_Serie

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber 
        integer                                             :: nProperties,pl
        integer                                             :: STAT_CALL, iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !Begin-----------------------------------------------------------------

        nProperties = 4
        
        if(Me%TurbOptions%MLD_Calc) then
            nProperties = 5
            if(Me%TurbOptions%MLD_Calc_Bot) nProperties = 6
        end if
        
        if (Me%TurbOptions%MODTURB == TurbulenceEquation_) nProperties = nProperties + 5 

        !Allocates PropertyList
        allocate(PropertyList(nProperties))

        !Fills up PropertyList
        PropertyList(1) = 'Visc' !Eddy viscosity
        PropertyList(2) = 'Diff' !Eddy diffusivity
        PropertyList(3) = 'NN'   !Brunt-Vaisalla frequency
        PropertyList(4) = 'SS'   !Prandtl frequency 
        pl              = 4
        
        if (Me%TurbOptions%MLD_Calc) then
        
            pl = pl + 1
            PropertyList(pl) = 'surf mld'       !Mixed layer depth 

            if (Me%TurbOptions%MLD_Calc_Bot) then 
                pl = pl + 1                
                PropertyList(pl) = 'bottom mld' !Bottom Mixed layer depth
            end if

        end if
        
        if (Me%TurbOptions%MODTURB == TurbulenceEquation_) then
            PropertyList(pl+1) = 'TKE' !Turbulent kinetic energy
            PropertyList(pl+2) = 'eps' !Dissipation Rate of TKE
            PropertyList(pl+3) = 'L'   !Macro length scale
            PropertyList(pl+4) = 'P'   !Production term in turbulence equations
            PropertyList(pl+5) = 'B'   !Buoyancy term in turbulence equations
            pl = pl + 5
        end if 
        
        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleWaterTurbulence',                            &
                     Default      = Me%Files%ConstructData,                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop "Construct_Time_Serie - Turbulence - ERR10" 

        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "srt",                                        &
                            WaterPoints3D = Me%ExternalVar%WaterPoints3D,               &
                            ModelName     = Me%ModelName,                               &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Construct_Time_Serie - Turbulence - ERR20" 

        !Deallocates PropertyList
        deallocate(PropertyList)

        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Construct_Time_Serie - Turbulence - ERR30" 

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop "Construct_Time_Serie - Turbulence - ERR40" 

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop "Construct_Time_Serie - Turbulence - ERR50" 

                    if (IgnoreOK) then
                        cycle
                    else
                        stop "Construct_Time_Serie - Turbulence - ERR60" 
                    endif

                endif

                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop "Construct_Time_Serie - Turbulence - ERR70" 
            endif

        enddo


        
    end subroutine Construct_Time_Serie

    !--------------------------------------------------------------------------

    subroutine Construct_Output_Profile

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL, iflag
        character(len=StringLength)                         :: TimeSerieLocationFile
        integer                                             :: nProperties
        character(len=StringLength), dimension(:,:), pointer:: PropertyList

        !----------------------------------------------------------------------
        
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_Output_Profile - ModuleTurbulence - ERR01' 

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleTurbulence',                                 &
                     Default      = Me%Files%ConstructData,                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_Output_Profile - ModuleTurbulence - ERR02' 
        
        nProperties = 0

        if(Me%TurbOptions%MODTURB .ne. Constant_ .or. &
           Me%TurbOptions%MODTURB .ne. File2D_)then
            
            nProperties = 5
            
            if (Me%TurbOptions%MODTURB == TurbulenceEquation_)then
                nProperties = nProperties + 5 
            end if

        end if

        !Allocates PropertyList
        allocate(PropertyList(1:nProperties,1:2))

        PropertyList(1,1) = "ViscosityH"
        PropertyList(2,1) = "ViscosityZ"
        PropertyList(3,1) = "Diffusivity"
        PropertyList(4,1) = "NN"
        PropertyList(5,1) = "SS"

        PropertyList(1,2) = "m2/s"
        PropertyList(2,2) = "log10(m2/s)"
        PropertyList(3,2) = "log10(m2/s)"
        PropertyList(4,2) = "s-2"
        PropertyList(5,2) = "s-2"

        
        if (Me%TurbOptions%MODTURB == TurbulenceEquation_)then
            PropertyList(6, 1) = "TKE"
            PropertyList(7, 1) = "eps"
            PropertyList(8, 1) = "L"
            PropertyList(9, 1) = "P"
            PropertyList(10,1) = "B"
            PropertyList(6, 2) = "log10(m2/s2)"
            PropertyList(7, 2) = "log10(m2/s3)"
            PropertyList(8, 2) = "log10(m)"
            PropertyList(9, 2) = "m2/s3"
            PropertyList(10,2) = "m2/s3"
        end if

        call StartProfile  (ProfileID       = Me%ObjProfile,                            &
                            ObjTime         = Me%ObjTime,                               &
                            ProfileDataFile = trim(TimeSerieLocationFile),              &
                            WaterPoints2D   = Me%ExternalVar%WaterPoints2D,             &
                            KUB             = Me%WorkSize%KUB,                          &
                            nProperties     = nProperties,                              &
                            PropertyList    = PropertyList,                             &
                            ClientName      = "Turbulence",                             &
                            STAT            = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_Output_Profile - ModuleTurbulence - ERR03' 
       
        deallocate(PropertyList)

        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_Output_Profile - ModuleTurbulence - ERR04' 

    end subroutine Construct_Output_Profile

    !--------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        real, pointer, dimension    (:, :)          :: Bathymetry
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        !Gets a pointer to Bathymetry
        call GetGridData      (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR01'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5,                                    &
                                 trim(Me%Files%OutPutFields)//"5",              &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR03'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5,              &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR04'
        
        !Sets Limits
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,          &
                              WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR07'


        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Bathymetry,                             &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR05'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",        &
                              Array3D = Me%ExternalVar%WaterPoints3D,           &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR06'

        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB+1, WorkJLB,          &
                              WorkJUB+1, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR07'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR08'


        !Ungets the Bathymetry
        call UngetGridData (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleTurbulence - ERR09'

    end subroutine Open_HDF5_OutPut_File
    
    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Local-----------------------------------------------------------------
        integer :: ILB, IUB
        integer :: JLB, JUB
        integer :: KLB, KUB

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB


        !Allocates
        allocate(Me%Viscosity%Vertical        (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Diffusivity%Vertical      (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Viscosity%HorizontalCenter(ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Viscosity%HorizontalCorner(ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%TurbVar%Richardson        (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%TurbVar%FPRANDTL          (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%TurbVar%FBRUNTV           (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%TurbVar%MixingLengthX     (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%TurbVar%MixingLengthY     (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%TurbVar%MixingLengthZ     (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%TurbVar%Ldownward         (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%TurbVar%VertPrandtlNumber (ILB:IUB, JLB:JUB, KLB:KUB))
       

        Me%Viscosity%Vertical           (:,:,:)   = null_real
        Me%Diffusivity%Vertical         (:,:,:)   = null_real
        Me%Viscosity%HorizontalCenter   (:,:,:)   = null_real
        Me%Viscosity%HorizontalCorner   (:,:,:)   = null_real
        Me%TurbVar%Richardson           (:,:,:)   = null_real
        Me%TurbVar%FPRANDTL             (:,:,:)   = null_real
        Me%TurbVar%FBRUNTV              (:,:,:)   = null_real       
        Me%TurbVar%MixingLengthX        (:,:,:)   = null_real
        Me%TurbVar%MixingLengthY        (:,:,:)   = null_real
        Me%TurbVar%MixingLengthZ        (:,:,:)   = null_real
        Me%TurbVar%Ldownward            (:,:,:)   = null_real
        Me%TurbVar%VertPrandtlNumber    (:,:,:)   = null_real


        if (Me%TurbOptions%MODTURB == backhaus_) then 
            allocate(Me%TurbVar%VMOD (ILB:IUB, JLB:JUB, KLB:KUB))
            Me%TurbVar%VMOD(:,:,:) = null_real
        end if


        if(Me%TurbOptions%MLD_Calc) then
            allocate(Me%TurbVar%MLD_Surf (ILB:IUB, JLB:JUB))
            Me%TurbVar%MLD_surf = null_real

            if(Me%TurbOptions%MLD_Calc_Bot) then
                allocate(Me%TurbVar%MLD_bot (ILB:IUB, JLB:JUB))
                Me%TurbVar%MLD_Bot(:,:) = null_real
            end if
        end if


        if(Me%OutPut%ON .or. Me%OutPut%ProfileON)then
            allocate(Me%OutPut%Aux3D(ILB:IUB, JLB:JUB, KLB:KUB))
            Me%OutPut%Aux3D(:,:,:) = null_real
        end if
        

        !----------------------------------------------------------------------

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    
!-------------------------------------------------------------------------
!        IST/MARETEC, Marine Modelling Group, Mohid2000 modelling system
!-------------------------------------------------------------------------

!BOP
!
! !ROUTINE: ConstructMLDStatistics

! !INTERFACE:

    subroutine ConstructMLDStatistics

! !INPUT PARAMETERS:
        
! !DESCRIPTION: 
!   This subroutine reads all the information needed to construct the output          
!   of each property.

 
! !REVISION HISTORY: 
!  20Nov2000   Paulo Chambel  the ConstructMLDStatistics routine calls from
!                             the ModuleStatistic the routine ConstructStatistic           

           
!EOP

!       Local ------------------------------------------------------------------------------
        integer                     :: STAT_CALL
        integer                     :: iflag
        !Begin----------------------------------------------------------------------------

        !<BeginKeyword>
            !Keyword          : STATISTICS_MLD
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends the statistics of Turbulence mixing length
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : IN_TURB                 
            !Multiple Options : Do not have
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%StatMLD%ON,                         &
                     Me%ObjEnterData, iflag,                &
                     Keyword    = 'STATISTICS_MLD',         &
                     Default    = .false.,                  &
                     SearchType = FromFile,                 &
                     ClientModule = 'ModuleTurbulence',     &
                     STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                             &
            call SetError(FATAL_, INTERNAL_, 'ConstructMLDStatistics - Turbulence - ERR01') 

        
cd1:    if (Me%StatMLD%ON) then

            !<BeginKeyword>
                !Keyword          : STATISTICS_MLD_FILE
                !<BeginDescription>       
                   ! 
                   ! The statistics definition file of this property
                   ! 
                !<EndDescription>
                !Type             : Character
                !Default          : Do not have
                !File keyword     : IN_TURB                 
                !Multiple Options : Do not have
                !Search Type      : FromFile
            !<EndKeyword>
            call GetData(   Me%StatMLD%File,                        &
                            Me%ObjEnterData, iflag,                &
                            Keyword    = 'STATISTICS_MLD_FILE',     &
                            SearchType = FromFile,                  &
                            ClientModule = 'ModuleTurbulence',      &
                            STAT = STAT_CALL)

            if (iflag /= 1)                                         &
                call SetError(FATAL_, INTERNAL_, 'ConstructMLDStatistics - Turbulence - ERR02') 

            if (STAT_CALL /= SUCCESS_)                                 &
                call SetError(FATAL_, INTERNAL_, 'ConstructMLDStatistics - Turbulence - ERR03') 

            call ConstructStatistic (Me%StatMLD%ID,                 &
                                     Me%ObjTime,                    &
                                     Me%ObjHDF5,                    &
                                     Me%Size,                       &
                                     Me%WorkSize,                   &
                                     Me%StatMLD%File,               &
                                    'MLD_Surf',                     &
                                     STAT = STAT_CALL)                                 
            if (STAT_CALL /= SUCCESS_)                                 &
                call SetError(FATAL_, INTERNAL_, 'ConstructMLDStatistics - Turbulence - ERR04') 

            if (Me%TurbOptions%MLD_Method == brunt_mld_)   then

                call ConstructStatistic (Me%StatMLD%FrBV_ID,        &
                                         Me%ObjTime,                &
                                         Me%ObjHDF5,                &
                                         Me%Size,                   &
                                         Me%WorkSize,               &
                                         Me%StatMLD%File,           &
                                         'MaxBruntVaisalla',        &
                                         STAT = STAT_CALL)                                 
                if (STAT_CALL /= SUCCESS_)                             &
                    call SetError(FATAL_, INTERNAL_, 'ConstructMLDStatistics - Turbulence - ERR05') 

            endif



        endif cd1


    end subroutine ConstructMLDStatistics

    !--------------------------------------------------------------------------


    subroutine InicHorizontalModels


case1 : select case  (Me%TurbOptions%MODVISH)

            case (File2D_     )

                call InicFileHorizontalModel
                 
            case (Constant_   )

                call InicConstantHorizontalModel
     
            case (Estuary_    )

                call InicEstuarySmagorinskyModel
        
                call InicEstuaryModel           
        
            case (Smagorinsky_)

                call InicEstuarySmagorinskyModel
       
                call InicSmagorinskyModel
        
            case default
                stop 'Subroutine InicHorizontalModels; module ModuleTurbulence. ERR07.'

        end select case1

        
        call TurbulentViscosity_CellCorner


    end subroutine InicHorizontalModels   



    !--------------------------------------------------------------------------

    
    
    subroutine InicVerticalModels

        !External--------------------------------------------------------------

        integer :: STAT_CALL, iflag

        !----------------------------------------------------------------------

case1 : select case  (Me%TurbOptions%MODTURB) 

            case (Constant_             )

                call InicConstantVerticalModel

            case (File2D_     )
                call InicFileVerticalModel

            case (Layers_               )

                call InicLayersVerticalModel
     
            case (leendertsee_          )

                call InicNihoulLeendertseeModel
        
            case (backhaus_             )

                continue
                        
            case (pacanowski_           )

                continue
       
            case (TurbulenceEquation_   )   


                call StartTurbGOTM( Me%ObjTurbGOTM,                     &
                                    Me%ObjTime,                         &
                                    Me%ObjGridData,                     &
                                    Me%ObjMap,                          &
                                    Me%ObjHorizontalMap,                &
                                    Me%ObjHorizontalGrid,               &
                                    Me%ObjGeometry,                     & 
                                    Me%TurbOptions%Continuous_Compute,  &
                                    STAT = STAT_CALL    )
                
                if (STAT_CALL .NE. SUCCESS_)                            &
                   stop 'Subroutine InicVerticalModels; module ModuleTurbulence. ERR03'
        

            case (nihoul_               )

                call InicNihoulLeendertseeModel
        
            case default

                stop 'Subroutine InicVerticalModels; module ModuleTurbulence. ERR05'

            end select case1


            !Reads background viscosity/diffusivity. By default, they are equal to molecular

            !Viscosity
            !By default it is equal to the molecular 
            call GetData(   Me%Viscosity%Background,                    &
                            Me%ObjEnterData, iflag,                     &
                            SearchType   = FromFile,                    &
                            keyword      = 'Background_Viscosity',      &
                            ClientModule = 'ModuleTurbulence',          &
                            default      = 1.3e-6,                      &
                            STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InicVerticalModels - ModuleTurbulence - ERR06'
                 

    end subroutine InicVerticalModels   

    !--------------------------------------------------------------------------

    subroutine InicNihoulLeendertseeModel

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------

        integer :: STAT_CALL
        integer :: flag

        !----------------------------------------------------------------------


        call GetData(   Me%TurbVar%MAXMixingLength,             &
                        Me%ObjEnterData, flag,                  &
                        SearchType   = FromFile,                &
                        ClientModule = 'ModuleTurbulence',      &
                        keyword      = 'MIXLENGTH_MAX',         &
                        Default      = 100.0,                   &
                        STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                          &
            stop 'Subroutine InicNihoulLeendertseeModel; module ModuleTurbulence - ERR01'

cd1 :   if (flag .EQ. 0) then           
            write(*,*) 
            write(*,*) 'Maximum alowed mixing length not defined.'
            write(*,*) 'Keyword - MIXLENGTH_MAX'
            write(*,*) 'value set to : ', Me%TurbVar%MAXMixingLength
            write(*,*) 'Subroutine InicNihoulLeendertseeModel; module ModuleTurbulence - WRN01'
            write(*,*) 
        end if cd1
        

    end subroutine InicNihoulLeendertseeModel   

    !--------------------------------------------------------------------------

    subroutine TurbulenceOptions

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: flag
        character(LEN = StringLength)       :: String
        !----------------------------------------------------------------------

        call GetData(String,                                &
                     Me%ObjEnterData, flag,                 &
                     Keyword      = 'MODTURB',              &
                     ClientModule = 'ModuleTurbulence',     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceOptions - ModuleTurbulence - ERR01'

cd1 :   if (flag .EQ. 0) then
            Me%TurbOptions%MODTURB = Constant_
        else
            String = trim(adjustl(String))
            call CheckVerticalTurbulenceOption(String, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceOptions - ModuleTurbulence - ERR02'
        end if cd1


        call GetData(String,                                &
                     Me%ObjEnterData, flag,                 &
                     Keyword      = 'MODVISH',              &
                     ClientModule = 'ModuleTurbulence',     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceOptions - ModuleTurbulence - ERR03'

cd2 :   if (flag .EQ. 0) then
            Me%TurbOptions%MODVISH = Constant_
        else
            String = trim(adjustl(String))
            call CheckHorizontalTurbulenceOption(String, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceOptions - ModuleTurbulence - ERR04'
        end if cd2


        if (Me%TurbOptions%MODTURB == TurbulenceEquation_) then

             call GetData(  Me%TurbOptions%Continuous_Compute,  &
                            Me%ObjEnterData, flag,              &
                            Keyword      = 'CONTINUOUS',        &
                            ClientModule = 'ModuleTurbulence',  &
                            Default      = .false.,             & 
                            STAT         = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR05'
                

        else

            Me%TurbOptions%Continuous_Compute = .false.

        endif

        call GetData (Me%Output%TimeSerie,                  &
                      Me%ObjEnterData, flag,                &
                      Keyword      = 'TIME_SERIE',          &
                      ClientModule = 'ModuleTurbulence',    &
                      Default      = .false.,               &
                      STAT         = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR06'
        
        
        call GetData (Me%Output%ProfileON,                  &
                      Me%ObjEnterData, flag,                &
                      Keyword      = 'OUTPUT_PROFILE',      &
                      ClientModule = 'ModuleTurbulence',    &
                      Default      = .false.,               &
                      STAT         = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR07'

        !Careful, this getdata is already made in the ModuleWaterProperties 
        call GetData (Me%TurbOptions%DensityMethod,                                      &
                     Me%ObjEnterData, flag,                                              &
                     SearchType = FromFile,                                              &
                     keyword    = 'DENSITY_METHOD',                                      &
                     Default    = UNESCOState_,                                          &
                     ClientModule = 'ModuleWaterProperties',                             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR08'

        call GetData(Me%TurbOptions%PressureCorrec,                                      &
                     Me%ObjEnterData, flag,                                              &
                     SearchType = FromFile,                                              &
                     keyword    = 'PRESSURE_CORRECTION',                                 &
                     Default    = .true.,                                                &
                     ClientModule = 'ModuleWaterProperties',                             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR09'

        !<BeginKeyword>
            !Keyword          : MLD
            !<BeginDescription>       
               ! 
               !Checks out if the user pretends to compute the mixed layer length
               ! 
            !<EndDescription>
            !Type             : Logical
            !Default          : .false.
            !File keyword     : MLD 
            !Multiple Options : .true. , .false.
            !Search Type      : From File
        !<EndKeyword>
        call GetData(Me%TurbOptions%MLD_Calc,               &
                     Me%ObjEnterData, flag,                 &
                     keyword      = 'MLD',                  &
                     ClientModule = 'ModuleTurbulence',     &
                     Default      = .false.,                &
                     STAT         = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR10'
        
        
        if(Me%TurbOptions%MLD_Calc) then

            !<BeginKeyword>
            !Keyword          : MLD_BOTTOM
            !<BeginDescription>       
               ! 
               !Checks out if the user pretends to compute the bottom mixed layer length
               ! 
            !<EndDescription>
            !Type             : Logical
            !Default          : .false.
            !File keyword     : MLD_BOTTOM 
            !Multiple Options : .true. , .false.
            !Search Type      : From File
            !<EndKeyword>
          call GetData(Me%TurbOptions%MLD_Calc_Bot,         &
                     Me%ObjEnterData, flag,                 &
                     keyword      = 'MLD_BOTTOM',           &
                     ClientModule = 'ModuleTurbulence',     &
                     Default      = .false.,                &
                     STAT         = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR11'
           
            !<BeginKeyword>
            !Keyword          : MLD_Method
            !<BeginDescription>       
               ! 
               !Checks out if the user pretends to compute the mixed layer length
               ! 
            !<EndDescription>
            !Type             : Integer
            !Default          : 2
            !File keyword     : MLD_Method
            !Multiple Options : 1 !TKE_min, 2 !Rich_crit 3 !
            !Search Type      : From File
            !<EndKeyword>

           call GetData(Me%TurbOptions%MLD_Method,              &
                        Me%ObjEnterData, flag,                  &
                        keyword      = 'MLD_METHOD',            &
                        ClientModule = 'ModuleTurbulence',      &
                        Default      = rich_mld_,               &
                        STAT         = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR12'

           !<BeginKeyword>
            !Keyword          : TKE_MLD
            !<BeginDescription>       
               ! 
               !Sets the minimum value of TKE for computation of Mixed layer depth
               ! 
            !<EndDescription>
            !Type             : Real
            !Default          : 1e-5
            !File keyword     : TKE_MLD
            !Search Type      : From File
        !<EndKeyword> 
           if(Me%TurbOptions%MLD_Method == TKE_MLD_) then
        
              if (Me%TurbOptions%MODTURB == TurbulenceEquation_) then

                call GetData(Me%TurbVar%TKE_MLD,            &
                     Me%ObjEnterData, flag,                 &
                     keyword      = 'TKE_MLD',              &
                     ClientModule = 'ModuleTurbulence',     &
                     Default      = 1.e-5,                  &
                     STAT         = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR13'
              
              else
              
                write(*,*)
                write(*,*) 'if MODTURB not turbulenceequation, MLD_Method cannot be set to TKE_MLD'
                write(*,*) 'Method set to RICH_MLD'
                write(*,*)  
                Me%TurbOptions%MLD_Method = Rich_MLD_

              end if

           end if

           !<BeginKeyword>
            !Keyword          : RICH_MLD
            !<BeginDescription>       
               ! 
               !Sets the minimum value of RICH for computation of Mixed layer depth
               ! 
            !<EndDescription>
            !Type             : Real
            !Default          : 0.5
            !File keyword     : RICH_MLD
            !Search Type      : From File
           !<EndKeyword> 
        
           if(Me%TurbOptions%MLD_Method == Rich_Mld_) then
        
              call GetData(Me%TurbVar%RICH_MLD,             &
                     Me%ObjEnterData, flag,                 &
                     keyword      = 'RICH_MLD',             &
                     ClientModule = 'ModuleTurbulence',     &
                     Default      = 0.5,                    &
                     STAT         = STAT_CALL)            
             if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR14'
        
           end if

         end if

         
        call GetOutPutTime(Me%ObjEnterData,                                              &
                           CurrentTime = Me%ExternalVar%Now,                             &
                           EndTime     = Me%ExternalVar%EndTime,                         &
                           keyword     = 'OUTPUT_TIME',                                  &
                           SearchType  = FromFile,                                       &
                           OutPutsTime = Me%OutPut%OutTime,                              &
                           OutPutsOn   = Me%OutPut%ON,                                   &
                           STAT        = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR15'

        if (Me%OutPut%ON) then

            Me%OutPut%NextOutPut = 1

        endif 

        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime = Me%ExternalVar%Now,                            &
                           EndTime     = Me%ExternalVar%EndTime,                        &
                           keyword     = 'RESTART_FILE_OUTPUT_TIME',                    &
                           SearchType  = FromFile,                                      &
                           OutPutsTime = Me%OutPut%RestartOutTime,                      &
                           OutPutsOn   = Me%OutPut%WriteRestartFile,                    &
                           STAT        = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)stop 'TurbulenceOptions - ModuleTurbulence - ERR15'

        if(Me%TurbOptions%MODTURB .ne. TurbulenceEquation_ .and. &
           Me%OutPut%WriteRestartFile)then

           write(*,*)'Request for restart file writing during simulation ignored!'
           write(*,*)'in ModuleTurbulence. This is only performed when using GOTM.'
           write(*,*)'Model ID', Me%InstanceID
           write(*,*)'TurbulenceOptions - ModuleTurbulence - WRN01'

           Me%OutPut%WriteRestartFile = .false.

        end if


        if(Me%OutPut%WriteRestartFile)then

            Me%OutPut%NextRestartOutput = 1

        end if 

        !<BeginKeyword>
            !Keyword          : RESTART_FILE_OVERWRITE
            !<BeginDescription>       
               ! 
               ! This option checks wether the restart file is to be overwritten or not
               ! 
            !<EndDescription>
            !Type             : logical 
            !Default          : .true.
            !Multiple Options : Do not have
            !Search Type      : FromFile
        !<EndKeyword>
        
        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'Turbulence',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, KEYWORD_, "TurbulenceOptions - Turbulence - ERR16")

    end subroutine TurbulenceOptions   

    !--------------------------------------------------------------------------

    subroutine CheckVerticalTurbulenceOption(String, STAT)

        !Arguments-------------------------------------------------------------
        character(LEN = *), intent(IN)           :: String
        integer                                  :: STAT

        !Local-----------------------------------------------------------------

        character(LEN = StringLength), parameter :: Char_constant           = trim(adjustl('constant'            ))
        character(LEN = StringLength), parameter :: Char_file2d             = trim(adjustl('file2D'              ))
        character(LEN = StringLength), parameter :: Char_leendertsee        = trim(adjustl('leendertsee'         ))
        character(LEN = StringLength), parameter :: Char_backhaus           = trim(adjustl('backhaus'            ))
        character(LEN = StringLength), parameter :: Char_pacanowski         = trim(adjustl('pacanowski'          ))
        character(LEN = StringLength), parameter :: Char_nihoul             = trim(adjustl('nihoul'              ))
        character(LEN = StringLength), parameter :: Char_turbulenceequation = trim(adjustl('turbulence_equation' ))
                                                                                             
        !----------------------------------------------------------------------

        STAT = UNKNOWN_

case1 : select case(String)

            case(Char_constant       )
                Me%TurbOptions%MODTURB = Constant_
                STAT                              = SUCCESS_
            case(Char_file2d       )
                Me%TurbOptions%MODTURB = file2d_
                STAT                              = SUCCESS_

            case(Char_leendertsee    )
                Me%TurbOptions%MODTURB = leendertsee_
                STAT                              = SUCCESS_

            case(Char_backhaus       )
                Me%TurbOptions%MODTURB = backhaus_
                STAT                              = SUCCESS_

            case(Char_pacanowski     )
                Me%TurbOptions%MODTURB = pacanowski_
                STAT                              = SUCCESS_

            case(Char_turbulenceequation)
                Me%TurbOptions%MODTURB = TurbulenceEquation_
                STAT                              = SUCCESS_

            case(Char_nihoul         )
                Me%TurbOptions%MODTURB = nihoul_
                STAT                              = SUCCESS_

            case default
                Me%TurbOptions%MODTURB = null_int
                STAT                              = NOT_FOUND_ERR_

        end select case1

    end subroutine CheckVerticalTurbulenceOption   


    !--------------------------------------------------------------------------


    subroutine CheckHorizontalTurbulenceOption(String, STAT)

        !Arguments-------------------------------------------------------------
        character(LEN = *), intent(IN)      :: String
        integer                             :: STAT

        !Local-----------------------------------------------------------------
        character(LEN = StringLength), parameter :: Char_constant     = trim(adjustl('constant'     ))
        character(LEN = StringLength), parameter :: Char_Estuary      = trim(adjustl('estuary'      ))
        character(LEN = StringLength), parameter :: Char_smagorinsky  = trim(adjustl('smagorinsky'  ))
        character(LEN = StringLength), parameter :: Char_File2D       = trim(adjustl('file2D'       ))
        character(LEN = StringLength), parameter :: Char_Boxes2D      = trim(adjustl('boxes2D'      ))
       
        !----------------------------------------------------------------------

        STAT = UNKNOWN_

case1 : select case(String)

            case(Char_File2D     )
                Me%TurbOptions%MODVISH = File2D_
                STAT                              = SUCCESS_

            case(Char_constant   )
                Me%TurbOptions%MODVISH = Constant_
                STAT                              = SUCCESS_

            case(Char_Estuary    )
                Me%TurbOptions%MODVISH = Estuary_
                STAT                              = SUCCESS_

            case(Char_smagorinsky)
                Me%TurbOptions%MODVISH = Smagorinsky_
                STAT                              = SUCCESS_

            case default
                Me%TurbOptions%MODVISH = null_int
                STAT                              = NOT_FOUND_ERR_

        end select case1

    end subroutine CheckHorizontalTurbulenceOption 

    !--------------------------------------------------------------------------

    subroutine InicFileHorizontalModel

        !Local-----------------------------------------------------------------
        integer                         :: ClientNumber
        logical                         :: BlockFound
        type(T_PropertyID)              :: ID
        integer                         :: STAT_CALL

        !----------------------------------------------------------------------

        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                    viscosity_h_begin, viscosity_h_end, BlockFound, &
                                    STAT = STAT_CALL)

        if(blockfound) then

            call ConstructFillMatrix  (PropertyID           = ID,                               &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill3D       = Me%ExternalVar%WaterPoints3D,     &
                                       Matrix3D             = Me%Viscosity%HorizontalCenter,    &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            
            if (STAT_CALL  /= SUCCESS_) stop 'InicFileHorizontalModel - ModuleTurbulence - ERR12a'

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'InicFileHorizontalModel - ModuleTurbulence - ERR12a1'

            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InicFileHorizontalModel - ModuleTurbulence - ERR12a2'

            call KillFillMatrix(ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InicFileHorizontalModel - ModuleTurbulence- ERR12a3'

        endif

    end subroutine InicFileHorizontalModel

    !--------------------------------------------------------------------------

    subroutine InicFileVerticalModel

        !Local-----------------------------------------------------------------
        integer                         :: ClientNumber
        logical                         :: BlockFound
        type(T_PropertyID)              :: ID
        integer                         :: STAT_CALL

        !----------------------------------------------------------------------

        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                              &
                                    viscosity_v_begin, viscosity_v_end, BlockFound,             &
                                    STAT = STAT_CALL)
        if (STAT_CALL  /= SUCCESS_) stop 'InicFileVerticalModel - ModuleTurbulence - ERR01'


        if(BlockFound) then

            call ConstructFillMatrix  (PropertyID           = ID,                               &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill3D       = Me%ExternalVar%WaterPoints3D,     &
                                       Matrix3D             = Me%Viscosity%Vertical,            &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'InicFileVerticalModel - ModuleTurbulence - ERR02'

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'InicFileVerticalModel - ModuleTurbulence - ERR03'

            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'InicFileVerticalModel - ModuleTurbulence - ERR04'

            call KillFillMatrix(ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'InicFileVerticalModel - ModuleTurbulence - ERR05'

        endif


    end subroutine InicFileVerticalModel

    !--------------------------------------------------------------------------

    subroutine InicConstantHorizontalModel

        !Local-----------------------------------------------------------------
        integer                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                     :: i, j, k
        integer                     :: iflag, STAT_CALL
        real                        :: HorizontalCenterViscosity
                                                                          
        !----------------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        !Searches for the initial horizontal viscosity
        call GetData(HorizontalCenterViscosity,                     &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'VISCOSITY_H',                  &
                     ClientModule = 'ModuleTurbulence',             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                  &
            stop 'InicConstantHorizontalModel - ModuleTurbulence - ERR10'

        if (iflag == 0) then

            write(*,*) 
            write(*,*) 'Initial Horizontal Viscosity not defined'
            write(*,*) 'Keyword - VISCOSITY_H'
            stop       'InicConstantHorizontalModel - ModuleTurbulence - ERR12'

        else
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExternalVar%WaterPoints3D (i, j, k) == WaterPoint) then
                    Me%Viscosity%HorizontalCenter(i, j, k) = HorizontalCenterViscosity
                end if

            enddo
            enddo
            enddo
        endif


    end subroutine InicConstantHorizontalModel

    !--------------------------------------------------------------------------

    subroutine InicEstuarySmagorinskyModel

        !Local-----------------------------------------------------------------
        integer                     :: iflag, STAT_CALL
        
        !----------------------------------------------------------------------

        !Maximum horizontal viscosity
        call GetData(Me%TurbVar%MINHorizontalViscosity,         &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromFile,                   &
                     keyword      = 'VISH_REF',                 &
                     ClientModule = 'ModuleTurbulence',         &
                     default      = 50.0,                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                              &
            stop 'InicEstuarySmagorinskyModel - ModuleTurbulence - ERR01'

cd1 :   if (iflag .EQ. 0) then
            write(*,*) 
            write(*,*) 'Maximum Horizontal Viscosity not defined'
            write(*,*) 'Keyword - VISH_REF'
            write(*,*) 'Value set to ', Me%TurbVar%MINHorizontalViscosity
            write(*,*) 'InicEstuarySmagorinskyModel - ModuleTurbulence - WRN01'
            write(*,*) 
        end if cd1

        !----------------------------------------------------------------------

    end subroutine InicEstuarySmagorinskyModel


    !--------------------------------------------------------------------------


    subroutine InicEstuaryModel

        !Local-----------------------------------------------------------------
        integer                 :: iflag, STAT_CALL
                                  
        !----------------------------------------------------------------------

                     
        !Reference depth
        call GetData(Me%TurbVar%ReferenceDepth,             &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'HREF_VIS',             &
                     ClientModule = 'ModuleTurbulence',     &
                     default      = 10.0,                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                          &
            stop 'InicEstuaryModel - ModuleTurbulence - ERR01'

cd1 :   if (iflag .EQ. 0) then
            write(*,*) 
            write(*,*) 'Reference depth not defined'
            write(*,*) 'Keyword - HREF_VIS'
            write(*,*) 'Value set to ', Me%TurbVar%ReferenceDepth
            write(*,*) 'InicEstuaryModel - ModuleTurbulence - WRN01'
            write(*,*) 
        end if cd1


        !Reference velocity
        call GetData(Me%TurbVar%ReferenceVelocity,              &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromFile,                   &
                     keyword      = 'VREF_VIS',                 &
                     ClientModule = 'ModuleTurbulence',         &
                     default      = 1.0,                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                              &
            stop 'InicEstuaryModel - ModuleTurbulence - ERR02'

cd2 :   if (iflag .EQ. 0) then
            write(*,*) 
            write(*,*) 'Reference depth not defined'
            write(*,*) 'Keyword - VREF_VIS'
            write(*,*) 'Value set to ', Me%TurbVar%ReferenceVelocity
            write(*,*) 'InicEstuaryModel - ModuleTurbulence - WRN02'
            write(*,*) 
        end if cd2

    end subroutine InicEstuaryModel

    !--------------------------------------------------------------------------


    subroutine InicSmagorinskyModel
        
        
        !Local-----------------------------------------------------------------
        integer                     :: iflag, STAT_CALL
                                  
        !----------------------------------------------------------------------
                      
        !Reference depth
        call GetData(Me%TurbVar%HORCON,                     &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'HORCON',               &
                     ClientModule = 'ModuleTurbulence',     &
                     default      = 0.2,                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                          &
            stop 'InicEstuaryModel - ModuleTurbulence - ERR01'

cd1 :   if (iflag .EQ. 0) then
            write(*,*) 
            write(*,*) 'Keyword - HORCON'
            write(*,*) 'Value set to ', Me%TurbVar%HORCON
            write(*,*) 'InicEstuaryModel - ModuleTurbulence - WRN01'
            write(*,*) 
        end if cd1

cd2 :   if ((Me%TurbVar%HORCON .LT. 0.0) .OR.               &
            (Me%TurbVar%HORCON .GT. 1.0)) then
            write(*,*) 
            write(*,*) 'Keyword - HORCON'
            write(*,*) 'Value should be between 0.0 and 1.0'
            write(*,*) 'Actual value is ', Me%TurbVar%HORCON
        end if cd2

    end subroutine InicSmagorinskyModel
 
    !--------------------------------------------------------------------------

    subroutine InicConstantVerticalModel

        !External--------------------------------------------------------------
        integer             :: iflag, STAT_CALL
        real                :: VerticalViscosity
        real                :: VerticalPrandtlNumber
        real                :: VerticalMixLength

        !Local-----------------------------------------------------------------
        integer             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer             :: i, j, k
                                                      
        !----------------------------------------------------------------------

        call GetGeometryWaterColumn(Me%ObjGeometry, WaterColumn = Me%ExternalVar%HT, &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                          &
            stop 'InicConstantVerticalModel - ModuleTurbulence - ERR01'

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        !Searches for the initial vertical Prandtl number
        call GetData(VerticalPrandtlNumber,                 &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'PRANDTL_0',            &
                     ClientModule = 'ModuleTurbulence',     &
                     default      = 1.,                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                          &
            stop 'InicConstantVerticalModel - ModuleTurbulence - ERR02'


        !Searches for the initial vertical viscosity
        call GetData(VerticalViscosity,                     &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     ClientModule = 'ModuleTurbulence',     &
                     keyword      = 'VISCOSITY_V',          &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                          &
            stop 'InicConstantVerticalModel - ModuleTurbulence - ERR04'

cd1:    if (iflag == 0) then

            write(*,*) 
            write(*,*) 'Initial Vertical Viscosity not defined'
            write(*,*) 'Keyword - VISCOSITY_V'
            stop       'InicConstantVerticalModel - ModuleTurbulence - ERR05'

        else cd1

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExternalVar%WaterPoints3D(i, j, k)   == WaterPoint)   then

                    Me%Viscosity%Vertical(i, j, k)          = VerticalViscosity

                    Me%TurbVar%VertPrandtlNumber(i, j, k)   = VerticalPrandtlNumber

                    Me%Diffusivity%Vertical(i, j, k)        = VerticalViscosity * VerticalPrandtlNumber
                                                              
                endif

            enddo
            enddo
            enddo

        endif cd1

        !Searches for the constant vertical mixing length
        call GetData(VerticalMixLength,                       &
                     Me%ObjEnterData, iflag,                  &
                     SearchType   = FromFile,                 &
                     default      = 10.,                      &
                     ClientModule = 'ModuleTurbulence',       &
                     keyword      = 'MIXLENGTH_V',            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                            &
            stop 'InicConstantVerticalModel - ModuleTurbulence - ERR06'

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint)   then

                Me%TurbVar%MixingLengthZ(i, j, k) = VerticalMixLength 

            endif

        enddo
        enddo
        enddo

        call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%HT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InicConstantVerticalModel - ModuleTurbulence - ERR08'

    end subroutine InicConstantVerticalModel


    !--------------------------------------------------------------------------


    subroutine InicLayersVerticalModel

        !Local-----------------------------------------------------------------
        integer                         :: iflag, STAT_CALL
        real, dimension(:), allocatable :: VerticalViscosity
        real, dimension(:), allocatable :: VerticalPrandtlNumber
        real, dimension(:), allocatable :: VerticalMixLength
        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                         :: i, j, k

        !----------------------------------------------------------------------
        
        call GetGeometryWaterColumn(Me%ObjGeometry, WaterColumn = Me%ExternalVar%HT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'InicConstantVerticalModel - ModuleTurbulence - ERR01'

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        allocate(VerticalViscosity    (KLB:KUB))
        allocate(VerticalPrandtlNumber(KLB:KUB))
        allocate(VerticalMixLength    (KLB:KUB))


        !Searches for the initial vertical Prandtl number
        call GetData(VerticalPrandtlNumber,                 &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'PRANDTL_0',            &
                     ClientModule = 'ModuleTurbulence',     &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                          &
            stop 'InicConstantVerticalModel - ModuleTurbulence - ERR03'

        if (iflag /= KUB) then

            write(*,*) 
            write(*,*) 'Prandtl number by layers not defined'
            write(*,*) 'Keyword - PRANDTL_0'
            stop       'InicConstantVerticalModel - ModuleTurbulence - ERR04'

        endif


        !Searches for the initial vertical viscosity
        call GetData(VerticalViscosity,                     &
                     Me%ObjEnterData, iflag,               &
                     SearchType   = FromFile,               &
                     ClientModule = 'ModuleTurbulence',     &
                     keyword      = 'VISCOSITY_V',          &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                          &
            stop 'InicLayersVerticalModel - ModuleTurbulence - ERR05'


        if (iflag /= KUB) then

            write(*,*) 
            write(*,*) 'vertical viscosity by layers not defined'
            write(*,*) 'Keyword - VISCOSITY_V'
            stop       'InicConstantVerticalModel - ModuleTurbulence - ERR06'

        endif

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExternalVar%WaterPoints3D(i, j, k)   == WaterPoint)   then

                Me%Viscosity%Vertical(i, j, k)          = VerticalViscosity(k)

                Me%TurbVar%VertPrandtlNumber(i, j, k)   = VerticalPrandtlNumber(k)

                Me%Diffusivity%Vertical(i, j, k)        = VerticalViscosity(k)   &
                                                                   * VerticalPrandtlNumber(k)

            endif

        enddo
        enddo
        enddo


        !Searches for the constant vertical mixing length
        call GetData(VerticalMixLength,                 &
                     Me%ObjEnterData, iflag,           &
                     SearchType   = FromFile,           &
                     ClientModule = 'ModuleTurbulence', &
                     keyword      = 'MIXLENGTH_V',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                      &
            stop 'InicLayersVerticalModel - ModuleTurbulence - ERR07'

        if (iflag /= KUB) then

            write(*,*) 
            write(*,*) 'vertical mixing length by layers not defined'
            write(*,*) 'Keyword - MIXLENGTH_V'
            stop       'InicConstantVerticalModel - ModuleTurbulence - ERR08'

        endif

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExternalVar%WaterPoints3D(i, j, k)   == WaterPoint)   then

                Me%TurbVar%MixingLengthZ(i, j, k) = VerticalMixLength(k) 

            endif

        enddo
        enddo
        enddo

        deallocate(VerticalViscosity)
        deallocate(VerticalPrandtlNumber)
        deallocate(VerticalMixLength)


        call UnGetGeometry (Me%ObjGeometry,         &
                            Me%ExternalVar%HT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                  &
            stop 'InicLayersVerticalModel; ModuleTurbulence. ERR10.'


        !----------------------------------------------------------------------

    end subroutine InicLayersVerticalModel

    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------

    subroutine InicMixingLengthHorizontal

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------

        integer :: iflag, STAT_CALL

        real    :: MixingLengthH
        real    :: Nyquist

        !Local-----------------------------------------------------------------

        logical :: aux1, aux2
!        logical :: CalcML 

        integer :: ILB, IUB, JLB, JUB, KLB, KUB
        integer :: I,J,K
                                  
        !----------------------------------------------------------------------


        aux1 = .FALSE.
        aux2 = .FALSE.

        Me%MixingLengthH = OFF


        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB


        call GetData(MixingLengthH,                                             &
                     Me%ObjEnterData, iflag,                         &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'ModuleTurbulence',                         &
                     keyword      ='CONST_MIXING_LENGTH_HORIZONTAL',            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'InicMixingLengthHorizontal - ModuleTurbulence - ERR02'

cd1 :   if (iflag .EQ. 1) then 
            aux1 = .TRUE.

cd3 :       if (.NOT. Me%MixingLengthH) then
                Me%TurbVar%MixingLengthX = MixingLengthH 
                Me%TurbVar%MixingLengthY = MixingLengthH

                Me%MixingLengthH = ON
            end if cd3
        end if cd1


        !Valor teorico que multiplicado pela                                                                             
        ! passo da malha obtem-se o comprimento de onda                                                              
        ! maior filtrado pelo modelo teoricamente.                                                                        
        ! Na realidade Nyquest = 4 e 5                                                                                    
        ! ler : Coastal Estuarial                                                                                    
        ! and Haurbour Engineers' Reference Book,                                                                         
        ! Edited Abbott, M. B. and Price, W. A., 1994                                                                
        ! " Diffusion, Dispersion and Sub-Grid parameterization"                                                          
        ! BEDFORD,K. W. (o Prof. Ramiro tem este livro)    
        call GetData(Nyquist,                                                 &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromFile,                                 &
                     keyword      ='NYQUIST',                                 &
                     Default      = 2.0,                                      &
                     ClientModule = 'ModuleTurbulence',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                            &
            stop 'InicMixingLengthHorizontal - ModuleTurbulence - ERR03'

        if (iflag .EQ. 1)                                                     &
            aux2 = .TRUE.

cd2 :   if (.NOT. Me%MixingLengthH) then
            call GetHorizontalGrid(Me%ObjHorizontalGrid,           &
                                   DUX = Me%ExternalVar%DUX,       &
                                   STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                        &
                stop 'InicMixingLengthHorizontal - ModuleTurbulence - ERR04'


            call GetHorizontalGrid(Me%ObjHorizontalGrid,           &
                                   DVY = Me%ExternalVar%DVY,       &
                                   STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                        &
                stop 'InicMixingLengthHorizontal - ModuleTurbulence - ERR05'



do1 :       do K = KLB, KUB
do2 :       do J = JLB, JUB
do3 :       do I = ILB, IUB
                Me%TurbVar%MixingLengthX(I,J,K) = Nyquist * Me%ExternalVar%DUX(I,J)
                Me%TurbVar%MixingLengthY(I,J,K) = Nyquist * Me%ExternalVar%DVY(I,J)
            end do do3
            end do do2
            end do do1

            Me%MixingLengthH = ON




            call UngetHorizontalGrid(Me%ObjHorizontalGrid,         &
                                     Me%ExternalVar%DUX,           & 
                                     STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_)                                      &
                stop 'InicMixingLengthHorizontal - ModuleTurbulence - ERR06'


            call UngetHorizontalGrid(Me%ObjHorizontalGrid,         &
                                     Me%ExternalVar%DVY,           &
                                     STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_)                                      &
                stop 'InicMixingLengthHorizontal - ModuleTurbulence - ERR07'
        end if cd2


        !keywords 
        ! CONST_MIXING_LENGTH_HORIZONTAL & NYQUIST
        ! found in the same file
        if (aux1 .AND. aux2)                                                  &
            stop 'InicMixingLengthHorizontal - ModuleTurbulence - ERR01'

        !----------------------------------------------------------------------

    end subroutine InicMixingLengthHorizontal



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine GetHorizontalViscosity(TurbulenceID, HorizontalCenterViscosity, &
                                                     HorizontalCornerViscosity, STAT)

        !Arguments-------------------------------------------------------------

        integer                                     :: TurbulenceID

        real, dimension(:,:,:), optional, pointer   :: HorizontalCenterViscosity
        real, dimension(:,:,:), optional, pointer   :: HorizontalCornerViscosity

        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------

        integer :: ready_              

        !Local-----------------------------------------------------------------

        integer :: STAT_             

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
        if ((ready_ == IDLE_ERR_     ) .OR.                                   &
            (ready_ == READ_LOCK_ERR_)) then

cd1 :       if (present(HorizontalCenterViscosity)) then
                call Read_Lock(mTURBULENCE_, Me%InstanceID)

                HorizontalCenterViscosity => Me%Viscosity%HorizontalCenter
            end if cd1
            

cd2 :       if (present(HorizontalCornerViscosity)) then
                call Read_Lock(mTURBULENCE_, Me%InstanceID)

                HorizontalCornerViscosity => Me%Viscosity%HorizontalCorner
            end if cd2
            
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHorizontalViscosity

    !--------------------------------------------------------------------------
            
     subroutine GetTurbulenceModelsList(Constant, Estuary,    Smagorinsky, Leendertsee,     &
                                        Backhaus, Pacanowski, Nihoul, TurbulenceEquation )
!                                        ,      TherryLacarrere, &
!                                        MellorYamada)
       
        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: constant      
        integer, optional, intent(OUT) :: Estuary       
        integer, optional, intent(OUT) :: smagorinsky       
        integer, optional, intent(OUT) :: leendertsee        
        integer, optional, intent(OUT) :: backhaus    
        integer, optional, intent(OUT) :: pacanowski    
        integer, optional, intent(OUT) :: nihoul      
        integer, optional, intent(OUT) :: TurbulenceEquation 
                 
        !----------------------------------------------------------------------
     
        if (present(constant       )) constant        = Constant_
        if (present(Estuary        )) Estuary         = Estuary_
        if (present(smagorinsky    )) smagorinsky     = Smagorinsky_
        if (present(leendertsee    )) leendertsee     = leendertsee_
        if (present(backhaus       )) backhaus        = backhaus_
        if (present(pacanowski     )) pacanowski      = pacanowski_
        if (present(nihoul         )) nihoul          = nihoul_
        if (present(TurbulenceEquation)) TurbulenceEquation = TurbulenceEquation_

        !----------------------------------------------------------------------

    end subroutine GetTurbulenceModelsList

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine GetVerticalViscosity(TurbulenceID, VerticalViscosityCenter, STAT)

        !Arguments-------------------------------------------------------------
              
        integer, optional, intent(OUT)  :: STAT
       
        real, dimension(:,:,:), pointer :: VerticalViscosityCenter

        integer                         :: TurbulenceID

        !External--------------------------------------------------------------

        integer :: ready_              
   
        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                   &
            (ready_ == READ_LOCK_ERR_)) then

                call Read_Lock(mTURBULENCE_, Me%InstanceID)

                VerticalViscosityCenter => Me%Viscosity%Vertical

            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetVerticalViscosity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetVerticalDiffusivity(TurbulenceID, VerticalDiffusivityCenter, STAT)

        !Arguments-------------------------------------------------------------
              
        integer, optional, intent(OUT) :: STAT
       
        real, dimension(:,:,:), pointer :: VerticalDiffusivityCenter
        integer                         :: TurbulenceID

        !External--------------------------------------------------------------

        integer :: ready_              

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                   &
            (ready_ == READ_LOCK_ERR_)) then

                call Read_Lock(mTURBULENCE_, Me%InstanceID)

                VerticalDiffusivityCenter => Me%Diffusivity%Vertical                        

            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetVerticalDiffusivity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetMixingLengthHorizontal(TurbulenceID, MixingLengthX, MixingLengthY, STAT)

        !Arguments-------------------------------------------------------------

        integer                                     :: TurbulenceID

        real, optional, pointer, dimension(:,:,:)   :: MixingLengthX, MixingLengthY

        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------

        integer :: ready_              
        
        !Local-----------------------------------------------------------------

        integer :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                   &
            (ready_ == READ_LOCK_ERR_)) then

cd4 :       if (Me%MixingLengthH) then

cd2 :           if (present(MixingLengthX)) then
                    call Read_Lock(mTURBULENCE_, Me%InstanceID)
                    MixingLengthX => Me%TurbVar%MixingLengthX
                end if cd2
            

cd3 :           if (present(MixingLengthY)) then
                    call Read_Lock(mTURBULENCE_, Me%InstanceID)
                    MixingLengthY => Me%TurbVar%MixingLengthY
                end if cd3
            
          
                STAT_ = SUCCESS_

            else cd4

                STAT_ = OFF_ERR_

            end if cd4
        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetMixingLengthHorizontal

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetMixingLengthVertical(TurbulenceID, Lupward, Ldownward, STAT)

        !Arguments-------------------------------------------------------------

        integer                                     :: TurbulenceID

        real, optional, pointer, dimension(:,:,:)   :: Lupward   
        real, optional, pointer, dimension(:,:,:)   :: Ldownward   

        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------

        integer :: ready_              
        integer :: STAT_CALL

        real, pointer, dimension(:,:,:) :: Aux

        !Local-----------------------------------------------------------------

        integer :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                   &
            (ready_ == READ_LOCK_ERR_)) then

            

cd3 :       if (present(Lupward     )) then
cd5 :           if (Me%TurbOptions%MODTURB == TurbulenceEquation_) then
                    call Read_Lock(mTURBULENCE_, Me%InstanceID)

                    call GetMixingLenghTurbGOTM(Me%ObjTurbGOTM, Lupward = Aux, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                              &
                        stop 'Subroutine GetMixingLengthVertical; module ModuleTurbulence. ERR01.'

                    Me%TurbVar%MixingLengthZ = Aux

                    Lupward => Me%TurbVar%MixingLengthZ

                    call UngetTurbGOTM(Me%ObjTurbGOTM, Aux, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                              &
                        stop 'Subroutine GetMixingLengthVertical; module ModuleTurbulence. ERR02.'
     
                else cd5
                    call Read_Lock(mTURBULENCE_, Me%InstanceID)

                    Lupward => Me%TurbVar%MixingLengthZ

                end if cd5
            end if cd3
            

cd4 :       if (present(Ldownward   )) then
cd6 :           if (Me%TurbOptions%MODTURB == TurbulenceEquation_) then
                    call Read_Lock(mTURBULENCE_, Me%InstanceID)

                    call GetMixingLenghTurbGOTM(Me%ObjTurbGOTM, Ldownward = Aux, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                              &
                        stop 'Subroutine GetMixingLengthVertical; module ModuleTurbulence. ERR03.'

                    Me%TurbVar%Ldownward = Aux
   
                    Ldownward => Me%TurbVar%Ldownward

                    call UngetTurbGOTM(Me%ObjTurbGOTM, Aux, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                              &
                        stop 'Subroutine GetMixingLengthVertical; module ModuleTurbulence. ERR04.'
     
                else cd6
                    call Read_Lock(mTURBULENCE_, Me%InstanceID)

                    Ldownward => Me%TurbVar%MixingLengthZ

                end if cd6
            end if cd4
            
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetMixingLengthVertical

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetMLD_Surf(TurbulenceID, MLD_Surf, STAT)

        !Arguments-------------------------------------------------------------

        integer                         :: TurbulenceID

        real, pointer, dimension(:,:)  :: MLD_Surf

        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------

        integer :: ready_              
        
        !Local-----------------------------------------------------------------

        integer :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

cd2 :       if (Me%TurbOptions%MLD_Calc) then
              
                call Read_Lock(mTURBULENCE_, Me%InstanceID)
                
                MLD_Surf => Me%TurbVar%MLD_Surf
       
                STAT_ = SUCCESS_

            else   cd2

                STAT_ = OFF_ERR_

            end if cd2

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine GetMLD_Surf

    !--------------------------------------------------------------------------

    subroutine GetContinuousGOTM(TurbulenceID, Continuous, ModelGOTM, STAT)

        !Arguments-------------------------------------------------------------

        integer                         ::TurbulenceID

        logical,           intent(OUT)  :: Continuous, ModelGOTM

        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------

        integer :: ready_              
        
        !Local-----------------------------------------------------------------

        integer :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

            call Read_Lock(mTURBULENCE_, Me%InstanceID)
            
            if (Me%TurbOptions%MODTURB == TurbulenceEquation_) then

                ModelGOTM = .true.

            else

                ModelGOTM = .false.

            endif

            Continuous = Me%TurbOptions%Continuous_Compute

            call Read_UnLock(mTurbulence_, Me%InstanceID, "GetContinuousGOTM")
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetContinuousGOTM

    !--------------------------------------------------------------------------
    
    subroutine GetTurbulenceOptions(TurbulenceID, NeedsTempSalinity, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: TurbulenceID
        logical,           intent(OUT)  :: NeedsTempSalinity
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                             &
            (ready_ == READ_LOCK_ERR_)) then

            if (Me%TurbOptions%MODTURB       == Constant_       .or.    &
                Me%TurbOptions%DensityMethod == ConstantDensity_) then

                NeedsTempSalinity = .false.

            else

                NeedsTempSalinity = .true.

            endif

            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetTurbulenceOptions

    !--------------------------------------------------------------------------

    subroutine UngetTurbulence3D(TurbulenceID, Array, STAT)

        !Arguments-------------------------------------------------------------

        real, pointer, dimension(:,:,:)     :: Array

        integer                             :: TurbulenceID  

        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_             

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)     

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mTurbulence_, Me%InstanceID, "UngetTurbulence3D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetTurbulence3D

    !--------------------------------------------------------------------------

    subroutine UngetTurbulence2D(TurbulenceID, Array, STAT)

        !Arguments-------------------------------------------------------------

        real, pointer, dimension(:,:)       :: Array

        integer                             :: TurbulenceID  

        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_             

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)     

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mTurbulence_, Me%InstanceID, "UngetTurbulence2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetTurbulence2D

    !----------------------------------------------------------------------

    subroutine SetTurbulenceBottomRugosity(TurbulenceID, BottomRugosity, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: TurbulenceID
        real, intent(IN)                :: BottomRugosity
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_              
        
        !Local-----------------------------------------------------------------
        integer                         :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then


            Me%ExternalVar%BottomRugosity = BottomRugosity            
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetTurbulenceBottomRugosity
    
    !----------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine Turbulence(TurbulenceID, VelocityX, VelocityY, VelocityZ, Chezy, &
                          Salinity, Temperature, STAT)

        !Arguments-------------------------------------------------------------
        integer, optional, intent(OUT)                  :: STAT
        integer                                         :: TurbulenceID
        real, dimension(:,:,:), pointer                 :: VelocityX
        real, dimension(:,:,:), pointer                 :: VelocityY
        real, dimension(:,:,:), pointer                 :: VelocityZ
        real, dimension(:,:,:), pointer                 :: Salinity
        real, dimension(:,:,:), pointer                 :: Temperature
        real, dimension(:,:  ), pointer                 :: Chezy

        !External--------------------------------------------------------------
        integer                                         :: ready_              
        integer                                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_         
        integer                                         :: ILB, IUB
        integer                                         :: JLB, JUB
        integer                                         :: KLB, KUB
        integer                                         :: I, J, K
        integer                                         :: CHUNK

        !----------------------------------------------------------------------                         

        if (MonitorPerformance) call StartWatch ("ModuleTurbulence", "Turbulence")

        STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_) 

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB
            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB

            !Time Properties - Actualises CurrentTime
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'Turbulence - ModuleTurbulence - ERR00'


            Me%ExternalVar%VelocityX          => VelocityX
            Me%ExternalVar%VelocityY          => VelocityY
            Me%ExternalVar%VelocityZ          => VelocityZ
            Me%ExternalVar%Chezy              => Chezy
            
i1 :        if      (.not.(Me%TurbOptions%MODTURB .EQ. Constant_   .or.       &
                           Me%TurbOptions%MODTURB .EQ. file2D_)) then            
                Me%ExternalVar%Salinity           => Salinity
                Me%ExternalVar%Temperature        => Temperature
                
            endif i1


            call GetGridData(Me%ObjGridData, Me%ExternalVar%Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR01'

            call GetGeometryDistances (Me%ObjGeometry,                                  &
                                       SZZ  = Me%ExternalVar%SZZ,                       &
                                       ZCellCenter = Me%ExternalVar%ZCellCenter,        &
                                       DZZ  = Me%ExternalVar%DZZ,                       &
                                       DWZ  = Me%ExternalVar%DWZ,                       &
                                       STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR02'

            call GetGeometryWaterColumn(Me%ObjGeometry,                                 &
                                        WaterColumn = Me%ExternalVar%HT,                &
                                        STAT        = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR03'

 
            call GetGeometryKFloor(Me%ObjGeometry, Z = Me%ExternalVar%KFloorZ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR04'



            call GetWaterPoints3D (Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR05'


            call GetOpenPoints3D (Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR06'


            call GetHorizontalGrid     (Me%ObjHorizontalGrid,                           &
                                        DVY  = Me%ExternalVar%DVY,                      &
                                        DUX  = Me%ExternalVar%DUX,                      &
                                        DYY  = Me%ExternalVar%DYY,                      &
                                        DZX  = Me%ExternalVar%DZX,                      &
                                        DZY  = Me%ExternalVar%DZY,                      &
                                        STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR07'


            !Gets pointers to ComputeFacesU3D, ComputeFacesV3D, ComputeFacesW3D
            call GetComputeFaces3D (Me%ObjMap,                                          &
                                    ComputeFacesU3D = Me%ExternalVar%ComputeFacesU3D,   &
                                    ComputeFacesV3D = Me%ExternalVar%ComputeFacesV3D,   &
                                    ComputeFacesW3D = Me%ExternalVar%ComputeFacesW3D,   &
                                    STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR08'

            call GetImposedTangentialFaces(Me%ObjMap,                                   &
                                           Me%ExternalVar%ImposedTangentialFacesU,      &
                                           Me%ExternalVar%ImposedTangentialFacesV,      &
                                           STAT = STAT_CALL)      
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR09'


            call GetImposedNormalFaces(Me%ObjMap,                                       &
                                       Me%ExternalVar%ImposedNormalFacesU,              &  
                                       Me%ExternalVar%ImposedNormalFacesV,              &
                                       STAT = STAT_CALL)      
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR10'

        
            !Vertical model

            !Computes Brunt-Vaisalla and Prandtl frequency and Richardson number
            if (Me%TurbOptions%MODTURB /= Constant_) call Richardson
               
              
            !Vertical model options
cd3 :       if      (Me%TurbOptions%MODTURB .EQ. Constant_   .or.       &
                     Me%TurbOptions%MODTURB .EQ. file2D_) then

                continue    !If turbulence is Constant_ it does not change in time and space
        
            elseif (Me%TurbOptions%MODTURB .EQ. leendertsee_       ) then cd3
                    
                call LeendertseeModel
        
            elseif (Me%TurbOptions%MODTURB .EQ. backhaus_          ) then cd3
                 
                call BackhausModel
            
            elseif (Me%TurbOptions%MODTURB .EQ. pacanowski_        ) then cd3
                
                call PacanowskiModel
        
            elseif (Me%TurbOptions%MODTURB .EQ. nihoul_            ) then cd3

                call NihoulModel

            elseif (Me%TurbOptions%MODTURB .EQ. TurbulenceEquation_) then cd3

                call TurbulenceEquationModel

            else cd3

                stop 'Turbulence - ModuleTurbulence - ERR11'

            end if cd3


            !Here the dependence on stratification is quantified

            !Manuel, 1999
            !A constant value is added to get viscosity. By default it is
            !the molecular viscosity , but another value can be given 
            !(keyword: BAckground_viscosity) to account for any kind of
            !unresolved mixing. A constant diffusivity could be added in waterproperties to
            !every property

cd2 :       if (Me%TurbOptions%MODTURB .ne. Constant_ .and. &
                Me%TurbOptions%MODTURB .ne. File2D_) then

                CHUNK = CHUNK_J(Me%Size%JLB, Me%Size%JUB)
                !$OMP PARALLEL SHARED(CHUNK) PRIVATE(I,J)

do1 :           do K = KLB, KUB
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :           do J = JLB, JUB
do3 :           do I = ILB, IUB

                    if(Me%ExternalVar%WaterPoints3D(i,j,k) == WaterPoint)then
                    
                        Me%Viscosity%Vertical(i,j,k) = Me%Viscosity%Vertical(i,j,k) + &
                                                       Me%Viscosity%Background
                    end if
                                                            
                end do do3
                end do do2
                !$OMP END DO NOWAIT
                end do do1

                !$OMP END PARALLEL

            end if cd2

                !Horizontal turbulence model
cd4 :       if     (Me%TurbOptions%MODVISH .EQ. Constant_   ) then 
                
                continue    !If turbulence is Constant_ it does not change in time and space
                
            elseif (Me%TurbOptions%MODVISH .EQ. File2D_     ) then cd4
            
                continue
            
            elseif (Me%TurbOptions%MODVISH .EQ. Estuary_    ) then cd4
            
                call EstuaryModel
            
            elseif (Me%TurbOptions%MODVISH .EQ. Smagorinsky_) then cd4
            
                call SmagorinskyModel
            
            else cd4
            
                stop 'Turbulence - ModuleTurbulence - ERR12'
            
            end if cd4

            call TurbulentViscosity_CellCorner

            call OutPut_Results_HDF

            call UngetGridData (Me%ObjGridData, Me%ExternalVar%Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR14'

            call UngetGeometry (Me%ObjGeometry, Me%ExternalVar%SZZ,     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR15'

            call UngetGeometry (Me%ObjGeometry, Me%ExternalVar%DZZ,     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR16'


            call UngetGeometry (Me%ObjGeometry, Me%ExternalVar%DWZ,     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR17'

            call UngetGeometry (Me%ObjGeometry, Me%ExternalVar%HT,      STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR18'

            call UngetGeometry (Me%ObjGeometry, Me%ExternalVar%KFloorZ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR19'

            call UngetGeometry (Me%ObjGeometry, Me%ExternalVar%ZCellCenter, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR20'


            call UngetHorizontalGrid(Me%ObjHorizontalGrid,         &
                                     Me%ExternalVar%DUX,           &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR21'


            call UngetHorizontalGrid(Me%ObjHorizontalGrid,         &
                                     Me%ExternalVar%DVY,           &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR22'

            call UngetHorizontalGrid(Me%ObjHorizontalGrid,         &
                                     Me%ExternalVar%DYY,           &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR23'


            call UngetHorizontalGrid(Me%ObjHorizontalGrid,         &
                                     Me%ExternalVar%DZX,           &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR24'

            
            call UngetHorizontalGrid(Me%ObjHorizontalGrid,         &
                                     Me%ExternalVar%DZY,           &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR25'

            call UngetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesU3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR26'

            call UngetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesV3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR27'

            call UngetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesW3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR28'

            call UngetMap(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR29'
                       
            call UngetMap(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR30'

            call UngetMap(Me%ObjMap, Me%ExternalVar%ImposedTangentialFacesU, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR31'
            
            call UngetMap(Me%ObjMap, Me%ExternalVar%ImposedTangentialFacesV, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR32'

            call UngetMap(Me%ObjMap, Me%ExternalVar%ImposedNormalFacesU, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR33'

            call UngetMap(Me%ObjMap, Me%ExternalVar%ImposedNormalFacesV, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR34'

            nullify(Me%ExternalVar%VelocityX)
            nullify(Me%ExternalVar%VelocityY)
            nullify(Me%ExternalVar%Chezy    )

            STAT_ = SUCCESS_
        else               

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleTurbulence", "Turbulence")

        !----------------------------------------------------------------------

    end subroutine Turbulence

    !--------------------------------------------------------------------------

    subroutine OutPut_Results_HDF

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: NextOutput
        type(T_Time)                                    :: NextProfileOutput
        real                                            :: Year, Month, Day
        real                                            :: Hour, Minute, Second

        !Begin-----------------------------------------------------------------


        call GetOpenPoints3D (Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR06'

            call GetGeometryDistances (Me%ObjGeometry,                                  &
                                       SZZ  = Me%ExternalVar%SZZ,                       &
                                       STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR02'

            call GetGeometryWaterColumn(Me%ObjGeometry,                                 &
                                        WaterColumn = Me%ExternalVar%HT,                &
                                        STAT        = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR03'

 
            call GetGeometryKFloor(Me%ObjGeometry, Z = Me%ExternalVar%KFloorZ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR04'



        if((Me%TurbOptions%MODTURB .eq. TurbulenceEquation_) .and.    &
           ((Me%Output%TimeSerie).or.(Me%OutPut%ON).or.(Me%OutPut%ProfileON)))then
                     
            call  GetTurbGOTM_TurbEq(Me%ObjTurbGOTM,                  &
                                     Me%TurbVar%TKE,                  &
                                     Me%TurbVar%eps,                  &
                                     Me%TurbVar%L,                    &
                                     Me%TurbVar%P,                    &
                                     Me%TurbVar%B,                    &
                                     STAT  = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'OutPut_Results_HDF - ModuleTurbulence - ERR10'

        end if

        if(Me%TurbOptions%MLD_Calc)  call ComputeMixedLayerDepth

        if(Me%Output%TimeSerie    )  call OutPut_TimeSeries

        if(Me%Output%ProfileON) then

            call GetProfileNextOutputTime(Me%ObjProfile, NextProfileOutput, STAT = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)stop 'OutPut_Results_HDF - ModuleTurbulence - ERR20'


            if (Me%ExternalVar%Now .ge. NextProfileOutput)then
                call Output_Profile
            endif

        end if

            !Output
        if (Me%OutPut%ON) then
    
            NextOutPut = Me%OutPut%NextOutPut

            if (Me%ExternalVar%Now >= Me%OutPut%OutTime(NextOutPut)) then

                if (Me%ExternalVar%EndTime == Me%OutPut%OutTime(NextOutPut)) then

                    Me%OutPut%Run_End = .true.

                else

                    Me%OutPut%Run_End = .false.

                endif

                !Output in HDF
                call Write_HDF5_Format

                Me%OutPut%NextOutPut = NextOutPut + 1

            end if

        end if

        !Ungets the variables from TurbGOTM needed for output
        if((Me%TurbOptions%MODTURB == TurbulenceEquation_).and.     &
           ((Me%Output%TimeSerie).or.(Me%OutPut%ON).or.(Me%OutPut%ProfileON)))  then
                
                call UnGetTurbGOTM_TurbEq (Me%ObjTurbGOTM,          &
                                           Me%TurbVar%TKE,          &
                                           Me%TurbVar%eps,          &
                                           Me%TurbVar%L,            &
                                           Me%TurbVar%P,            &
                                           Me%TurbVar%B,            &
                                           STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'OutPut_Results_HDF - ModuleTurbulence - ERR30'

        end if 

        !Write restart files
        if(Me%OutPut%WriteRestartFile .and. .not. Me%OutPut%Run_End)then

            if(Me%ExternalVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then

                call Write_Final_Turbulence_File(Me%ObjTurbGOTM,                            &
                                                 Overwrite = Me%Output%RestartOverwrite,    &
                                                 STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'OutPut_Results_HDF - ModuleTurbulence - ERR40'

                Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1

                call ExtractDate(Me%ExternalVar%Now, Year = Year, Month  = Month,  Day    = Day, &
                                                    Hour = Hour, Minute = Minute, Second = Second)

                call SetError(WARNING_, INTERNAL_, "Turbulence restart file saved          : ", &
                              Year, Month, Day, Hour, Minute, Second)

            end if

        endif

        call UngetMap(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR30'

            call UnGetGeometry(Me%ObjGeometry,                                  &
                              Me%ExternalVar%SZZ,                       &
                              STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR02'

            call UnGetGeometry(Me%ObjGeometry,                                 &
                              Me%ExternalVar%HT,                &
                                        STAT        = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR03'

 
            call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%KFloorZ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Turbulence - ModuleTurbulence - ERR04'


    
    end subroutine OutPut_Results_HDF
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! Vertical viscosities computed according to Leendertsee and Liu 
    ! This parameterization was designed for coastal zones
    ! Depends on depth, velocity gradient and stratification 
    ! (the strat. dependence is taken from Nihoul). Based on Prandtl mixing length...
    ! Last modified:
    ! 2001: Manuel RV 
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine LeendertseeModel

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        real    :: Z_H 
        real    :: CMIST
        real    :: VISC_V
        real    :: Aux

        integer :: ILB, IUB
        integer :: JLB, JUB
        integer :: KLB, KUB

        integer :: I, J, K
        integer                                         :: CHUNK

        !----------------------------------------------------------------------


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB


!        call Richardson

        if (MonitorPerformance) call StartWatch ("ModuleTurbulence", "LeendertseeModel")
        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,aux,Z_H,CMIST,VISC_V)

do1 :   do K = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do J = JLB, JUB
do3 :   do I = ILB, IUB
            
cd1 :       if (Me%ExternalVar%ComputeFacesW3D(i, j, k) == Compute) then


!Manuel !Vertical Prandtl number was BEFORE calculated in Richardson PRandtl = gama* e**(-gam2*Rich)
!       If changed it should be moved to a subroutine where VErtPRandtl number is calculated from different parametrizations

                Aux = max (0., Me%TurbVar%Richardson(I,J,K)) !Nihoul formula only makes sense for positive Rich

                Me%TurbVar%VertPrandtlNumber(I,J,K) = GAMA * Exp(-GAMA2 *  Aux) ! Nihoul, 1984 

                Z_H =Me%ExternalVar%Bathymetry(I,J) -     &
                      Me%ExternalVar%SZZ(I,J,K)                              
!                     &- Me%ExternalVar%DWZ(I,J,K)/2.)        & !Z_H nas faces , nao e?
!                    / Me%ExternalVar%HT(I,J)

                CMIST = CVK * (Z_H + Me%ExternalVar%BottomRugosity) &
                              * SQRT(1.0 - Z_H/ Me%ExternalVar%HT(I,J))  
                                                
                CMIST = CMIST * EXP(CNH * Aux) 

                Me%TurbVar%MixingLengthZ(I,J,K) = MIN(CMIST, Me%TurbVar%MAXMixingLength)

                VISC_V =(Me%TurbVar%MixingLengthZ(I,J,K) * Me%TurbVar%MixingLengthZ(I,J,K)) &
                       * SQRT(Me%TurbVar%FPRANDTL(I,J,K))
                
                Me%Viscosity%Vertical(I,J,K) = VISC_V 

                Me%Diffusivity%Vertical(i, j, k)      = Me%Viscosity%Vertical(i,j,k)   &
                                                                       *Me%TurbVar%VertPrandtlNumber(i,j,k)
                
! Manolo. This is done from file 
! Me%Viscosity%Vertical(I,J,K) = VISC_V + MolecularViscosity     !correccao pela visc. molec

            end if cd1
        end do do3
        end do do2
        !$OMP END DO
        end do do1

        !$OMP END PARALLEL
        if (MonitorPerformance) call StopWatch ("ModuleTurbulence", "LeendertseeModel")

        !----------------------------------------------------------------------

    end subroutine LeendertseeModel

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! Vertical viscosities computed according to Backhaus and Hainbucher (1987)
    ! This parameterization was designed for the North Sea.
    ! Depends on depth, velocity and stratification.
    ! Last modified:
    ! 2001: Manuel RV 
    !--------------------------------------------------------------------------

    subroutine BackhausModel

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        real    :: Rich
 
        integer :: ILB, IUB
        integer :: JLB, JUB
        integer :: KLB, KUB

        integer :: I, J, K
        integer                                         :: CHUNK

        !----------------------------------------------------------------------


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB


!        call Richardson


        !ATENÇÃO: Se RICH<0. pode-se obter uma potência impossivel de calcular porque
        !         o expoente CBK2 é real.Por este motivo para a correlação de Backhaus
        !         (e apenas neste caso) impôe-se um valor minimo de RICH o que impede
        !         com esta correlação o estudo de zonas com um gradiente de densidade
        !         instável(salt fingers,zonas onde há formação de gelo,etc...).
        !   We have imposed rich = max (0., Me%TurbVar%Richardson(I,J,K)) to allow computations,
        !   but what's written must be taken into account!

        if (MonitorPerformance) call StartWatch ("ModuleTurbulence", "BackhausModel")
        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,rich)
        
do1 :   do K = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do J = JLB, JUB
do3 :   do I = ILB, IUB
            
cd1 :       if (Me%ExternalVar%ComputeFacesW3D(i, j, k) == Compute) then

!cd2 :           if (Me%TurbVar%Richardson(I,J,K) .LT. (-0.99/CBK1)) then
!                    write(*,*           ) 
!                    write(*,*           ) 'Backhaus vertical viscosity model can not be applied in areas with'
!                    write(*,*           ) 'instable vertical gradient.'
!                    write(*,*           ) 'Richardson number is ', Me%TurbVar%Richardson(I,J,K)
!                    write(*,*           ) 'AT point (I,J,K) ', I, J, K
!                    stop                  'Subroutine BackhausModel; module ModuleTurbulence. ERR01.'
!                end if cd2
                 
                 rich = max (0., Me%TurbVar%Richardson(I,J,K))

               
!                Me%TurbVar%VertPrandtlNumber(I,J,K) = GAMA * Exp(-GAMA2 *     &
!                                                         Me%TurbVar%Richardson(I,J,K) ) ! Nihoul, 1984 
                
                Me%Viscosity%Vertical(I,J,K) =  &
!                CBK3  + ! CBK3 is only a background viscosity                                       
                          CBK4 * Me%ExternalVar%HT(I,J) * ABS(Me%TurbVar%VMOD(I,J,K)) &
                          * (1.0 + CBK1 * rich)**CBK2

               Me%Diffusivity%Vertical(i, j, k)      = Me%Viscosity%Vertical(i,j,k)   &
                          * (1.0 + rich)**CBK5


            end if cd1
        end do do3
        end do do2
        !$OMP END DO
        end do do1

        !$OMP END PARALLEL
        if (MonitorPerformance) call StopWatch ("ModuleTurbulence", "BackhausModel")

        !----------------------------------------------------------------------

    end subroutine BackhausModel

    !----------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! Vertical viscosities computed according to Pacanowski and Philander (1981)
    ! This parameterization was designed for the tropical ocean! Only depends on startification
    ! Last modified:
    ! 2001: Manuel RV 
    !--------------------------------------------------------------------------

    subroutine PacanowskiModel

        !Arguments-------------------------------------------------------------
        real    :: RICH
        integer :: ILB, IUB
        integer :: JLB, JUB
        integer :: KLB, KUB

        integer :: I, J, K
        integer                                         :: CHUNK
        !----------------------------------------------------------------------


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        if (MonitorPerformance) call StartWatch ("ModuleTurbulence", "PacanowskiModel")
        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,RICH)
        
do1 :   do K = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do J = JLB, JUB
do3 :   do I = ILB, IUB
cd1 :       if (Me%ExternalVar%ComputeFacesW3D(i, j, k) == Compute) then

               RICH = max(0., Me%TurbVar%Richardson(I,J,K))
              
               Me%Viscosity%Vertical(I,J,K)           =             &
!                CPW3 +                &   !CPW3 is the molecular viscosity
                    CPW4 * (1.0 + CPW1 * RICH)**ICPW2

               Me%Diffusivity%Vertical(i, j, k)      =    &
                   Me%Viscosity%Vertical(I,J,K)  / (1.0 + CPW1 * RICH)           

            end if cd1
        end do do3
        end do do2
        !$OMP END DO
        end do do1

        !$OMP END PARALLEL
        if (MonitorPerformance) call StopWatch ("ModuleTurbulence", "PacanowskiModel")

        !----------------------------------------------------------------------

    end subroutine PacanowskiModel

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! Vertical viscosities computed according to Leendertsee and Liu 
    ! This parameterization was designed for coastal zones
    ! Depends on depth, vel. gradient and stratification. Based on Prandtl mixing length... 
    ! Last modified:
    ! 2001: Manuel RV 
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine NihoulModel

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        real    :: Z_H 
        real    :: CMIST
        real    :: VISC_V
        real    :: Aux    
        integer :: ILB, IUB
        integer :: JLB, JUB
        integer :: KLB, KUB

        integer :: I, J, K
        integer :: CHUNK

        !----------------------------------------------------------------------


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
     
        if (MonitorPerformance) &
            call StartWatch ("ModuleTurbulence", "NihoulModel")

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,Aux,Z_H,CMIST,VISC_V)

do1 :   do K = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do J = JLB, JUB
do3 :   do I = ILB, IUB
            
cd1 :       if (Me%ExternalVar%ComputeFacesW3D(i, j, k) == Compute)        then
                !Vertical Prandtl number was BEFORE calculated in subroutine Richardson PRandtl = gama* e**(-gam2*Rich)
                ! If changed it should be moved to a subroutine where VErtPRandtl number is calculated from different 
                ! parametrizations
                Aux = max(0., Me%TurbVar%Richardson(I,J,K))
                Me%TurbVar%VertPrandtlNumber(I,J,K) = GAMA * Exp(-GAMA2 *  Aux) ! Nihoul, 1984 


                !Z_H computed at centers
                !Z_H distance to bottom
                Z_H =(Me%ExternalVar%Bathymetry(I,J) -    &
                      Me%ExternalVar%SZZ(I,J,K))
!                      Me%ExternalVar%DWZ(I,J,K)/2.)       &   !Manolo. Comentar? Z_H deberia estar nas faces??
!                    / Me%ExternalVar%HT(I,J)

                CMIST = CVK * (Z_H+ Me%ExternalVar%BottomRugosity) * (1.0 - DELTA * Z_H / Me%ExternalVar%HT(I,J))
                CMIST = CMIST * EXP(CNH * Aux)  
                Me%TurbVar%MixingLengthZ(I,J,K) = MIN(CMIST, Me%TurbVar%MAXMixingLength)


                VISC_V =(Me%TurbVar%MixingLengthZ(I,J,K) * Me%TurbVar%MixingLengthZ(I,J,K)) &
                       * SQRT(Me%TurbVar%FPRANDTL(I,J,K))
                
                Me%Viscosity%Vertical(I,J,K) = VISC_V

                 Me%Diffusivity%Vertical(i, j, k)      = VISC_V                        &
                                                                       *Me%TurbVar%VertPrandtlNumber(i,j,k)

!Manolo                Me%Viscosity%Vertical(I,J,K) = VISC_V + MolecularViscosity     !correccao pela visc. molec

            end if cd1
        end do do3
        end do do2
        !$OMP END DO
        end do do1
        
        !$OMP END PARALLEL
        
        if (MonitorPerformance) &
            call StopWatch ("ModuleTurbulence", "NihoulModel")

        !----------------------------------------------------------------------

    end subroutine NihoulModel


    !--------------------------------------------------------------------------
    ! Subroutine Richardson.
    ! 
    ! Computes Brunt-Vaisala and Prandtl frequencies and the Richardson number
    !
    ! Last modified:
    ! 2001: Manuel RV 
    !--------------------------------------------------------------------------

    subroutine Richardson

        !Local-----------------------------------------------------------------
        integer                                 :: ILB, IUB
        integer                                 :: JLB, JUB
        integer                                 :: KLB, KUB
        integer                                 :: i, j, k, kbottom
        real                                    :: U1, V1
        real                                    :: U2, V2
        real                                    :: VMODK1, VMODK2
        real                                    :: DRODZ, RO_PERT, RO, Depth
        real,    pointer, dimension(:,:,:)      :: DWZ
        real,    pointer, dimension(:,:,:)      :: DZZ
        real,    pointer, dimension(:,:,:)      :: ZCellCenter
        real,    pointer, dimension(:,:,:)      :: SZZ
        real,    dimension(:,:,:), pointer      :: VelocityX, VelocityY
        real,    dimension(:,:,:), pointer      :: S, T
        integer, dimension(:,:,:), pointer      :: ComputeFacesU3D, ComputeFacesV3D
        integer, dimension(:,:  ), pointer      :: KFloorZ
        real                                    :: RICH
        integer                                 :: CHUNK
        
        !----------------------------------------------------------------------

        if (MonitorPerformance) &
            call StartWatch ("ModuleTurbulence", "Richardson")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        DZZ             => Me%ExternalVar%DZZ
        DWZ             => Me%ExternalVar%DWZ
        SZZ             => Me%ExternalVar%SZZ
        ZCellCenter     => Me%ExternalVar%ZCellCenter
        KFloorZ         => Me%ExternalVar%KFloorZ
        T               => Me%ExternalVar%Temperature
        S               => Me%ExternalVar%Salinity
        VelocityX       => Me%ExternalVar%VelocityX
        VelocityY       => Me%ExternalVar%VelocityY
        ComputeFacesU3D => Me%ExternalVar%ComputeFacesU3D
        ComputeFacesV3D => Me%ExternalVar%ComputeFacesV3D        

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,RICH,DRODZ,RO,RO_PERT,U1,V1,U2,V2,Depth,VMODK1,VMODK2)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do J = JLB, JUB
do3 :   do I = ILB, IUB

cd1 :   if (Me%ExternalVar%WaterPoints3D(i, j, KUB)   == WaterPoint)         then

            kbottom = KFloorZ(i, j)
            
do1 :       do K = kbottom, KUB-1                         

                U1  =(VelocityX(I,  J+1,K  ) * ComputeFacesU3D(I,  J+1,K  )  + &
                      VelocityX(I,  J,  K  ) * ComputeFacesU3D(I,  J,  K  )) / 2.0

                V1  =(VelocityY(I+1,J,  K  ) * ComputeFacesV3D(I+1,J,  K  )  + &
                      VelocityY(I,  J,  K  ) * ComputeFacesV3D(I,  J,  K  )) / 2.0    

                U2  =(VelocityX(I,  J+1,K+1) * ComputeFacesU3D(I,  J+1,K+1)  + &
                      VelocityX(I,  J,  K+1) * ComputeFacesU3D(I,  J,  K+1)) / 2.0

                V2  =(VelocityY(I+1,J,  K+1) * ComputeFacesV3D(I+1,J,  K+1)  + &
                      VelocityY(I,  J,  K+1) * ComputeFacesV3D(I,  J,  K+1)) / 2.0

                select case (Me%TurbOptions%DensityMethod)                   

                    case (UNESCOState_)

                        RO      = SigmaUNESCO(T(i, j, k+1), S(i, j, k+1))
                        RO_PERT = SigmaUNESCO(T(i, j, k  ), S(i, j, k  ))

                        if (Me%TurbOptions%PressureCorrec) then

                            !calculo de Dro/Dz
                            Depth = -1.0*ZCellCenter(i,j,k+1)

                            RO      = SigmaUNESCOPressureCorrection (T(i, j, k+1), S(i, j, k+1), Depth, RO     )
                            RO_PERT = SigmaUNESCOPressureCorrection (T(I, J, K  ), S(i, j, k  ), Depth, RO_PERT)

                        end if

                    case (Mel96State_)

                            RO      = SigmaUNESCO(T(i, j, k+1), S(i, j, k+1))
                            RO_PERT = SigmaUNESCO(T(i, j, k  ), S(i, j, k  ))
                         
                            if (Me%TurbOptions%PressureCorrec) then
                            
                               !calculo de Dro/Dz
                               Depth = -1.0*ZCellCenter(i,j,k+1)

                               RO      = SigmaMel96PressureCorrection (T(i, j, k+1), S(i, j, k+1), Depth, RO      )
                               RO_PERT = SigmaMel96PressureCorrection (T(i, j, k  ), S(i, j, k  ), Depth, RO_PERT )
                        
                            end if

                    case (JMD95State_)

                        RO    = SigmaUNESCO(T(i, j, k+1), S(i, j, k+1))
                        RO_PERT = SigmaUNESCO(T(i, j, k), S(i, j, k))

                        if (Me%TurbOptions%PressureCorrec) then

                           !calculo de Dro/Dz
                           Depth = -1.0*ZCellCenter(i,j,k+1)

                           RO      = SigmaJMD95PressureCorrection  (T(i, j, k+1), S(i, j, k+1), Depth, RO)
                           RO_PERT = SigmaJMD95PressureCorrection  (T(i, j, k  ), S(i, j, k  ), Depth, RO_PERT)

                        end if

                    case (WangState_)

                        RO      = SigmaWang  (T(i, j, k+1), S(i, j, k+1))
                        RO_PERT = SigmaWang  (T(i, j, k  ), S(i, j, k  ))

                    case (LeendertseState_)

                        RO      = SigmaLeendertse  (T(i, j, k+1), S(i, j, k+1))
                        RO_PERT = SigmaLeendertse  (T(i, j, k  ), S(i, j, k  ))

                    case (Linear_)

                        RO      = 1025. -  dble(SigmaDensityReference) + 0.78 * (S(i, j, k+1) - 33.75)
                        RO_PERT = 1025. -  dble(SigmaDensityReference) + 0.78 * (S(i, j, k  ) - 33.75)
                                    

                    case(ConstantDensity_)
                        
                        !If density effects are not considered then DRODZ is null as well
                        !as the Brunt-Vaisalla frequency
                        RO      = SigmaDensityReference
                        RO_PERT = SigmaDensityReference 
                        !write(*,*)'WARNING: ConstantDensity_ method used in ModuleTurbulence'
                        !write(*,*)'This method is inconsistent if S and T are calculated in ModuleWaterProperties'

                    case default
                        
                        write(*,*)'Unknown method to compute density.'
                        stop 'Richardson - ModuleTurbulence - ERR00'

                 end select

                 DRODZ = (RO - RO_PERT) / DZZ(i,j,k)

!MAnolo          AUXDENS1 = (Density(I,J,K) * DWZ(I,J,K+1) + Density(I,J,K+1) * DWZ(I,J,K))/(DWZ(I,J,K+1)+DWZ(I,J,K))

                 !Calculo do quadrado da frequencia de Prandtl. 
                 ! Fprandtl in the bottom face has index Kbottom.
                 Me%TurbVar%FPRANDTL(I,J,K+1) = ((U2 - U1) / DZZ(I,J,K))**2.0 + &
                                                ((V2 - V1) / DZZ(I,J,K))**2.0
                                                         
                 !Calculo do quadrado da frequencia de Brunt-Vaisalla
                 Me%TurbVar%FBRUNTV(I,J,K+1) = -1.0 * G * DRODZ / SigmaDensityReference
                   
                 RICH = Me%TurbVar%FBRUNTV(I,J,K+1) / (Me%TurbVar%FPRANDTL(I,J,K+1)+1.e-10)

                 Me%TurbVar%Richardson(I,J,K+1) = MIN(RICH,RichardsonMAX)

                 if (Me%TurbOptions%MODTURB == backhaus_) then 
                     
                     VMODK1 = ABS(CMPLX(U1*U1,V1*V1))  
                     VMODK2 = ABS(CMPLX(U2*U2,V2*V2))  
                     !Calculo do modulo da velocidade
                     Me%TurbVar%VMOD(I,J,K) = (VMODK1 * DWZ(I,J,K+1) + &
                                               VMODK2 * DWZ(I,J,K )) / &
                                              (DWZ(I,J,K+1) + DWZ(I,J,K))
                 
                 end if
 
             end do do1
                
           !Necessary if turbulence equations are solved            
           Me%TurbVar%FBRUNTV (I,J,KBottom) = Me%TurbVar%FBRUNTV(I,J,KBottom+1)
           Me%TurbVar%FPRANDTL(I,J,KBottom) = Me%TurbVar%FPRANDTL(I,J,KBottom+1)

           Me%TurbVar%FBRUNTV (I,J,KUB+1)   = Me%TurbVar%FBRUNTV(I,J,KUB)
           Me%TurbVar%FPRANDTL(I,J,KUB+1)   = Me%TurbVar%FPRANDTL(I,J,KUB)
           
           RICH = Me%TurbVar%FBRUNTV(I,J,KBottom) / (Me%TurbVar%FPRANDTL(I,J,KBottom)+1.e-10)
           Me%TurbVar%Richardson(I,J,KBottom) = MIN(RICH,RichardsonMAX)

           RICH = Me%TurbVar%FBRUNTV(I,J,KUB+1) / (Me%TurbVar%FPRANDTL(I,J,KUB+1)+1.e-10)
           Me%TurbVar%Richardson(I,J,KUB+1)   = MIN(RICH,RichardsonMAX)

        end if cd1          

        end do do3
        end do do2
        !$OMP END DO
        !$OMP END PARALLEL

        nullify(DZZ    )
        nullify(DWZ    )
        nullify(KfloorZ)
        nullify(VelocityX   )
        nullify(VElocityY    )
        nullify(ComputeFacesU3D)
        nullify(ComputeFacesV3D) 

        if (MonitorPerformance) &
            call StopWatch ("ModuleTurbulence", "Richardson")

    end subroutine Richardson

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine TurbulenceEquationModel

        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer     :: DiffVertical, ViscVertical
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleTurbulence", "TurbulenceEquationModel")

        call TurbGOTM(Me%ObjTurbGOTM,                                &
                      Me%TurbVar%FBRUNTV,                            &
                      Me%TurbVar%FPRANDTL,                           &
                      STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceEquationModel - ModuleTurbulence - ERR01'
            
        call GetViscosityTurbGOTM(Me%ObjTurbGOTM, ViscVertical, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceEquationModel - ModuleTurbulence - ERR02'

        call SetMatrixValue(Me%Viscosity%Vertical, Me%Size, ViscVertical)

        call UngetTurbGOTM(Me%ObjTurbGOTM, ViscVertical, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceEquationModel - ModuleTurbulence - ERR03'
            
        call GetDiffusivityTurbGOTM(Me%ObjTurbGOTM, DiffVertical, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceEquationModel - ModuleTurbulence - ERR04'
        
        call SetMatrixValue(Me%Diffusivity%Vertical, Me%Size, DiffVertical)

        call UngetTurbGOTM(Me%ObjTurbGOTM, DiffVertical, STAT = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_) stop 'TurbulenceEquationModel - ModuleTurbulence - ERR05'
                  
        if (MonitorPerformance) call StopWatch ("ModuleTurbulence", "TurbulenceEquationModel")

    end subroutine TurbulenceEquationModel


    !-------------------------------------------------------------------------- 
    ! The depth of the mixed layer is computed in this subroutine
    !--------------------------------------------------------------------------

    subroutine ComputeMixedLayerDepth

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
!        real, dimension(:,:,:), pointer :: DiffVertical, ViscVertical
        integer :: ILB, IUB
        integer :: JLB, JUB
        integer :: KUB

        integer :: I, J, k, kbottom

        real, pointer, dimension(:,:,:) :: TKE
        real, pointer, dimension(:,:)   :: MLD_Surf, FBruntVmax
        real, pointer, dimension(:,:)   :: MLD_Bot
        integer, pointer, dimension(:,:):: KFloorZ

        integer :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) &
            call StartWatch ("ModuleTurbulence", "ComputeMixedLayerDepth")
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KUB = Me%WorkSize%KUB
        
        MLD_Surf => Me%TurbVar%MLD_Surf
        KFloorZ  => Me%ExternalVar%KFloorZ
        TKE      => Me%TurbVar%TKE
        
        if(Me%TurbOptions%MLD_Calc_Bot)                                       &
            MLD_Bot  => Me%TurbVar%MLD_Bot

        if (Me%TurbOptions%MLD_Method == brunt_mld_) then

            allocate (FBruntVmax(Me%Size%ILB:Me%Size%IUB,          &
                                 Me%Size%JLB:Me%Size%JUB))

            FBruntVmax(:,:) = FillValueReal

        endif

            
        !In the computation of mixed layer depth, care has to be taken on the indexing of szz
        !szz(0) is the water depth, so the distance from surface to face k(lower face of control volume k) is szz(k-1)

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,k,kbottom)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :   do J = JLB, JUB
do2 :   do I = ILB, IUB

            if (Me%ExternalVar%OpenPoints3D(i,j,KUB) == OpenPoint) then
                            
                MLD_Surf(i,j) = 0.
                k=KUB

                select case(Me%TurbOptions%MLD_Method)
      
                case(TKE_MLD_)          ! MLD according to TKE criterium
         
                    do while (TKE(i,j,k) > Me%TurbVar%TKE_MLD)

                        MLD_Surf(i,j) = Me%ExternalVar%SZZ(i,j,k-1)
                        k=k-1
                        if(k == KFloorZ(i, j)) then
                            MLD_Surf(i,j) = Me%ExternalVar%SZZ(i,j,k-1)
                            exit
                        end if

                    end do   

                    if(Me%TurbOptions%MLD_Calc_Bot) then

                        MLD_Bot(i,j) = 0.
                        k = KFloorZ(i,j)+1

                        do while (TKE(i,j,k) > Me%TurbVar%TKE_MLD)

                            MLD_Bot(i,j) = Me%ExternalVar%HT(i,j)-Me%ExternalVar%SZZ(i,j,k-1)
                            k=k+1
                            if(k == KUB+1) then
                                MLD_Bot(i,j) = Me%ExternalVar%HT(i,j)
                                exit
                            end if

                        end do   
                 
                    end if

                
                case(Rich_MLD_)         ! MLD according to critical Ri number
 
                    do while (Me%TurbVar%Richardson(i,j,k) < Me%TurbVar%Rich_MLD)

                        MLD_Surf (i,j) = Me%ExternalVar%SZZ(i,j,k-1)
                        k=k-1
                        if(k == KFloorZ(i, j)) then
                            MLD_Surf(i,j) = Me%ExternalVar%SZZ(i,j,k-1)
                            exit
                        end if

                    end do    
                 
                    if(Me%TurbOptions%MLD_Calc_Bot) then

                        MLD_Bot(i,j) = 0.
                        k = KFloorZ(i,j)+1

                        do while (Me%TurbVar%Richardson(i,j,k) < Me%TurbVar%Rich_MLD)

                            MLD_Bot(i,j) = Me%ExternalVar%HT(i,j)-Me%ExternalVar%SZZ(i,j,k-1)
                            k=k+1
                            if(k == KUB+1) then
                                MLD_Bot(i,j) = Me%ExternalVar%HT(i,j)
                                exit
                            end if

                        end do   
                 
                    end if           
                    
                case(brunt_mld_)         ! MLD according to maximum Brunt-vaisalla number
               
                    kbottom = KFloorZ(i, j)

                    do k = KUB, kbottom, -1
 
                        if (Me%TurbVar%FBRUNTV(i,j,k) >= FBruntVmax(i, j)) then
                 
                            FBruntVmax(i, j)     = Me%TurbVar%FBRUNTV(i,j,k)

                            MLD_Surf (i,j) = Me%ExternalVar%SZZ(i,j, k - 1) -          &
                                             Me%ExternalVar%SZZ(i,j, KUB  )  
                
                        endif

                    end do    
                 
                    if(Me%TurbOptions%MLD_Calc_Bot) then

                        call SetError (FATAL_, INTERNAL_, "ComputeMixedLayerDepth - Turbulence - ERR01")

                    end if      
                    
                end select
      
            end if
        
        end do do2
        end do do3
        !$OMP END DO
        !$OMP END PARALLEL

        if (Me%StatMLD%ON) then

            call OutPut_Statistics (MLD_Surf, Me%StatMLD%ID)


            if (Me%TurbOptions%MLD_Method == brunt_mld_)                      &
                call OutPut_Statistics (FBruntVmax, Me%StatMLD%FrBV_ID)

        endif
        
        if (Me%TurbOptions%MLD_Method == brunt_mld_) then
           
            deallocate (FBruntVmax)
            nullify    (FBruntVmax)

        endif

        nullify    (MLD_Surf  )
        nullify    (KFloorZ   )
        nullify    (TKE       )
        
        if(Me%TurbOptions%MLD_Calc_Bot)                                       &
            nullify    (MLD_Bot)

        if (MonitorPerformance) &
            call StopWatch ("ModuleTurbulence", "ComputeMixedLayerDepth")    
            
    !----------------------------------------------------------------------

    end subroutine ComputeMixedLayerDepth

    
    !--------------------------------------------------------------------------


    subroutine OutPut_Statistics (Value2D, StatisticsID)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer :: Value2D
        integer                       :: StatisticsID

        !Local-----------------------------------------------------------------
        integer                       :: MethodStatistic, Value2DStat2D
        integer                       :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetStatisticMethod (StatisticsID, MethodStatistic, STAT = STAT_CALL)                                     
        if (STAT_CALL /= SUCCESS_)                                                      & 
            call SetError (FATAL_, INTERNAL_, 'OutPut_Statistics - Turbulence - ERR01')
                                                                                    
        call GetStatisticParameters (StatisticsID,                                   &
                                     Value2DStat2D = Value2DStat2D,                  &
                                     STAT          = STAT_CALL)                        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, INTERNAL_, 'OutPut_Statistics - Turbulence - ERR02')
                                                                                    
                                                                                    
        if (MethodStatistic /= Value2DStat2D)                                        &
            call SetError (FATAL_, INTERNAL_, 'OutPut_Statistics - Turbulence - ERR03')


        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Statistics - Turbulence - ERR04'
                                                                                   
        call ModifyStatistic (StatisticsID,                                              &
                              Value2D       = Value2D,                                   &
                              WaterPoints2D = Me%ExternalVar%WaterPoints2D,              &
                              STAT          = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Statistics - Turbulence - ERR05'

        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Statistics - Turbulence - ERR06'



    end subroutine OutPut_Statistics

    
    !--------------------------------------------------------------------------


    subroutine OutPut_TimeSeries

        !Local-----------------------------------------------------------------
        real                                    :: DepthLevel
        integer                                 :: STAT_CALL, TimeSerieNumber, dn, id, jd, kd
        logical                                 :: DepthON, IgnoreOK

        !Begin-----------------------------------------------------------------

        !Corrects if necessary the cell of the time serie based in the time serie depth
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR10'  

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      LocalizationI = id,                               &
                                      LocalizationJ = jd,                               &
                                      DepthLevel    = DepthLevel,                       &
                                      DepthON       = DepthON,                          & 
                                      STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR20'

            if (DepthON) then

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR30'

                    if (IgnoreOK) then
                        cycle
                    else
                        stop 'OutPut_TimeSeries - ModuleTurbulence - ERR40'
                    endif

                endif

                kd = GetLayer4Level(Me%ObjGeometry, id, jd, DepthLevel, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR50'


                if (Me%ExternalVar%WaterPoints3D(id, jd, kd) == WaterPoint) then

                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn,  k = kd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR60'

                else
                    
                    stop 'OutPut_TimeSeries - ModuleTurbulence - ERR70'

                endif
            endif


        enddo


        !Eddy viscosity
        call WriteTimeSerie(Me%ObjTimeSerie,                                  &
                            Data3D = Me%Viscosity%Vertical,                   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR80'

        !Eddy diffusivity
        call WriteTimeSerie(Me%ObjTimeSerie,                                  &
                            Data3D = Me%Diffusivity%Vertical,                 &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR90'

        !Brunt-Vaisalla frequency
        call WriteTimeSerie(Me%ObjTimeSerie,                                  &
                            Data3D = Me%TurbVar%FBRUNTV,                      &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR100'

        !Prandtl frequency
        call WriteTimeSerie(Me%ObjTimeSerie,                                  &
                            Data3D = Me%TurbVar%FPRANDTL,                     &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR110'

        !Mixed layer depth
        if(Me%TurbOptions%MLD_Calc) then

          call WriteTimeSerie(Me%ObjTimeSerie,                                &
                            Data2D = Me%TurbVar%MLD_Surf,                     &
                            STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR120'

          if(Me%TurbOptions%MLD_Calc_Bot) then
            call WriteTimeSerie(Me%ObjTimeSerie,                              &
                            Data2D = Me%TurbVar%MLD_Bot,                      &
                            STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR130'
          end if
         
        
        end if
        
        if (Me%TurbOptions%MODTURB == TurbulenceEquation_) then

          !TKE
          call WriteTimeSerie(Me%ObjTimeSerie,                                &
                            Data3D = Me%TurbVar%TKE,                          &
                            STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR140'

          !Dissipation rate of TKE
          call WriteTimeSerie(Me%ObjTimeSerie,                                &
                            Data3D = Me%TurbVar%eps,                          &
                            STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR150'


          !Mixing length
          call WriteTimeSerie(Me%ObjTimeSerie,                                &
                            Data3D = Me%TurbVar%L,                            &
                            STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR160'

          !Production term in turbulence equations
          call WriteTimeSerie(Me%ObjTimeSerie,                                &
                            Data3D = Me%TurbVar%P,                            &
                            STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR170'


          !Buoyancy term in turbulence equations
          call WriteTimeSerie(Me%ObjTimeSerie,                                &
                            Data3D = Me%TurbVar%B,                            &
                            STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleTurbulence - ERR180'

        end if
 
   
    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------------------

    subroutine OutPut_Profile

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: i, j, k, kbottom
        integer                                 :: WILB, WIUB
        integer                                 :: WJLB, WJUB
        integer                                 :: WKUB

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleTurbulence", "OutPut_Profile")

        WILB = Me%WorkSize%ILB 
        WIUB = Me%WorkSize%IUB 
        WJLB = Me%WorkSize%JLB 
        WJUB = Me%WorkSize%JUB 
        WKUB = Me%WorkSize%KUB 

        call WriteProfile(Me%ObjProfile,                                        &
                          Me%Viscosity%HorizontalCenter,                        &
                          SZZ    = Me%ExternalVar%SZZ,                          &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'
        
        call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

        !ViscosityZ
        do j = WJLB, WJUB
        do i = WILB, WIUB
    
           if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
              kbottom = Me%ExternalVar%KFloorZ (i,j)

              do k = kbottom, WKUB

                 Me%OutPut%Aux3D(i,j,k) = (Me%Viscosity%Vertical(i,j,k)       &
                                          +Me%Viscosity%Vertical(i,j,k+1))/2.

                 Me%OutPut%Aux3D(i,j,k) = log10(abs(Me%OutPut%Aux3D(i,j,k))+1.e-14)

              end do
            
           end if
           
        end do
        end do     
        

        !ViscosityZ
        call WriteProfile(Me%ObjProfile,                                        &
                          Me%OutPut%Aux3D,                                      &
                          SZZ    = Me%ExternalVar%SZZ,                          &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'



        !Diffusivity
        call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

        do j=WJLB,WJUB
        do i=WILB,WIUB
     
            if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
                Kbottom = Me%ExternalVar%KFloorZ (i,j)

                do k=Kbottom,WKUB

                    Me%OutPut%Aux3D (i,j,k) = (Me%Diffusivity%Vertical(i,j,k)       &
                                              +Me%Diffusivity%Vertical(i,j,k+1))/2.

                    Me%OutPut%Aux3D(i,j,k) = log10(abs(Me%OutPut%Aux3D(i,j,k))+1.e-14)

                end do
            
            end if

        end do
        end do

        !Diffusivity
        call WriteProfile(Me%ObjProfile,                                        &
                          Me%OutPut%Aux3D,                                      &
                          SZZ    = Me%ExternalVar%SZZ,                          &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'

                
        !Brunt-Vaisalla frequency
        call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

        do j=WJLB,WJUB
        do i=WILB,WIUB
     
            if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
                Kbottom = Me%ExternalVar%KFloorZ (i,j)

                do k=Kbottom,WKUB

                    Me%OutPut%Aux3D(i,j,k) = (Me%TurbVar%FBruntv(i,j,k)        &
                                             +Me%TurbVar%FBruntv(i,j,k+1))/2.

                end do
            
            end if

        end do
        end do

        !Brunt-Vaisalla frequency
        call WriteProfile(Me%ObjProfile,                                        &
                          Me%OutPut%Aux3D,                                      &
                          SZZ    = Me%ExternalVar%SZZ,                          &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'

      

        !Prandtl frequency
        call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

        do j=WJLB,WJUB
        do i=WILB,WIUB
     
            if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
                Kbottom = Me%ExternalVar%KFloorZ (i,j)

                do k=Kbottom,WKUB

                    Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%FPrandtl(i,j,k)   &
                                              +Me%TurbVar%FPrandtl(i,j,k+1))/2.

                end do
            
            end if

        end do
        end do

        !Prandtl frequency 
        call WriteProfile(Me%ObjProfile,                                        &
                          Me%OutPut%Aux3D,                                      &
                          SZZ    = Me%ExternalVar%SZZ,                          &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'

                
   

        if (Me%TurbOptions%MODTURB == TurbulenceEquation_)then


            !TKE
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

            do j=WJLB,WJUB
            do i=WILB,WIUB
     
                if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
                    Kbottom = Me%ExternalVar%KFloorZ (i,j)

                    do k=Kbottom,WKUB

                        Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%TKE(i,j,k)   &
                                                  +Me%TurbVar%TKE(i,j,k+1))/2.

                        Me%OutPut%Aux3D(i,j,k)  = log10(abs(Me%OutPut%Aux3D(i,j,k))+1.e-14)

                    end do
            
                end if

            end do
            end do

            !TKE
            call WriteProfile(Me%ObjProfile,                                        &
                              Me%OutPut%Aux3D,                                      &
                              SZZ    = Me%ExternalVar%SZZ,                          &
                              STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'


            !eps
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

            do j=WJLB,WJUB
            do i=WILB,WIUB
     
                if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
                    Kbottom = Me%ExternalVar%KFloorZ (i,j)

                    do k=Kbottom,WKUB

                        Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%eps(i,j,k)   &
                                                  +Me%TurbVar%eps(i,j,k+1))/2.

                        Me%OutPut%Aux3D(i,j,k)  = log10(abs(Me%OutPut%Aux3D(i,j,k))+1.e-14)

                    end do
            
                end if

            end do
            end do

            !eps
            call WriteProfile(Me%ObjProfile,                                        &
                              Me%OutPut%Aux3D,                                      &
                              SZZ    = Me%ExternalVar%SZZ,                          &
                              STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'

            !L
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

            do j=WJLB,WJUB
            do i=WILB,WIUB
     
                if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
                    Kbottom = Me%ExternalVar%KFloorZ (i,j)

                    do k=Kbottom,WKUB

                        Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%L(i,j,k)   &
                                                  +Me%TurbVar%L(i,j,k+1))/2.

                        Me%OutPut%Aux3D(i,j,k)  = log10(abs(Me%OutPut%Aux3D(i,j,k))+1.e-14)

                    end do
            
                end if

            end do
            end do

            !L
            call WriteProfile(Me%ObjProfile,                                        &
                              Me%OutPut%Aux3D,                                      &
                              SZZ    = Me%ExternalVar%SZZ,                          &
                              STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'


            !P
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

            do j=WJLB,WJUB
            do i=WILB,WIUB
     
                if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
                    Kbottom = Me%ExternalVar%KFloorZ (i,j)

                    do k=Kbottom,WKUB

                        Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%P(i,j,k)   &
                                                  +Me%TurbVar%P(i,j,k+1))/2.

                    end do
            
                end if

            end do
            end do

            !P
            call WriteProfile(Me%ObjProfile,                                        &
                              Me%OutPut%Aux3D,                                      &
                              SZZ    = Me%ExternalVar%SZZ,                          &
                              STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'



            !B
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

            do j=WJLB,WJUB
            do i=WILB,WIUB
     
                if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint) then
           
                    Kbottom = Me%ExternalVar%KFloorZ (i,j)

                    do k=Kbottom,WKUB

                        Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%B(i,j,k)   &
                                                  +Me%TurbVar%B(i,j,k+1))/2.

                    end do

                end if

            end do
            end do

            !B
            call WriteProfile(Me%ObjProfile,                                        &
                              Me%OutPut%Aux3D,                                      &
                              SZZ    = Me%ExternalVar%SZZ,                          &
                              STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Profile - ModuleTurbulence - ERR01'
        end if


        if (MonitorPerformance) call StopWatch ("ModuleTurbulence", "OutPut_Profile")
    
    end subroutine OutPut_Profile
    
    !--------------------------------------------------------------------------------------

    subroutine Write_HDF5_Format

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL

        real, dimension(6), target          :: AuxTime
        real, dimension(:), pointer         :: TimePtr
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                             :: WorkKLB, WorkKUB, Index
        integer                             :: i, j, k, kbottom
        real                                :: TotalSeconds

        !Begin-------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleTurbulence", "Write_HDF5_Format")


        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 
        KLB = Me%Size%KLB 
        KUB = Me%Size%KUB 

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 


        !Current output index
        Index        = Me%OutPut%NextOutPut
        TotalSeconds = Me%ExternalVar%Now - Me%OutPut%OutTime(1)

        !Sets the attributes of the instant after HW
!        HouresAfterHW  = Me%OutPut%TimeAfterHW(Index)/3600.
!        MinutesAfterHW = (Me%OutPut%TimeAfterHW(Index) - (HouresAfterHW * 3600.0)) / 60.0

        !Writes current time
        call ExtractDate   (Me%ExternalVar%Now, AuxTime(1), AuxTime(2), AuxTime(3), &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime
        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR01'

        call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",&
                             Array1D = TimePtr, OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR02'



        !Writes SZZ
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,           &
                             WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR03'

        call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",       &
                             "m", Array3D = Me%ExternalVar%SZZ,               &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR04'


        !Writes OpenPoints
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                    &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB,              &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR05'
        
        call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",    &
                             "-", Array3D = Me%ExternalVar%OpenPoints3D,      &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR06'


        !Writes the horizontal turbulent viscosity compute for the center of the Z cells
        call HDF5WriteData  (Me%ObjHDF5, "/Results/ViscosityH",               &
                             "ViscosityH", "m2/s",                            &
                             Array3D = Me%Viscosity%HorizontalCenter,         &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR06a'


 
        !Viscosity. We have to interpolate it to the center of the cells because of the interface
        call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)


        do j=JLB,JUB
        do i=ILB,IUB
     
           if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
              Kbottom = Me%ExternalVar%KFloorZ (i,j)

              do k=Kbottom,WorkKUB

                 Me%OutPut%Aux3D (i,j,k) = (Me%Viscosity%Vertical(i,j,k)   &
                                           +Me%Viscosity%Vertical(i,j,k+1))/2.

              end do
            
           end if
           
        end do
        end do


        Me%OutPut%Aux3D(:,:,:) = log10(abs(Me%OutPut%Aux3D(:,:,:))+1.e-14)

        call HDF5WriteData  (Me%ObjHDF5, "/Results/ViscosityZ",                 &
                             "ViscosityZ", "log10(m2/s)",                       &
                             Array3D = Me%OutPut%Aux3D,                         &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR08'


        !Diffusivity
        call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)


        do j=JLB,JUB
        do i=ILB,IUB
     
           if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
              Kbottom = Me%ExternalVar%KFloorZ (i,j)

              do k=Kbottom,WorkKUB

                 Me%OutPut%Aux3D (i,j,k) = (Me%Diffusivity%Vertical(i,j,k)   &
                                           +Me%Diffusivity%Vertical(i,j,k+1))/2.

              end do
            
           end if
           
        end do
        end do

        Me%OutPut%Aux3D(:,:,:) = log10(abs(Me%OutPut%Aux3D(:,:,:))+1.e-14)

        call HDF5WriteData  (Me%ObjHDF5, "/Results/Diffusivity",                        &
                             "Diffusivity", "log10(m2/s)", Array3D = Me%OutPut%Aux3D,   &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR09'

        !NN
        call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

        do j=JLB,JUB
        do i=ILB,IUB
     
           if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
              Kbottom = Me%ExternalVar%KFloorZ (i,j)

              do k=Kbottom,WorkKUB

                 Me%OutPut%Aux3D(i,j,k) = (Me%TurbVar%FBruntv(i,j,k)        &
                                          +Me%TurbVar%FBruntv(i,j,k+1))/2.

              end do
            
           end if
           
        end do
        end do

!        Aux = log10(abs(Aux)+1.e-14)
        call HDF5WriteData  (Me%ObjHDF5, "/Results/NN",                         &
                             "NN", "s-2", Array3D = Me%OutPut%Aux3D,            &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR10'
        
        !SS
        call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

        do j=JLB,JUB
        do i=ILB,IUB
     
           if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
              Kbottom = Me%ExternalVar%KFloorZ (i,j)

              do k=Kbottom,WorkKUB

                 Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%FPrandtl(i,j,k)   &
                                           +Me%TurbVar%FPrandtl(i,j,k+1))/2.

              end do
            
           end if
           
        end do
        end do

!        Aux = log10(abs(Aux)+1.e-14)
        call HDF5WriteData  (Me%ObjHDF5, "/Results/SS",                         &
                             "SS", "s-2", Array3D = Me%OutPut%Aux3D,            &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR11'

ifTKE:  if (Me%TurbOptions%MODTURB .EQ. TurbulenceEquation_ ) then
        
            !TKE
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

            do j=JLB,JUB
            do i=ILB,IUB
     
               if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
                  Kbottom = Me%ExternalVar%KFloorZ (i,j)

                  do k=Kbottom,WorkKUB

                     Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%TKE(i,j,k)   &
                                               +Me%TurbVar%TKE(i,j,k+1))/2.

                  end do
            
               end if
           
            end do
            end do

            Me%OutPut%Aux3D(:,:,:) = log10(abs(Me%OutPut%Aux3D(:,:,:))+1.e-14)
            call HDF5WriteData  (Me%ObjHDF5, "/Results/TKE",                        &
                                 "TKE", "log10(m2/s2)", Array3D = Me%OutPut%Aux3D,  &
                                 OutputNumber = Index, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR12'

            !eps
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)
            
            do j=JLB,JUB
            do i=ILB,IUB
     
               if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
                  Kbottom = Me%ExternalVar%KFloorZ (i,j)

                  do k=Kbottom,WorkKUB

                     Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%eps(i,j,k)   &
                                                +Me%TurbVar%eps(i,j,k+1))/2.

                  end do
            
               end if
           
            end do
            end do

            Me%OutPut%Aux3D(:,:,:) = log10(abs(Me%OutPut%Aux3D)+1.e-14)
            call HDF5WriteData  (Me%ObjHDF5, "/Results/eps",                      &
                                 "eps", "log10(m2/s3)", Array3D = Me%OutPut%Aux3D,&
                                 OutputNumber = Index, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR13'
            
            !L
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)
            
            do j=JLB,JUB
            do i=ILB,IUB
     
               if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
                  Kbottom = Me%ExternalVar%KFloorZ (i,j)

                  do k=Kbottom,WorkKUB

                     Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%L(i,j,k)   &
                                                +Me%TurbVar%L(i,j,k+1))/2.

                  end do
            
               end if
           
            end do
            end do

            Me%OutPut%Aux3D(:,:,:) = log10(abs(Me%OutPut%Aux3D(:,:,:))+1.e-14)
            call HDF5WriteData  (Me%ObjHDF5, "/Results/L",                          &
                                 "L", "log10(m)", Array3D = Me%OutPut%Aux3D,        &
                                 OutputNumber = Index, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR14'

            !P
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)
            
            do j=JLB,JUB
            do i=ILB,IUB
     
               if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
                  Kbottom = Me%ExternalVar%KFloorZ (i,j)

                  do k=Kbottom,WorkKUB

                     Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%P(i,j,k)   &
                                                +Me%TurbVar%P(i,j,k+1))/2.

                  end do
            
               end if
           
            end do
            end do

    !        Aux = log10(abs(Aux)+1.e-14)
            call HDF5WriteData  (Me%ObjHDF5, "/Results/P",                        &
                                 "P", "m2/s3", Array3D = Me%OutPut%Aux3D,         &
                                 OutputNumber = Index, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR15'

            !B
            call SetMatrixValue(Me%OutPut%Aux3D, Me%Size, FillValueReal)

            do j=JLB,JUB
            do i=ILB,IUB
     
               if (Me%ExternalVar%WaterPoints3D(i, j, WorkKUB)   == WaterPoint) then
           
                  Kbottom = Me%ExternalVar%KFloorZ (i,j)

                  do k=Kbottom,WorkKUB

                     Me%OutPut%Aux3D (i,j,k) = (Me%TurbVar%B(i,j,k)   &
                                                +Me%TurbVar%B(i,j,k+1))/2.

                  end do
            
               end if
           
            end do
            end do

    !        Aux = log10(abs(Aux)+1.e-14)
            call HDF5WriteData  (Me%ObjHDF5, "/Results/B",                        &
                                 "B", "m2/s3", Array3D = Me%OutPut%Aux3D,         &
                                 OutputNumber = Index, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR16'
            
        endif ifTKE
                  
ifMLD:  if(Me%TurbOptions%MLD_Calc) then

            !MLD_Surf
            call HDF5WriteData  (Me%ObjHDF5, "/Results/MLD_Surf",                 &
                                 "MLD_Surf", "m", Array2D = Me%TurbVar%MLD_Surf,  &
                                 OutputNumber = Index, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR17'

            !MLD_bot
            if(Me%TurbOptions%MLD_Calc_Bot) then

                call HDF5WriteData  (Me%ObjHDF5, "/Results/MLD_Bot",               &
                                     "MLD_Bot", "m", Array2D = Me%TurbVar%MLD_Bot, &
                                     OutputNumber = Index, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR18'

            end if
      

        end if ifMLD
 
 
cd2:    if (Me%OutPut%Run_End) then

            call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR20'

        else

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleTurbulence - ERR21'

        endif cd2

        if (MonitorPerformance) call StopWatch ("ModuleTurbulence", "Write_HDF5_Format")

    end subroutine Write_HDF5_Format

    !--------------------------------------------------------------------------

    subroutine EstuaryModel

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB
        integer                             :: JLB, JUB
        integer                             :: KLB, KUB
        integer                             :: I, J, K
        integer, pointer, dimension(:,:,:) :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:) :: ComputeFacesV3D
        real                                :: VELMOD
        real                                :: U, V
        real                                :: VISHAUX
        real, pointer, dimension(:,:,:)     :: DWZ
        real, pointer, dimension(:,:,:)     :: VelocityX
        real, pointer, dimension(:,:,:)     :: VelocityY
        real, pointer, dimension(:,:  )     :: HT
        real                                :: MINHorizontalViscosity !Estuary & Smagorinsky
        real                                :: ReferenceDepth         !Estuary
        real                                :: ReferenceVelocity      !Estuary
        logical                             :: calc 

        integer                             :: CHUNK

        !----------------------------------------------------------------------


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        DWZ                     => Me%ExternalVar%DWZ
        VelocityX               => Me%ExternalVar%VelocityX
        VelocityY               => Me%ExternalVar%VelocityY
        ComputeFacesU3D         => Me%ExternalVar%ComputeFacesU3D
        ComputeFacesV3D         => Me%ExternalVar%ComputeFacesV3D
        HT                      => Me%ExternalVar%HT

        MINHorizontalViscosity  = Me%TurbVar%MINHorizontalViscosity
        ReferenceDepth          = Me%TurbVar%ReferenceDepth
        ReferenceVelocity       = Me%TurbVar%ReferenceVelocity

        if (MonitorPerformance) &
            call StartWatch ("ModuleTurbulence", "EstuaryModel")

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,K,calc,VELMOD,U,V,VISHAUX)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do J = JLB, JUB
do3 :   do I = ILB, IUB
            calc = .FALSE.

            !Vertical integration of SQRT(U**2+V**2)
            VELMOD = 0.0

do1 :       do K = KLB, KUB
cd1 :       if (Me%ExternalVar%WaterPoints3D(I,J,K) .EQ. WaterPoint) then

                calc = .TRUE.

                U      = 0.5*(VelocityX(I,J,K) * ComputeFacesU3D(I,J,K) + VelocityX(I,J+1,K) * ComputeFacesU3D(I,J+1,K))

                V      = 0.5*(VelocityY(I,J,K) * ComputeFacesV3D(I,J,K) + VelocityY(I+1,J,K) * ComputeFacesV3D(I+1,J,K))

                VELMOD = VELMOD + ABS(CMPLX(U*U,V*V)) * DWZ(I,J,K)
            end if cd1
            end do do1

cd3:        if (calc) then

                !V. media integrada na coluna de agua
                VELMOD = VELMOD / HT(I,J)

                VISHAUX = MINHorizontalViscosity * VELMOD * VELMOD * HT(I,J) / ReferenceDepth &
                        / (ReferenceVelocity * ReferenceVelocity)

                VISHAUX = MAX(VISHAUX,MINHorizontalViscosity)


do4 :           do K = KLB, KUB
cd2 :           if (Me%ExternalVar%WaterPoints3D(I,J,K) .EQ. WaterPoint) then

                    Me%Viscosity%HorizontalCenter(I,J,K) = VISHAUX

                end if cd2
                end do do4

            endif cd3

        end do do3
        end do do2
        !$OMP END DO
        !$OMP END PARALLEL
        
        nullify(DWZ)

        if (MonitorPerformance) &
            call StopWatch ("ModuleTurbulence", "EstuaryModel")

    end subroutine EstuaryModel

    !--------------------------------------------------------------------------

    subroutine SmagorinskyModel

        !Local-----------------------------------------------------------------

        integer :: ILB, IUB
        integer :: JLB, JUB
        integer :: KLB, KUB
        integer :: I, J, K

        integer, pointer, dimension(:,:,:) :: ComputeFacesU3D, ImposedTangentialFacesU, ImposedNormalFacesU
        integer, pointer, dimension(:,:,:) :: ComputeFacesV3D, ImposedTangentialFacesV, ImposedNormalFacesV

        real    :: AUX1,  AUX2
        real    :: UMED1, UMED2
        real    :: VMED1, VMED2
        real    :: dUdY,  dVdX
        real    :: dUdX,  dVdY
        real    :: DXDY
        real    :: ViscSmagorinsky

        real,    pointer, dimension(:,:  ) :: DUX
        real,    pointer, dimension(:,:  ) :: DVY
        real,    pointer, dimension(:,:,:) :: VelocityX
        real,    pointer, dimension(:,:,:) :: VelocityY
        integer                            :: FaceU1, FaceU2, FaceU3, FaceU4, FaceU5, FaceU6
        integer                            :: FaceV1, FaceV2, FaceV3, FaceV4, FaceV5, FaceV6

        integer :: CHUNK

        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB


        VelocityX               => Me%ExternalVar%VelocityX
        VelocityY               => Me%ExternalVar%VelocityY
        DUX                     => Me%ExternalVar%DUX
        DVY                     => Me%ExternalVar%DVY
        ComputeFacesU3D         => Me%ExternalVar%ComputeFacesU3D
        ComputeFacesV3D         => Me%ExternalVar%ComputeFacesV3D
        ImposedTangentialFacesU => Me%ExternalVar%ImposedTangentialFacesU
        ImposedNormalFacesU     => Me%ExternalVar%ImposedNormalFacesU
        ImposedTangentialFacesV => Me%ExternalVar%ImposedTangentialFacesV
        ImposedNormalFacesV     => Me%ExternalVar%ImposedNormalFacesV

        if (MonitorPerformance) &
            call StartWatch ("ModuleTurbulence", "SmagorinskyModel")

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,FaceU1,FaceU2,FaceU3,FaceU4,FaceU5,AUX1,AUX2,UMED1,UMED2,FaceV1,FaceV2,FaceV3,FaceV4,FaceV5,VMED1,VMED2,DXDY,dUdY,dVdX,dUdX,dVdY,ViscSmagorinsky)

do1 :   do K = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do J = JLB, JUB
do3 :   do I = ILB, IUB
cd1 :       if (Me%ExternalVar%WaterPoints3D(I,J,K) .EQ. WaterPoint) then
                ! Velocidade U na face de baixo "UMED1"

                FaceU1 = ImposedTangentialFacesU(I-1,J,  K) + ImposedNormalFacesU(I-1,J  ,K) + ComputeFacesU3D(I-1,J  ,K)
                if (FaceU1 > 0) FaceU1 = 1
                FaceU2 = ImposedTangentialFacesU(I-1,J+1,K) + ImposedNormalFacesU(I-1,J+1,K) + ComputeFacesU3D(I-1,J+1,K)
                if (FaceU2 > 0) FaceU2 = 1
                FaceU3 = ImposedTangentialFacesU(I  ,J,  K) + ImposedNormalFacesU(I  ,J,  K) + ComputeFacesU3D(I  ,J,  K)
                if (FaceU3 > 0) FaceU3 = 1
                FaceU4 = ImposedTangentialFacesU(I,  J+1,K) + ImposedNormalFacesU(I,  J+1,K) + ComputeFacesU3D(I,  J+1,K)
                if (FaceU4 > 0) FaceU4 = 1
                FaceU5 = ImposedTangentialFacesU(I+1,J,  K) + ImposedNormalFacesU(I+1,J,  K) + ComputeFacesU3D(I+1,J,  K)
                if (FaceU5 > 0) FaceU5 = 1
                FaceU6 = ImposedTangentialFacesU(I+1,J+1,K) + ImposedNormalFacesU(I+1,J+1,K) + ComputeFacesU3D(I+1,J+1,K)
                if (FaceU6 > 0) FaceU6 = 1

                AUX1 =(VelocityX(I-1,J,  K) * FaceU1 + VelocityX(I-1,J+1,K) * FaceU2) / 2.0

                AUX2 =(VelocityX(I,  J,  K) * FaceU3 + VelocityX(I,  J+1,K) * FaceU4) / 2.0

cd2 :           IF (I .EQ. ILB) THEN
                    UMED1 = AUX2 
                ELSE    
                    UMED1 = (AUX1 * DVY(I,J) + AUX2 * DVY(I-1,J)) / (DVY(I-1,J) + DVY(I,J))
                END IF cd2 



                ! Velocidade U na face de cima "UMED2"
                AUX1 =(VelocityX(I,  J,  K) * FaceU3 + VelocityX(I,  J+1,K) * FaceU4) / 2.0

                AUX2 =(VelocityX(I+1,J,  K) * FaceU5 + VelocityX(I+1,J+1,K) * FaceU6) / 2.0

cd3 :           IF (I .EQ. IUB) THEN
                    UMED2 = AUX1
                ELSE
                    UMED2 = (AUX1 * DVY(I+1,J) + AUX2 * DVY(I,J)) / (DVY(I+1,J) + DVY(I,J))
                END IF cd3

                FaceV1 = ImposedTangentialFacesV(I,  J-1,K) + ImposedNormalFacesV(I,  J-1,K) + ComputeFacesV3D(I,  J-1,K)
                if (FaceV1 > 0) FaceV1 = 1
                FaceV2 = ImposedTangentialFacesV(I+1,J-1,K) + ImposedNormalFacesV(I+1,J-1,K) + ComputeFacesV3D(I+1,J-1,K)
                if (FaceV2 > 0) FaceV2 = 1
                FaceV3 = ImposedTangentialFacesV(I  ,J,  K) + ImposedNormalFacesV(I  ,J,  K) + ComputeFacesV3D(I  ,J,  K)
                if (FaceV3 > 0) FaceV3 = 1
                FaceV4 = ImposedTangentialFacesV(I+1,J,  K) + ImposedNormalFacesV(I+1,J,  K) + ComputeFacesV3D(I+1,J,  K)
                if (FaceV4 > 0) FaceV4 = 1
                FaceV5 = ImposedTangentialFacesV(I,  J+1,K) + ImposedNormalFacesV(I,  J+1,K) + ComputeFacesV3D(I,  J+1,K)
                if (FaceV5 > 0) FaceV5 = 1
                FaceV6 = ImposedTangentialFacesV(I+1,J+1,K) + ImposedNormalFacesV(I+1,J+1,K) + ComputeFacesV3D(I+1,J+1,K)
                if (FaceV6 > 0) FaceV6 = 1

                ! Velocidade V na face esquerda "VMED1"
                AUX1 =(VelocityY(I,  J-1,K) * FaceV1 + VelocityY(I+1,J-1,K) * FaceV2) / 2.0

                AUX2 =(VelocityY(I,  J,  K) * FaceV3 + VelocityY(I+1,J,  K) * FaceV4) / 2.0

cd4 :           IF (J .EQ. JLB) THEN
                    VMED1 = AUX2
                ELSE
                    VMED1 = (AUX1 * DUX(I,J) + AUX2 * DUX(I,J-1)) / (DUX(I,J) + DUX(I,J-1))
                ENDIF cd4 

                ! Velocidade V na face da direita "VMED2"
                AUX1 =(VelocityY(I,  J,  K) * FaceV3 + VelocityY(I+1,J,  K) * FaceV4) / 2.0

                AUX2 =(VelocityY(I,  J+1,K) * FaceV5 + VelocityY(I+1,J+1,K) * FaceV6) / 2.0

cd5 :           IF (J .EQ. JUB) THEN
                    VMED2 = AUX1
                ELSE
                    VMED2 = (AUX1 * DUX(I,J+1) + AUX2 * DUX(I,J))/ (DUX(I,J) + DUX(I,J+1))
                END IF cd5


                DXDY = DUX(I,J) * DVY(I,J)
                dUdY = (UMED2 - UMED1) / DVY(I,J)
                dVdX = (VMED2 - VMED1) / DUX(I,J)

                dUdX =(VelocityX(I,J+1,K) * FaceU4 - VelocityX(I,J,K) * FaceU3) / DUX(I,J)

                dVdY =(VelocityY(I+1,J,K) * FaceV4 - VelocityY(I,J,K) * FaceV3) / DVY(I,J)


                ViscSmagorinsky = Me%TurbVar%HORCON * DXDY    &
                                * SQRT(dUdX * dUdX + 0.5 * (dVdX + dUdY) * (dVdX + dUdY) + dVdY * dVdY)

                Me%Viscosity%HorizontalCenter(I,J,K) =             &
                    MAX(Me%TurbVar%MINHorizontalViscosity, ViscSmagorinsky)

            end if cd1

        end do do3
        end do do2
        !$OMP END DO
        end do do1

        !$OMP END PARALLEL
        if (MonitorPerformance) &
            call StopWatch ("ModuleTurbulence", "SmagorinskyModel")
      
        nullify(VelocityX      )
        nullify(VelocityY      )
        nullify(DUX            )
        nullify(DVY            )
        nullify(ComputeFacesU3D, ImposedTangentialFacesU, ImposedNormalFacesU)
        nullify(ComputeFacesV3D, ImposedTangentialFacesV, ImposedNormalFacesV)


    end subroutine SmagorinskyModel


    !----------------------------------------------------------------------------
    !This subroutine computes the horizontal turbulent viscosity at Z cell corners

    subroutine TurbulentViscosity_CellCorner

        !Local-------------------------------------------------------------------
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k
        real                                    :: ViscFace1, ViscFace2
        integer, parameter                      :: Covered  = 1
        integer                                 :: WPT, WP1, WP2, WP3, WP4

        integer                                 :: CHUNK
        
        !------------------------------------------------------------------------  
        
        if (MonitorPerformance) &
            call StartWatch ("ModuleTurbulence", "TurbulentViscosity_CellCorner")
      
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        !Interpolates

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J,WP1,WP2,WP3,WP4,WPT,ViscFace1,ViscFace2)
        
do1 :   do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do j = JLB, JUB
do3 :   do i = ILB, IUB

            WP1 = Me%ExternalVar%WaterPoints3D(i  , j  , k)
            WP2 = Me%ExternalVar%WaterPoints3D(i  , j-1, k)
            WP3 = Me%ExternalVar%WaterPoints3D(i-1, j  , k)
            WP4 = Me%ExternalVar%WaterPoints3D(i-1, j-1, k)

            WPT = WP1 + WP2 + WP3 + WP4

            if (WPT > 0) then

                if (WPT == 4) then

                    ViscFace1 =(Me%Viscosity%HorizontalCenter(i,  j,  k) * Me%ExternalVar%DUX(i,  j-1)      &
                               +Me%Viscosity%HorizontalCenter(i,  j-1,k) * Me%ExternalVar%DUX(i,  j  ))     &
                              /(Me%ExternalVar%DUX(i,  j-1) + Me%ExternalVar%DUX(i,  j))


                    ViscFace2 =(Me%Viscosity%HorizontalCenter(i-1,j,  k) * Me%ExternalVar%DUX(i-1,j-1)      &
                               +Me%Viscosity%HorizontalCenter(i-1,j-1,k) * Me%ExternalVar%DUX(i-1,j  ))     &
                              /(Me%ExternalVar%DUX(i-1,j-1) + Me%ExternalVar%DUX(i-1,j))


                    Me%Viscosity%HorizontalCorner(i,j,k) =                              &
                                            (ViscFace1 * Me%ExternalVar%DYY(i-1,j)    + &
                                             ViscFace2 * Me%ExternalVar%DYY(i,j))     / &
                                            (Me%ExternalVar%DYY(i-1,j) + Me%ExternalVar%DYY(i,j))

                else

                    Me%Viscosity%HorizontalCorner(i,j,k) =                              &
                            (Me%Viscosity%HorizontalCenter(i,  j,  k) * WP1 +           &
                             Me%Viscosity%HorizontalCenter(i,  j-1,k) * WP2 +           &
                             Me%Viscosity%HorizontalCenter(i-1,j,  k) * WP3 +           &
                             Me%Viscosity%HorizontalCenter(i-1,j-1,k) * WP4) / WPT
                endif


            endif

        end do do3
        end do do2
        !$OMP END DO
        end do do1

        !$OMP END PARALLEL
        if (MonitorPerformance) call StopWatch ("ModuleTurbulence", "TurbulentViscosity_CellCorner")

    end subroutine TurbulentViscosity_CellCorner

    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillTurbulence(TurbulenceID, STAT)

        !Arguments-------------------------------------------------------------

        integer                         :: TurbulenceID
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------

        integer                            :: ready_              

        !Local-----------------------------------------------------------------

        integer                            :: STAT_, nUsers
        integer                        :: STAT_CALL

        !----------------------------------------------------------------------                         

       STAT_ = UNKNOWN_

        call Ready(TurbulenceID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTURBULENCE_,  Me%InstanceID)

            if (nUsers == 0) then

                if(Me%Output%ProfileON) then

                    call GetWaterPoints3D(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillTurbulence - ModuleTurbulence - ERR000'

                    call GetGeometryDistances (Me%ObjGeometry,                                  &
                                               SZZ  = Me%ExternalVar%SZZ,                       &
                                               STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillTurbulence - ModuleTurbulence - ERR001'
                    
                    
                    call GetGeometryKFloor(Me%ObjGeometry, Z = Me%ExternalVar%KFloorZ, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillTurbulence - ModuleTurbulence - ERR002'


                    if((Me%TurbOptions%MODTURB .eq. TurbulenceEquation_) .and.    &
                       ((Me%OutPut%ProfileON)))then
                             
                        call  GetTurbGOTM_TurbEq(Me%ObjTurbGOTM,                  &
                                                 Me%TurbVar%TKE,                  &
                                                 Me%TurbVar%eps,                  &
                                                 Me%TurbVar%L,                    &
                                                 Me%TurbVar%P,                    &
                                                 Me%TurbVar%B,                    &
                                                 STAT)
                    end if


                    call Output_Profile

                    if((Me%TurbOptions%MODTURB == TurbulenceEquation_).and.     &
                       ((Me%OutPut%ProfileON)))  then
                    
                            call UnGetTurbGOTM_TurbEq (Me%ObjTurbGOTM,          &
                                                       Me%TurbVar%TKE,          &
                                                       Me%TurbVar%eps,          &
                                                       Me%TurbVar%L,            &
                                                       Me%TurbVar%P,            &
                                                       Me%TurbVar%B,            &
                                                       STAT) 
                    end if


                    call UngetMap(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillTurbulence - ModuleTurbulence - ERR003'
                    
                    call UngetGeometry (Me%ObjGeometry, Me%ExternalVar%SZZ,     STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillTurbulence - ModuleTurbulence - ERR004'

                    call UngetGeometry (Me%ObjGeometry, Me%ExternalVar%KFloorZ, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillTurbulence - ModuleTurbulence - ERR005'

                end if
                

                !Deallocates Instance
                
cd7 :           if (Me%TurbOptions%MODTURB .EQ. TurbulenceEquation_ ) then
              
                    call KillTurbGOTM(Me%ObjTurbGOTM, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                &
                        stop 'Kill_Turbulence - ModuleTurbulence - ERR01'

                end if cd7


                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime           )
                if (nUsers == 0) stop 'Kill_Turbulence - ModuleTurbulence - ERR02'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap  )
                if (nUsers == 0) stop 'Kill_Turbulence - ModuleTurbulence - ERR03'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid )
                if (nUsers == 0) stop 'Kill_Turbulence - ModuleTurbulence - ERR04'

                nUsers = DeassociateInstance (mGEOMETRY_,       Me%ObjGeometry       )
                if (nUsers == 0) stop 'Kill_Turbulence - ModuleTurbulence - ERR05'

                nUsers = DeassociateInstance (mGRIDDATA_,       Me%ObjGridData       )
                if (nUsers == 0) stop 'Kill_Turbulence - ModuleTurbulence - ERR06'

                nUsers = DeassociateInstance (mMAP_,            Me%ObjMap            )
                if (nUsers == 0) stop 'Kill_Turbulence - ModuleTurbulence - ERR07'


                !Deallocates
                deallocate(Me%Viscosity%HorizontalCenter, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR09'


                deallocate(Me%Viscosity%Vertical,         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR10'


                deallocate(Me%Diffusivity%Vertical,         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR11'

                deallocate(Me%Viscosity%HorizontalCorner, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR12'


                deallocate(Me%TurbVar%Richardson,         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR13'


                deallocate(Me%TurbVar%FPRANDTL,            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR14'


                deallocate(Me%TurbVar%FBRUNTV,            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR15'


                deallocate(Me%TurbVar%MixingLengthX,       STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR16'


                deallocate(Me%TurbVar%MixingLengthY,       STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR17'


                deallocate(Me%TurbVar%MixingLengthZ,       STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR18'


                deallocate(Me%TurbVar%Ldownward,          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR19'


                if (Me%TurbOptions%MODTURB == backhaus_) then 

                   deallocate(Me%TurbVar%VMOD,               STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR20'
                
                end if

                deallocate(Me%TurbVar%VertPrandtlNumber,  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR21'

                if(Me%TurbOptions%MLD_Calc) then

                  deallocate(Me%TurbVar%MLD_Surf, STAT = STAT_CALL)
                  if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR22'

                  if (Me%StatMLD%ON) then

                    call KillStatistic (Me%StatMLD%ID, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                           &
                        stop 'Kill_Turbulence - ModuleTurbulence - ERR23'

                    if (Me%TurbOptions%MLD_Method == brunt_mld_) then  
                    
                        call KillStatistic (Me%StatMLD%FrBV_ID, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                       &
                            stop 'Kill_Turbulence - ModuleTurbulence - ERR24'

                    endif
 

                  endif
                  
                  if(Me%TurbOptions%MLD_Calc_Bot) then
                     deallocate(Me%TurbVar%MLD_Bot, STAT = STAT_CALL)
                     if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR25'
                  end if

                endif
                 

                if(Me%OutPut%ON .or. Me%OutPut%ProfileON)then
                    deallocate(Me%OutPut%Aux3D, STAT = STAT_CALL)
                     if (STAT_CALL /= SUCCESS_) stop 'Kill_Turbulence - ModuleTurbulence - ERR26'
                end if


                !Kills the TimeSerie
                if (Me%Output%TimeSerie) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                                 &
                        stop 'Kill_Turbulence - ModuleTurbulence - ERR27'
                endif
                
                
                if(Me%Output%ProfileON)then

                    call KillProfile(Me%ObjProfile, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                                 &
                        stop 'Kill_Turbulence - ModuleTurbulence - ERR28'

                end if


                !Nullifies
                nullify(Me%Viscosity%Vertical         )
                nullify(Me%Diffusivity%Vertical       )
                nullify(Me%Viscosity%HorizontalCenter )
                nullify(Me%Viscosity%HorizontalCorner )
      
                nullify(Me%TurbVar%Richardson         )
                nullify(Me%TurbVar%FPRANDTL           )           !Prandtl freq.
                nullify(Me%TurbVar%FBRUNTV            )           !Brunt-Vaisalla freq.
                nullify(Me%TurbVar%MixingLengthX      )
                nullify(Me%TurbVar%MixingLengthY      )
                nullify(Me%TurbVar%MixingLengthZ      )
                nullify(Me%TurbVar%Ldownward          )
                nullify(Me%TurbVar%VMOD               )
                nullify(Me%TurbVar%VertPrandtlNumber  )

                if(Me%TurbOptions%MLD_Calc) then
                   nullify(Me%TurbVar%MLD_surf)
                   if(Me%TurbOptions%MLD_Calc_Bot) nullify(Me%TurbVar%MLD_Bot)
                end if

                call DeallocateInstance

                TurbulenceID    = 0
                STAT_           = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
        !----------------------------------------------------------------------
    
    end subroutine KillTurbulence

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Turbulence), pointer          :: AuxObjTurbulence
        type (T_Turbulence), pointer          :: PreviousObjTurbulence

        !Updates pointers
        if (Me%InstanceID == FirstObjTurbulence%InstanceID) then
            FirstObjTurbulence => FirstObjTurbulence%Next
        else
            PreviousObjTurbulence => FirstObjTurbulence
            AuxObjTurbulence      => FirstObjTurbulence%Next
            do while (AuxObjTurbulence%InstanceID /= Me%InstanceID)
                PreviousObjTurbulence => AuxObjTurbulence
                AuxObjTurbulence      => AuxObjTurbulence%Next
            enddo

            !Now update linked list
            PreviousObjTurbulence%Next => AuxObjTurbulence%Next

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

    subroutine Ready(ObjTurbulence_Id, ready_) 

        !Arguments-------------------------------------------------------------

        integer :: ready_              !Auxiliar local variable

        integer :: ObjTurbulence_Id

        !----------------------------------------------------------------------
 
        nullify (Me)

cd1:    if (ObjTurbulence_Id > 0) then
            call LocateObjTurbulence (ObjTurbulence_Id)
            ready_ = VerifyReadLock (mTURBULENCE_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjTurbulence (ObjTurbulence_Id)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTurbulence_Id

        !Local-----------------------------------------------------------------

        Me => FirstObjTurbulence
        
        do while (associated (Me))
            if (Me%InstanceID == ObjTurbulence_Id) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleTurbulence - LocateObjTurbulence - ERR01'

    end subroutine LocateObjTurbulence

    !--------------------------------------------------------------------------

end module ModuleTurbulence

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
