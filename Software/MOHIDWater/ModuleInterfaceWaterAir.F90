!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : InterfaceWaterAir
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Pedro Chambel Leitao, Luis Fernandes- v4.0
! DESCRIPTION   : Module to compute water-air interface fluxes
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

Module ModuleInterfaceWaterAir

    use ModuleGlobalData
    use ModuleTime
    use ModuleHDF5             
    use ModuleFunctions,            only: SaturatedVaporPressure, ConstructPropertyID,          &
                                          LatentHeat, SensibleHeat, LongWaveDownward,           &
                                          LongWaveUpward, AerationFlux, AerationFlux_CO2,       &
                                          SetMatrixValue, CHUNK_J
    use ModuleEnterData,            only: ConstructEnterData, GetData, ExtractBlockFromBuffer,  &
                                          Block_Unlock, GetOutPutTime, ReadFileName,            &
                                          KillEnterData             
    use ModuleGridData,             only: GetGridData, UnGetGridData
    use ModuleHorizontalGrid,       only: GetHorizontalGridSize, WriteHorizontalGrid,           &
                                          RotateVectorFieldToGrid, RotateVectorGridToField,     &
                                          GetGridCellArea, UnGetHorizontalGrid, GetXYCellZ         
    use ModuleHorizontalMap,        only: GetWaterPoints2D, UnGetHorizontalMap          
    use ModuleGeometry,             only: GetGeometrySize
    use ModuleMap,                  only: GetWaterPoints3D, UnGetMap                   
    use ModuleBoxDif,               only: StartBoxDif, BoxDif, KillBoxDif               
    use ModuleTimeSerie,            only: StartTimeSerie, WriteTimeSerie, KillTimeSerie,        &
                                          GetTimeSerieLocation, CorrectsCellsTimeSerie,         &
                                          GetNumberOfTimeSeries, TryIgnoreTimeSerie
    use ModuleWaterProperties,      only: GetDensity, GetConcentration, UnGetWaterProperties,   &
                                          GetWaterPropertiesAirOptions, SetSurfaceFlux,         &
                                          GetPropertySurfaceFlux
    use ModuleHydrodynamic,         only: GetHydrodynamicAirOptions,SetAtmosphericPressure,     &
                                          SetSurfaceWaterFlux, SetWindStress,                   &
                                          GetHorizontalVelocity, UngetHydrodynamic

#ifndef _LAGRANGIAN_                                         
#ifdef  _LAGRANGIAN_GLOBAL_                                         
    use ModuleLagrangianGlobal,     only: SetLagrangianShearGlobal, GetLagrangianAirOptionsGlobal,   &
                                          SetLagrangianWindGlobal, SetLagSolarRadiationGlobal,  &
                                          SetLagrangianAtmPressureGlobal   
#else
    use ModuleLagrangian,           only: SetLagrangianShear, GetLagrangianAirOptions,          &
                                          SetLagrangianWind, SetLagrangianSolarRadiation,       &
                                          SetLagrangianAtmPressure
#endif    
#endif

    use ModuleTurbGOTM,             only: SetTurbGOTMSurfaceRugosity, SetTurbGOTMWindShearVelocity
    use ModuleAtmosphere,           only: GetAtmosphereProperty, AtmospherePropertyExists, UngetAtmosphere
    use ModuleFillMatrix,           only: ConstructFillMatrix, ModifyFillMatrix,                &
                                          GetDefaultValue, GetIfMatrixRemainsConstant, KillFillMatrix

#ifndef _WAVES_
    use ModuleWaves,                only: SetWavesWind, GetWaves, UnGetWaves
#endif
    use ModuleStopWatch,           only : StartWatch, StopWatch         


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartInterfaceWaterAir
    private ::      AllocateInstance
    private ::      ReadWaterAirFilesName
    private ::      ConstructGlobalVariables
    private ::      ConstructRugosity
    private ::      Construct_PropertyList
    private ::          Construct_Property
    private ::              Construct_PropertyValues
    private ::              Construct_PropertyEvolution
    private ::              Construct_PropertyOutPut
    private ::          Add_Property
    private ::      CheckInternalOptions
    private ::      CheckOptionsWater
    private ::      CheckOptionsAir
    private ::      Construct_Sub_Modules
    private ::          Construct_Time_Serie
    private ::          StartOutputBoxFluxes
    
    private ::      SetSubModulesConstructor

    !Selector
                     
    !Modifier
    public  :: ModifyInterfaceWaterAir
    private ::      ModifyRugosity
    private ::          ComputeWavesRugosity
    private ::      ModifyLocalAtmVariables
    private ::      ModifyWaterAirFluxes
    private ::          ModifySensibleHeat
    private ::              ComputeSensibleHeat
    private ::          ModifyLatentHeat
    private ::              ComputeLatentHeat
    private ::          ModifyEvaporation
    private ::              ComputeEvaporation
    private ::          ModifyNetLongWaveRadiation
    private ::              ComputeUpLongWaveRad
    private ::              ComputeDownLongWaveRad
    private ::          ModifyOxygenFlux
    private ::          ModifyCarbonDioxideFlux
    private ::              ModifyAerationFlux
    private ::              ModifyCO2AerationFlux
    private ::          ModifySurfaceRadiation
    private ::              ComputeSurfaceRadiation
    private ::          ModifyNonSolarFlux
    private ::              ComputeNonSolarFlux
    private ::          ModifyWindShearVelocity
    private ::          ModifyTurbulentKE
    private ::              ComputeTKEWind
    private ::          ModifySurfaceWaterFlux
    private ::      Output_TimeSeries
    private ::      Output_BoxTimeSeries  
    private ::      OutPut_Results_HDF

    !Destructor
    public  :: KillInterfaceWaterAir                                                     
    private ::      DeAllocateInstance


    !Management
    private ::      Ready
    private ::          LocateObjInterfaceWaterAir
    
    private ::              ReadLockExternalGlobal
    private ::              ReadLockExternalWater
    private ::              ReadUnlockExternalGlobal
    private ::              ReadUnlockExternalWater
    
    !Interfaces----------------------------------------------------------------
    
    !Parameters----------------------------------------------------------------
    character(LEN = StringLength), parameter        :: prop_block_begin = '<beginproperty>'
    character(LEN = StringLength), parameter        :: prop_block_end   = '<endproperty>'

    character(LEN = StringLength), parameter        :: rugosity_begin   = '<begin_rugosity>'
    character(LEN = StringLength), parameter        :: rugosity_end     = '<end_rugosity>'

    !Coefficients
    real,    parameter                              :: WaterEmissivity          = 0.96
    real,    parameter                              :: StefanBoltzmann          = 5.669e-08     ![W/m2/K4]
    real,    parameter                              :: BowenCoefficient         = 0.47          ![mmHg/ºC]
    real,    parameter                              :: ReferenceDensity         = 1000.         ![kg/m3]
    real,    parameter                              :: LatentHeatOfVaporization = 2.5e6         ![J/kg]

    !Constant for Air-Sea kinetic energy transfer        
    real,    parameter                              :: CDE = 0.63E-06        
    
    !Types---------------------------------------------------------------------
    private :: T_Files
    type       T_Files 
         character(len=PathLength)                  :: InputData
         character(len=PathLength)                  :: Results
         character(len=PathLength)                  :: BoxesFile
    end type T_Files

    private :: T_OutPut
    type       T_OutPut
         type (T_Time), pointer, dimension(:)       :: OutTime
         integer                                    :: NextOutPut
         logical                                    :: Yes                  =.false.
    end type T_OutPut
    
    private :: T_Ext_Global
    type       T_Ext_Global
        type(T_Time)                                :: Now
        real,    pointer, dimension(:,:)            :: GridCellArea
    end type T_Ext_Global


    private :: T_Ext_Options
    type       T_Ext_Options
        logical                                     :: HeatFluxYes                  = .false.
        logical                                     :: OxygenFluxYes                = .false.
        logical                                     :: CarbonDioxideFluxYes         = .false.       
        logical                                     :: WQMYes                       = .false.
        logical                                     :: SurfaceWaterFluxYes          = .false.
        logical                                     :: HydrodynamicWindYes          = .false.
        logical                                     :: HydrodynamicAtmPressureYes   = .false.
        logical                                     :: OilYes                       = .false.
        logical                                     :: WavesWindYes                 = .false.
        logical                                     :: GOTMWindShearVelocityYes     = .false.
        logical                                     :: Irrigation                   = .false.
        logical                                     :: Precipitation                = .false.
        logical                                     :: LagrangianWindYes            = .false.
        logical                                     :: LagrangianWQMYes             = .false.
        logical                                     :: LagrangianT90Yes             = .false.
        logical                                     :: T90VariableYes               = .false.
        logical                                     :: HydrodynamicMslpYes          = .false.
    end type    T_Ext_Options

    private :: T_Int_Options
    type       T_Int_Options
        logical                                     :: Evaporation                  = .false.
        logical                                     :: WindStress                   = .false.
        logical                                     :: SurfaceWaterFlux             = .false.
        logical                                     :: LatentHeat                   = .false.
        logical                                     :: SensibleHeat                 = .false.
        logical                                     :: NetLongWaveRadiation         = .false.
        logical                                     :: UpwardLongWaveRadiation      = .false.
        logical                                     :: DownwardLongWaveRadiation    = .false.
        logical                                     :: NonSolarFlux                 = .false.
        logical                                     :: OxygenFlux                   = .false.
        logical                                     :: CarbonDioxideFlux            = .false.
        logical                                     :: SpecificOxygenFlux           = .false.
        logical                                     :: SpecificCarbonDioxideFlux    = .false.        
        logical                                     :: WindShearVelocity            = .false.
        logical                                     :: TurbulentKineticEnergy       = .false.
        logical                                     :: SurfaceRadiation             = .false.
    end type   T_Int_Options

    
    private :: T_Ext_Water
    type       T_Ext_Water
        real,    pointer, dimension(:,:,:)          :: WaterTemperature
        real,    pointer, dimension(:,:,:)          :: WaterSalinity        
        real,    pointer, dimension(:,:,:)          :: WaterVelocity
        real,    pointer, dimension(:,:,:)          :: WaterDepth
        real,    pointer, dimension(:,:,:)          :: Density
        integer, pointer, dimension(:,:  )          :: WaterPoints2D
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D
        real,    pointer, dimension(:,:  )          :: Bathymetry
    end type T_Ext_Water

    private :: T_ExtField
    type       T_ExtField
        real,    pointer, dimension(:,:  )          :: Field
        logical                                     :: Yes                  = .false.      
    end type  T_ExtField


    private :: T_Ext_Atm
    type       T_Ext_Atm
        type(T_ExtField)                            :: AtmosphericPressure
        type(T_ExtField)                            :: Mslp     !Mean Sea Level Pressure
        type(T_ExtField)                            :: SolarRadiation
        type(T_ExtField)                            :: RelativeHumidity     
        type(T_ExtField)                            :: AirTemperature
        type(T_ExtField)                            :: CloudCover
        type(T_ExtField)                            :: WindVelocityX
        type(T_ExtField)                            :: WindVelocityY
        type(T_ExtField)                            :: WindDirection
    end type T_Ext_Atm

    private :: T_Local_Atm
    type       T_Local_Atm
        type(T_ExtField)                            :: AtmosphericPressure
        type(T_ExtField)                            :: Mslp     !Mean Sea Level Pressure
        type(T_ExtField)                            :: WindVelocityX
        type(T_ExtField)                            :: WindVelocityY
        type(T_ExtField)                            :: WindDirection
    end type   T_Local_Atm

    
    private :: T_Evolution
        type   T_Evolution
        logical                                     :: Variable             = .false.
        real                                        :: DTInterval
        type(T_Time)                                :: LastCompute
        type(T_Time)                                :: NextCompute
    end type T_Evolution

    private :: T_Property
    type       T_Property
         type(T_PropertyID)                         :: ID
         real, dimension(:,:), pointer              :: Field, FieldGrid
         type(T_Evolution)                          :: Evolution
         logical                                    :: TimeSerie            = .false.
         logical                                    :: BoxTimeSerie         = .false.
         logical                                    :: CEQUALW2             = .false.
         logical                                    :: OutputHDF            = .false.
         logical                                    :: Constant             = .false.
         type(T_Property), pointer                  :: Next
         type(T_Property), pointer                  :: Prev
    end type T_Property

    private :: T_Coupling
    type       T_Coupling
         type(T_Time)                               :: NextCompute
         real                                       :: DT_Compute           = FillValueReal
         logical                                    :: Yes                  = .false.
         integer                                    :: NumberOfProperties   = 0
    end type T_Coupling  

    private :: T_Coupled
    type       T_Coupled
         type(T_Coupling)                           :: TimeSerie
         type(T_Coupling)                           :: BoxTimeSerie
         type(T_Coupling)                           :: CEQUALW2
    end type T_Coupled

    type       T_Rugosity
        type(T_PropertyID)                          :: ID
        real, pointer, dimension (:,:)              :: Field 
        real                                        :: Scalar, WavesRelation
        logical                                     :: Constant, ON, WavesFunction
    end type T_Rugosity
    
    private :: T_InterfaceWaterAir
    type       T_InterfaceWaterAir
        integer                                     :: InstanceID
        character(PathLength)                       :: ModelName
        type(T_Time       )                         :: BeginTime
        type(T_Time       )                         :: EndTime
        type(T_Time       )                         :: ActualTime
        type(T_Size2D     )                         :: Size2D, WorkSize2D
        type(T_Size3D     )                         :: Size3D, WorkSize3D
        type(T_Files      )                         :: Files
        type(T_OutPut     )                         :: Output
        type(T_Property   ), pointer                :: FirstProperty
        type(T_Property   ), pointer                :: LastProperty
        type(T_Coupled    )                         :: Coupled
        type(T_Ext_Global )                         :: ExternalVar
        type(T_Ext_Water  )                         :: ExtWater
        type(T_Ext_Atm    )                         :: ExtAtm
        type(T_Ext_Options)                         :: ExtOptions
        type(T_Int_Options)                         :: IntOptions
        type(T_Local_Atm  )                         :: LocalAtm

        type(T_Rugosity   )                         :: Rugosity                 
        real                                        :: ReflectionCoef           = FillValueReal
        real                                        :: CDWIND                   = FillValueReal
        logical                                     :: DefineCDWIND             = .false.
        real(8), pointer, dimension(:,:)            :: Scalar2D
        real   , pointer, dimension(:,:)            :: WindShearVelocity
        integer                                     :: AerationEquation         = FillValueInt
        integer                                     :: CO2AerationEquation      = FillValueInt
        real                                        :: Altitude                 = FillValueReal
        real                                        :: AltitudeCorrection       = FillValueReal

        integer                                     :: PropertiesNumber         = 0

        !Instance of ModuleBoxDif                   
        integer                                     :: ObjBoxDif                = 0
                                                                            
        !Instance of ModuleEnterData                                        
        integer                                     :: ObjEnterData             = 0
                                                                            
        !Instance of ModuleTimeSerie                                        
        integer                                     :: ObjTimeSerie             = 0
                                                                                
        !Instance of ModuleHDF5                                                 
        integer                                     :: ObjHDF5                  = 0
                                                                                
        !Instance of ModuleTime                                                 
        integer                                     :: ObjTime                  = 0
                                                                                
        !Instance of ModuleHorizontalGrid                                       
        integer                                     :: ObjHorizontalGrid        = 0
                                                                                
        !Instance of ModuleGridData                                             
        integer                                     :: ObjGridData              = 0
                                                                                
        !Instance of ModuleGeometry                                             
        integer                                     :: ObjGeometry              = 0
        
        !Instance of ModuleHorizontalMap                                    
        integer                                     :: ObjHorizontalMap         = 0
                                                                                
        !Instance of ModuleMap
        integer                                     :: ObjMap                   = 0
        
        !Instance of ModuleHydrodynamic
        integer                                     :: ObjHydrodynamic          = 0

        !Instance of ModuleTurbGOTM
        integer                                     :: ObjTurbGOTM              = 0

        !Instance of ModuleWaves
        integer                                     :: ObjWaves                 = 0
        
        !Instance of ModuleInterfaceWaterAir
        integer                                     :: ObjWaterProperties       = 0

        !Instance of ModuleAtmosphere
        integer                                     :: ObjAtmosphere            = 0
        
        !Instance of ModuleLagrangian 
        integer                                     :: ObjLagrangian            = 0
        
        
        type(T_InterfaceWaterAir), pointer          :: Next
    end type  T_InterfaceWaterAir

    !Global Module Variables
    type (T_InterfaceWaterAir), pointer             :: FirstObjInterfaceWaterAir
    type (T_InterfaceWaterAir), pointer             :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartInterfaceWaterAir(ModelName,                      &
                                      ObjInterfaceWaterAirID,         &
                                      TimeID,                         &
                                      HorizontalGridID,               &
                                      WaterGridDataID,                &
                                      WaterHorizontalMapID,           &
                                      WaterMapID,                     &
                                      WaterGeometryID,                &
                                      HydrodynamicID,                 &
                                      TurbGOTMID,                     &
                                      WavesID,                        &
                                      WaterPropertiesID,              &
                                      AtmosphereID,                   &
                                      LagrangianID,                   &
                                      STAT)

        !Arguments---------------------------------------------------------------
        character(Len=*)                                :: ModelName
        integer                                         :: ObjInterfaceWaterAirID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: WaterGridDataID
        integer                                         :: WaterHorizontalMapID
        integer                                         :: WaterMapID
        integer                                         :: WaterGeometryID
        integer                                         :: HydrodynamicID        
        integer                                         :: TurbGOTMID
        integer                                         :: WavesID           
        integer                                         :: WaterPropertiesID
        integer                                         :: AtmosphereID
        integer                                         :: LagrangianID
        integer, optional, intent(OUT)                  :: STAT    

        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL       

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mInterfaceWaterAir_)) then
            nullify (FirstObjInterfaceWaterAir)
            call RegisterModule (mInterfaceWaterAir_) 
        endif

        call Ready(ObjInterfaceWaterAirID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ModelName = ModelName

            nullify (Me%FirstProperty)
            nullify (Me%LastProperty )

            nullify (Me%LocalAtm%WindVelocityX%Field)
            nullify (Me%LocalAtm%WindVelocityY%Field)
            nullify (Me%LocalAtm%WindDirection%Field)
            nullify (Me%LocalAtm%AtmosphericPressure%Field)
            nullify (Me%LocalAtm%Mslp%Field)    !Mean Sea Level Pressure
            
            !Associates External Instances
            Me%ObjTime                 = AssociateInstance(mTIME_,           TimeID                     )
            Me%ObjHorizontalGrid       = AssociateInstance(mHORIZONTALGRID_, HorizontalGridID           )
            
            !Water column
            Me%ObjGridData             = AssociateInstance(mGRIDDATA_,          WaterGridDataID         )
            Me%ObjHorizontalMap        = AssociateInstance(mHORIZONTALMAP_,     WaterHorizontalMapID    )
            Me%ObjGeometry             = AssociateInstance(mGEOMETRY_,          WaterGeometryID         )
            Me%ObjMap                  = AssociateInstance(mMAP_,               WaterMapID              )
            Me%ObjHydrodynamic         = AssociateInstance(mHYDRODYNAMIC_,      HydrodynamicID          )
            Me%ObjWaterProperties      = AssociateInstance(mWATERPROPERTIES_,   WaterPropertiesID       )
            Me%ObjAtmosphere           = AssociateInstance(mATMOSPHERE_,        AtmosphereID            )

            if(LagrangianID /= 0)then
                Me%ObjLagrangian       = AssociateInstance(mLAGRANGIAN_,        LagrangianID            )
            endif

            if(TurbGOTMID /= 0)then
                Me%ObjTurbGOTM         = AssociateInstance(mTURBGOTM_,          TurbGOTMID              )
            endif

            if(WavesID /= 0)then
                Me%ObjWaves            = AssociateInstance(mWAVES_,             WavesID                 )
            endif



            call ReadLockExternalGlobal
            
            call ReadLockExternalWater

            call ReadWaterAirFilesName

            call ConstructEnterData(Me%ObjEnterData, Me%Files%InPutData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                 &
                stop 'StartInterfaceWaterAir - InterfaceWaterAir - ERR01'

            call ConstructGlobalVariables

            call ConstructRugosity

            call Construct_PropertyList

            call CheckInternalOptions

            call CheckOptionsWater
            
            call CheckOptionsAir

            call Construct_Sub_Modules

            call ConstructGlobalOutput
            
            if (Me%OutPut%Yes) call Open_HDF5_OutPut_File

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                 &
                stop 'StartInterfaceWaterAir - InterfaceWaterAir - ERR02'

            call ReadUnlockExternalWater

            call ReadUnlockExternalGlobal

            call SetSubModulesConstructor

            call SetSubModulesModifier

            !Returns ID
            ObjInterfaceWaterAirID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleInterfaceWaterAir - StartInterfaceWaterAir - ERR99' 

        endif cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartInterfaceWaterAir
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance
                                                    
        !Local-----------------------------------------------------------------
        type (T_InterfaceWaterAir), pointer    :: NewObjWaterSedInterface
        type (T_InterfaceWaterAir), pointer    :: PreviousObjWaterSedInterface


        !Allocates new instance
        allocate (NewObjWaterSedInterface)
        nullify  (NewObjWaterSedInterface%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjInterfaceWaterAir)) then
            FirstObjInterfaceWaterAir           => NewObjWaterSedInterface
            Me                                  => NewObjWaterSedInterface
        else
            PreviousObjWaterSedInterface        => FirstObjInterfaceWaterAir
            Me                                  => FirstObjInterfaceWaterAir%Next
            do while (associated(Me))
                PreviousObjWaterSedInterface    => Me
                Me                              => Me%Next
            enddo
            Me                                  => NewObjWaterSedInterface
            PreviousObjWaterSedInterface%Next   => NewObjWaterSedInterface
        endif

        Me%InstanceID = RegisterNewInstance (mINTERFACEWATERAIR_)


    end subroutine AllocateInstance
    
    
   
    !--------------------------------------------------------------------------

    
    subroutine ReadWaterAirFilesName

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        character(len = StringLength)       :: Message

        !----------------------------------------------------------------------

        Message ='ASCII file used to construct new water-air properties.'
        Message = trim(Message)

        call ReadFileName('AIRW_DAT', Me%Files%InPutData, Message = Message, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadWaterAirFilesName - ModuleInterfaceWaterAir - ERR01'

        Message   ='Instant fields of water-air properties in HDF format.'
        Message   = trim(Message)

        call ReadFileName('AIRW_HDF', Me%Files%Results, Message = Message, TIME_END = Me%EndTime, &
                           Extension = 'arw', STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadWaterAirFilesName - ModuleInterfaceWaterAir - ERR02'


    end subroutine ReadWaterAirFilesName
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructGlobalVariables
                                                    
        !External--------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        call GetHorizontalGridSize(Me%ObjHorizontalGrid,                        &
                                   Size        = Me%Size2D,                     &
                                   WorkSize    = Me%WorkSize2D,                 &
                                   STAT        = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR10'

        call GetGeometrySize(Me%ObjGeometry,                                    &
                             Size       = Me%Size3D,                            &
                             WorkSize   = Me%WorkSize3D,                        &
                             STAT       = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR20'

        call GetComputeTimeLimits(Me%ObjTime,                                   &
                                  BeginTime = Me%BeginTime,                     &
                                  EndTime   = Me%EndTime,                       &
                                  STAT      = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR30'

        !Actualize the time
        Me%ActualTime = Me%BeginTime


        !Altitude in km
        call GetData(Me%Altitude,                                           &
                     Me%ObjEnterData, iflag,                                &
                     Keyword      ='ALTITUDE',                              &
                     SearchType   = FromFile,                               &
                     ClientModule = 'ModuleInterfaceWaterAir',              &
                     Default      = 0.0,                                    &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                         &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR40'

        !Mortimer's altitude correction
        Me%AltitudeCorrection = (1. - (Me%Altitude/1000.)/44.3)**5.25

        call GetData(Me%AerationEquation,                                   &
                     Me%ObjEnterData, iflag,                                &
                     Keyword      ='AERATION_METHOD',                       &
                     SearchType   = FromFile,                               &
                     ClientModule = 'ModuleInterfaceWaterAir',              &
                     Default      = Gelda_et_al_1996,                       &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                         &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR50'
        
        call GetData(Me%CO2AerationEquation,                                &
                     Me%ObjEnterData, iflag,                                &
                     Keyword      ='CO2_AERATION_METHOD',                   &
                     SearchType   = FromFile,                               &
                     ClientModule = 'ModuleInterfaceWaterAir',              &
                     Default      = Borges_et_al_2004,                      &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                         &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR60'
        

    end subroutine ConstructGlobalVariables
    
    !--------------------------------------------------------------------------

    subroutine ConstructRugosity

        !Local------------------------------------------------------------------
        integer                        :: ClientNumber, ILB, IUB, JLB, JUB 
        integer                        :: STAT_CALL, iflag
        logical                        :: BlockFound
        
        !Begin------------------------------------------------------------------

        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        Me%Rugosity%ON            = .false.

        call GetData(Me%Rugosity%Scalar,                                &
                     Me%ObjEnterData, iflag,                            &
                     keyword      ='RUGOSITY',                          &
                     SearchType   = FromFile,                           &
                     ClientModule = 'ModuleInterfaceWaterAir',          &
                     Default      = 0.0025,                             &
                     STAT         = STAT_CALL)

        if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR05'

i2:     if (iflag == 1) then

            nullify (Me%Rugosity%Field)
            allocate(Me%Rugosity%Field(ILB:IUB, JLB:JUB))
            Me%Rugosity%Field(:,:) = Me%Rugosity%Scalar
            Me%Rugosity%Constant      = .true.
            Me%Rugosity%ON            = .true.
            Me%Rugosity%WavesFunction = .false.
        
        else i2

            !Reads the Atmosphere rugosity
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                      &
                                        rugosity_begin, rugosity_end, BlockFound,           &
                                        STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR10'

i1:         if(BlockFound)then

                nullify (Me%Rugosity%Field)
                allocate(Me%Rugosity%Field(ILB:IUB, JLB:JUB))
                Me%Rugosity%Field(:,:) = FillValueReal

                Me%Rugosity%ON = .true.

                call ConstructFillMatrix  (PropertyID           = Me%Rugosity%ID,           &
                                           EnterDataID          = Me%ObjEnterData,          &
                                           TimeID               = Me%ObjTime,               &
                                           HorizontalGridID     = Me%ObjHorizontalGrid,     &
                                           ExtractType          = FromBlock,                &
                                           PointsToFill2D       = Me%ExtWater%WaterPoints2D,&
                                           Matrix2D             = Me%Rugosity%Field,        &
                                           TypeZUV              = TypeZ_,                   &
                                           STAT                 = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR20'


                call GetDefaultValue(Me%Rugosity%ID%ObjFillMatrix, Me%Rugosity%Scalar, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR30'

                call GetIfMatrixRemainsConstant(FillMatrixID    = Me%Rugosity%ID%ObjFillMatrix,&
                                                RemainsConstant = Me%Rugosity%Constant,        &
                                                STAT            = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR40'


                call GetData(Me%Rugosity%WavesFunction,                             &
                             Me%ObjEnterData, iflag,                                &
                             keyword      ='WAVES_FUNCTION',                        &
                             SearchType   = FromBlock,                              &
                             ClientModule = 'ModuleInterfaceWaterAir',              &
                             Default      = .false.,                                &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR60'

                if (Me%Rugosity%WavesFunction) then

                    call GetData(Me%Rugosity%WavesRelation,                         &
                                 Me%ObjEnterData, iflag,                            &
                                 keyword      ='WAVES_RELATION',                    &
                                 SearchType   = FromBlock,                          &
                                 ClientModule = 'ModuleInterfaceWaterAir',          &
                                 Default      = 1.,                                 &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR70'

                    if (Me%ObjWaves == 0) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR80'

                    if (Me%Rugosity%Constant) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR90'

                    call ModifyRugosity
            
                endif

                if(.not. Me%Rugosity%ID%SolutionFromFile) then

                    call KillFillMatrix(Me%Rugosity%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR50'
                endif
          

            else i1

                Me%Rugosity%ON = .false.

            endif i1

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR100'

        endif i2


    end subroutine ConstructRugosity

    !--------------------------------------------------------------------------
    
    subroutine Construct_PropertyList

        !External----------------------------------------------------------------
        integer                                 :: ClientNumber
        integer                                 :: STAT_CALL
        logical                                 :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_Property), pointer              :: NewProperty

        !------------------------------------------------------------------------

        !Initialize the properties number   
        Me%PropertiesNumber = 0


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = prop_block_begin,     &
                                        block_end       = prop_block_end,       &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Property 
                    Call Construct_Property(NewProperty)

                    ! Add new Property to the WaterProperties List 
                    Call Add_Property(NewProperty)
                else cd2
                    call Block_Unlock(Me%ObjEnterData , ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Construct_PropertyList - ModuleInterfaceWaterAir - ERR01'

                    exit do1    !No more blocks
                endif cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Construct_PropertyList - ModuleInterfaceWaterAir - ERR02'
            else cd1
                stop 'Construct_PropertyList - ModuleInterfaceWaterAir - ERR03'
            endif cd1
        end do do1

        !------------------------------------------------------------------------

    end subroutine Construct_PropertyList

    !----------------------------------------------------------------------------


    subroutine Construct_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        integer                         :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB 
        IUB = Me%Size2D%IUB 
        JLB = Me%Size2D%JLB 
        JUB = Me%Size2D%JUB 
             
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_Property - ModuleInterfaceWaterAir - ERR01' 

        nullify(NewProperty%Prev,NewProperty%Next)

        call ConstructPropertyID        (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call Construct_PropertyEvolution(NewProperty)

        call Construct_PropertyValues   (NewProperty)

        call Construct_PropertyOutPut   (NewProperty)

    end subroutine Construct_Property

    
    !--------------------------------------------------------------------------


    subroutine Construct_PropertyValues(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),   pointer                 :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------
        integer                                     :: SizeILB, SizeIUB
        integer                                     :: SizeJLB, SizeJUB 
        
        !----------------------------------------------------------------------

        SizeILB = Me%Size2D%ILB
        SizeIUB = Me%Size2D%IUB
        SizeJLB = Me%Size2D%JLB
        SizeJUB = Me%Size2D%JUB

        allocate (NewProperty%Field (SizeILB:SizeIUB, SizeJLB:SizeJUB))

        NewProperty%Field(:,:) = FillValueReal

        call ConstructFillMatrix  (PropertyID           = NewProperty%ID,           &
                                   EnterDataID          = Me%ObjEnterData,          &
                                   TimeID               = Me%ObjTime,               &
                                   HorizontalGridID     = Me%ObjHorizontalGrid,     &
                                   ExtractType          = FromBlock,                &
                                   PointsToFill2D       = Me%ExtWater%WaterPoints2D,&
                                   Matrix2D             = NewProperty%Field,        &
                                   TypeZUV              = TypeZ_,                   &
                                   STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR01'

        call GetIfMatrixRemainsConstant(FillMatrixID    = NewProperty%ID%ObjFillMatrix,&
                                        RemainsConstant = NewProperty%Constant,        &
                                        STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR02'

        if(.not. NewProperty%ID%SolutionFromFile)then

            call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)&
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR03'
        endif

        if (NewProperty%ID%IDNumber == SurfaceRadiation_) then

            call GetData(Me%ReflectionCoef,                                             &
                         Me%ObjEnterData, iflag,                                        &
                         keyword      = 'ALBEDO',                                       &  
                         default      = 0.0,                                            &
                         ClientModule = 'ModuleInterfaceWaterAir',                      &
                         SearchType   = FromBlock,                                      &
                         STAT         = STAT_CALL)
       
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR04'

        endif


        if (NewProperty%ID%IDNumber == WindStressX_) then

            !Shear Coefficient at the Atmosphere
            call GetData(Me%DefineCDWIND,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         keyword      = 'DEFINE_CDWIND',                                &  
                         default      = .false.,                                        &
                         ClientModule = 'ModuleInterfaceWaterAir',                      &
                         SearchType   = FromBlock,                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR05'

            !Shear Coefficient at the Atmosphere
            ! Mellor, Introduction to Physical Oceanography, p52 (1996)
            call GetData(Me%CDWIND,                                                     &
                         Me%ObjEnterData, iflag,                                        &
                         keyword      = 'CDWIND',                                       &  
                         ClientModule = 'ModuleInterfaceWaterAir',                      &
                         default      = 0.0015,                                         &
                         SearchType   = FromBlock,                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR06'

            allocate(NewProperty%FieldGrid(SizeILB:SizeIUB, SizeJLB:SizeJUB))

            NewProperty%FieldGrid(:,:) = FillValueReal

        elseif(NewProperty%ID%IDNumber == WindStressY_) then
                        
            allocate(NewProperty%FieldGrid(SizeILB:SizeIUB, SizeJLB:SizeJUB))

            NewProperty%FieldGrid(:,:) = FillValueReal

        endif


    end subroutine Construct_PropertyValues


   !--------------------------------------------------------------------------


    subroutine Construct_PropertyEvolution(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        real                                    :: ModelDT

        !Local-----------------------------------------------------------------
        integer                                 :: iflag
        real                                    :: ErrorAux, auxFactor, DTaux
        logical                                 :: VariableDT

        !----------------------------------------------------------------------


        !Time Step if the property field is variable in time
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR08'

            call GetVariableDT (Me%ObjTime, VariableDT, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR09'
    
            if (VariableDT) then
                
                NewProperty%Evolution%DTInterval       = ModelDT

            else

            !<BeginKeyword>
                !Keyword          : DT_INTERVAL
                !<BeginDescription>       
                   ! In the future a GET_DT_Hydro to a hydrodynamic module should be done to know the DT value  
                   ! that is been used to computethe the velocity field
                   ! By default the DTinterval is equal to the time step of the hydrodynamic model. In this case 
                   ! is admitted that the hydrodynamic model is computing the velocity field using a ADI method.
                !<EndDescription>
                !Type             : Real 
                !Default          : Time step of hydrodynamic model
                !File keyword     : DISPQUAL
                !Multiple Options : Do not have
                !Search Type      : FromBlock
                !Begin Block      : <beginproperty>
                !End Block        : <endproperty>
            !<EndKeyword>

            call GetData(NewProperty%Evolution%DTInterval,                          &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType     = FromBlock,                                &
                         keyword        = 'DT_INTERVAL',                            &
                         Default        = ModelDT,                                  &
                         ClientModule   = 'ModuleInterfaceWaterAir',                &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR10'


            if (NewProperty%evolution%DTInterval < (ModelDT)) then
                write(*,*) 
                write(*,*) ' Time step error.'
                stop       'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR11'

            elseif (NewProperty%evolution%DTInterval > (ModelDT)) then

                !Property DT  must be a multiple of the ModelDT
                auxFactor = NewProperty%evolution%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) ' Time step error.'
                    stop       'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR12'
                endif

                    ! Run period in seconds
                    DTaux = Me%EndTime - Me%ExternalVar%Now

                    !The run period   must be a multiple of the Property DT
                    auxFactor = DTaux / NewProperty%evolution%DTInterval

                    ErrorAux = auxFactor - int(auxFactor)
                    if (ErrorAux /= 0) then
                        write(*,*) 
                        write(*,*) ' Time step error.'
                        stop       'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR13'
                    endif

                endif

            endif

            NewProperty%Evolution%NextCompute = Me%ExternalVar%Now + NewProperty%Evolution%DTInterval
                                                
        else 

            call null_time(NewProperty%Evolution%NextCompute)

            NewProperty%evolution%DTInterval = FillValueReal

        endif

    end subroutine Construct_PropertyEvolution
    
    
    !-------------------------------------------------------------------------
    
    
    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: iflag

        !Begin----------------------------------------------------------------

        !<BeginKeyword>
            !Keyword          : OUTPUT_HDF
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends to write a HDF format file for this property
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%OutputHDF,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'OUTPUT_HDF',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceWaterAir',                        &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceWaterAir - ERR01'
           
        !<BeginKeyword>
            !Keyword          : TIME_SERIE
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends to write a time serie for this property
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%TimeSerie,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TIME_SERIE',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceWaterAir',                        &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceWaterAir - ERR02'


        !<BeginKeyword>
            !Keyword          : BOX_TIME_SERIE
            !<BeginDescription>       
                ! Checks out if the user pretends to write a time serie inside each box for this property
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : DISPQUAL
            !Multiple Options : 1 (.true.) , 0 (.false.) 
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%BoxTimeSerie,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'BOX_TIME_SERIE',                                 &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceWaterAir',                        &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceWaterAir - ERR03'

    end subroutine Construct_PropertyOutPut
    
    !--------------------------------------------------------------------------
    
    
    subroutine Add_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer     :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstProperty)) then
            Me%PropertiesNumber  = 1
            Me%FirstProperty     => NewProperty
            Me%LastProperty      => NewProperty
        else
            NewProperty%Prev     => Me%LastProperty
            Me%LastProperty%Next => NewProperty
            Me%LastProperty      => NewProperty
            Me%PropertiesNumber  = Me%PropertiesNumber + 1
        endif 

        !----------------------------------------------------------------------

    end subroutine Add_Property 

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: HDF5_CREATE
        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize2D%ILB 
        WorkIUB = Me%WorkSize2D%IUB 
        WorkJLB = Me%WorkSize2D%JLB 
        WorkJUB = Me%WorkSize2D%JUB 

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, trim(Me%Files%Results)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR02'
        
        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR03'

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",               &
                              Array2D = Me%ExtWater%Bathymetry,                     &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR04'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",            &
                              Array2D = Me%ExtWater%WaterPoints2D,                  &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR05'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR06'
       
        !----------------------------------------------------------------------

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty
        logical                                     :: OutputON
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        nullify(Me%OutPut%OutTime)

        OutputON = OFF

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))
            
            if(CurrentProperty%OutputHDF) OutputON = ON
            CurrentProperty => CurrentProperty%Next

        enddo

        if(OutputON)then

            call GetOutPutTime(Me%ObjEnterData,                              &
                               CurrentTime = Me%ExternalVar%Now,             &
                               EndTime     = Me%EndTime,                     &
                               keyword     = 'OUTPUT_TIME',                  &
                               SearchType  = FromFile,                       &
                               OutPutsTime = Me%OutPut%OutTime,              &
                               OutPutsOn   = Me%OutPut%Yes,                  &
                               STAT        = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'ConstructGlobalOutput - ModuleInterfaceWaterAir - ERR01' 

            if (Me%OutPut%Yes) then

                Me%OutPut%NextOutPut = 1

            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructGlobalOutput - ModuleInterfaceWaterAir - ERR02'
            endif 

        endif

    end subroutine ConstructGlobalOutput

    !------------------------------------------------------------------------
    

    subroutine Construct_Sub_Modules

        !Local-----------------------------------------------------------------
        type (T_Property),           pointer                 :: PropertyX            

        !----------------------------------------------------------------------

        Me%Coupled%TimeSerie%NumberOfProperties            = 0
        Me%Coupled%BoxTimeSerie%NumberOfProperties         = 0

        PropertyX => Me%FirstProperty

do1 :   do while (associated(PropertyX))

            if (PropertyX%TimeSerie) then
                Me%Coupled%TimeSerie%NumberOfProperties             = &
                Me%Coupled%TimeSerie%NumberOfProperties             + 1
                Me%Coupled%TimeSerie%Yes                            = ON
            endif

            if (PropertyX%BoxTimeSerie) then
                Me%Coupled%BoxTimeSerie%NumberOfProperties          = &
                Me%Coupled%BoxTimeSerie%NumberOfProperties          + 1
                Me%Coupled%BoxTimeSerie%Yes                         = ON
            endif

            PropertyX=>PropertyX%Next

        end do do1

        if(Me%Coupled%TimeSerie%Yes)then

            call Construct_Time_Serie
             
        endif
        
        if(Me%Coupled%BoxTimeSerie%Yes)then

            call StartOutputBoxFluxes
               
        endif

    end subroutine Construct_Sub_Modules


    !--------------------------------------------------------------------------


    subroutine Construct_Time_Serie

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
         
        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber  
        type(T_Property), pointer                           :: PropertyX
        integer                                             :: nProperties
        character(len=StringLength)                         :: TimeSerieLocationFile
        character(len=StringLength), dimension(:), pointer  :: PropertyList


        !----------------------------------------------------------------------

        !First checks out how many properties will have time series
        PropertyX   => Me%FirstProperty
        nProperties =  0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) nProperties = nProperties + 1
            PropertyX=>PropertyX%Next
        enddo

        if (nProperties > 0) then

            !Allocates PropertyList
            allocate(PropertyList(nProperties), STAT = STAT_CALL)
            if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR10'

            !Fills up PropertyList
            PropertyX   => Me%FirstProperty
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
                         Default      = Me%Files%InputData,                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR20' 


            !Constructs TimeSerie
            call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                            &
                                TimeSerieLocationFile,                                  &
                                PropertyList, "sri",                                    &
                                WaterPoints3D = Me%ExtWater%WaterPoints3D,              &
                                ModelName     = Me%ModelName,                           & 
                                STAT          = STAT_CALL)
            if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR30'

            !Deallocates PropertyList
            deallocate(PropertyList, STAT = STAT_CALL)
            if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR40'


            !Corrects if necessary the cell of the time serie based in the time serie coordinates
            call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR50'

            do dn = 1, TimeSerieNumber

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          CoordX   = CoordX,                                &
                                          CoordY   = CoordY,                                & 
                                          CoordON  = CoordON,                               &
                                          STAT     = STAT_CALL)
                if (CoordON) then
                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR60'

                    if (Id < 0 .or. Jd < 0) then
                
                        call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR70'

                        if (IgnoreOK) then
                            cycle
                        else
                            stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR80'
                        endif

                    endif


                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR90'
                endif

            enddo

        endif
        
    end subroutine Construct_Time_Serie

    !----------------------------------------------------------------------

    subroutine StartOutputBoxFluxes

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
        integer                                             :: ILB, IUB, JLB, JUB
        logical                                             :: Exist, Opened
 
        !Local-----------------------------------------------------------------
        type(T_Property),                           pointer :: PropertyX
        character(len=StringLength), dimension(:),  pointer :: ScalarOutputList
        integer                                             :: nScalars, n

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        !<BeginKeyword>
            !Keyword          : BOXFLUXES
            !<BeginDescription>       
            ! This keyword have two functions if exist fluxes between boxes are compute 
            ! and the value read is the name file where the boxes are defined
            !
            !<EndDescription>
            !Type             : Character 
            !Default          : Do not have
            !File keyword     : SEDPROP
            !Multiple Options : Do not have
            !Search Type      : From File
        !<EndKeyword>
        
        call GetData(Me%Files%BoxesFile,                                            &
                     Me%ObjEnterData, iflag,                                        &
                     keyword      = 'BOXFLUXES',                                    &
                     ClientModule = 'ModuleInterfaceWaterAir',                      &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR01'
        if (iflag .EQ. 0)                                                           &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR02'    
        
        inquire(File = Me%Files%BoxesFile, Exist = exist)
        if (exist) then
            inquire(File = Me%Files%BoxesFile, Opened  = Opened)
            if (opened) then
                write(*,*    ) 
                write(*,'(A)') 'BoxesFile = ',trim(adjustl(Me%Files%BoxesFile))
                write(*,*    ) 'Already opened.'
                stop           'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR03'    
            endif
        else
            write(*,*) 
            write(*,*)     'Could not find the boxes file.'
            write(*,'(A)') 'BoxFileName = ', Me%Files%BoxesFile
            stop           'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR04'    
        endif

        nScalars = Me%Coupled%BoxTimeSerie%NumberOfProperties 

        allocate(ScalarOutputList(nScalars), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR06'

        PropertyX  => Me%FirstProperty
        n = 0
        do while (associated(PropertyX))
            if (PropertyX%BoxTimeSerie) then
                n = n + 1
                ScalarOutputList(n) = trim(PropertyX%ID%name)
            endif 

            PropertyX=>PropertyX%Next
        end do

        call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                 &
                         TimeID             = Me%ObjTime,                   &
                         HorizontalGridID   = Me%ObjHorizontalGrid,         &
                         BoxesFilePath      = Me%Files%BoxesFile,           &
                         ScalarOutputList   = ScalarOutputList,             &
                         WaterPoints2D      = Me%ExtWater%WaterPoints2D,    &
                         STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR07'

        deallocate(ScalarOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR07'

        allocate(Me%Scalar2D(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR12'
        Me%Scalar2D(:,:) = 0.

    end subroutine StartOutputBoxFluxes

    
    !----------------------------------------------------------------------

    
    subroutine CheckInternalOptions

        !External--------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        type (T_Property), pointer                          :: PropertyY

        !----------------------------------------------------------------------
        
        call Search_Property(PropertyX, LatentHeat_,   STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%LatentHeat             = ON

        call Search_Property(PropertyX, SensibleHeat_, STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SensibleHeat           = ON

        call Search_Property(PropertyX, Evaporation_,  STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%Evaporation            = ON

        call Search_Property(PropertyX, NetLongWaveRadiation_,         STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%NetLongWaveRadiation     = ON

        call Search_Property(PropertyX, UpwardLongWaveRadiation_,      STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%UpwardLongWaveRadiation   = ON

        call Search_Property(PropertyX, DownwardLongWaveRadiation_,    STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%DownwardLongWaveRadiation = ON

        call Search_Property(PropertyX, NonSolarFlux_,           STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%NonSolarFlux           = ON

        call Search_Property(PropertyX, OxygenFlux_,             STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%OxygenFlux             = ON

        call Search_Property(PropertyX, CarbonDioxideFlux_,      STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%CarbonDioxideFlux      = ON
        
        call Search_Property(PropertyX, SpecificOxygenFlux_,     STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SpecificOxygenFlux     = ON

        call Search_Property(PropertyX, SpecificCarbonDioxideFlux_,  STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SpecificCarbonDioxideFlux  = ON
        
        call Search_Property(PropertyX, TurbulentKineticEnergy_, STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%TurbulentKineticEnergy = ON

        call Search_Property(PropertyX, SurfaceRadiation_,       STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SurfaceRadiation       = ON

        call Search_Property(PropertyX, SurfaceWaterFlux_,       STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SurfaceWaterFlux       = ON
        
        call Search_Property(PropertyX, WindStressX_,            STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            call Search_Property(PropertyY, WindStressY_,        STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_)  Me%IntOptions%WindStress        = ON
        endif


    end subroutine CheckInternalOptions


    !--------------------------------------------------------------------------
    !Module interface Water Air must provide to all other modules the interface fluxes.
    !If Interface fluxes are not specified in this module but are specified in another
    !module, model stops and writes a warning to the screen
    !--------------------------------------------------------------------------
    subroutine CheckOptionsWater

        !External--------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        type (T_Property), pointer                          :: PropertyY

        !----------------------------------------------------------------------
        
        !Check WaterProperties  
        call GetWaterPropertiesAirOptions(WaterPropertiesID    = Me%ObjWaterProperties,                & 
                                          TemperatureFluxYes   = Me%ExtOptions%HeatFluxYes,            &
                                          OxygenFluxYes        = Me%ExtOptions%OxygenFluxYes,          &
                                          CarbonDioxideFluxYes = Me%ExtOptions%CarbonDioxideFluxYes,   &
                                          WQMYes               = Me%ExtOptions%WQMYes,                 &
                                          T90VariableYes       = Me%ExtOptions%T90VariableYes,         &   
                                          STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR10'

        !Checks Hydrodynamic 
        call GetHydrodynamicAirOptions(HydrodynamicID      = Me%ObjHydrodynamic,                &
                                       SurfaceWaterFluxYes = Me%ExtOptions%SurfaceWaterFluxYes, &
                                       WindYes             = Me%ExtOptions%HydrodynamicWindYes, &
                                       AtmPressureYes      = Me%ExtOptions%HydrodynamicAtmPressureYes, &
                                       MslpYes             = Me%ExtOptions%HydrodynamicMslpYes, &
                                       STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR20'

#ifndef _LAGRANGIAN_
        if(Me%ObjLagrangian /= 0)then
#ifdef  _LAGRANGIAN_GLOBAL_                                         
        
            !Checks Lagrangian 
            call GetLagrangianAirOptionsGlobal(LagrangianID  = Me%ObjLagrangian,         &
                                         Oil           = Me%ExtOptions%OilYes,           &
                                         Wind          = Me%ExtOptions%LagrangianWindYes,&
                                         WaterQuality  = Me%ExtOptions%LagrangianWQMYes, &
                                         T90Variable   = Me%ExtOptions%LagrangianT90Yes, &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR30'
#else
            call GetLagrangianAirOptions(LagrangianID  = Me%ObjLagrangian,               &
                                         Oil           = Me%ExtOptions%OilYes,           &
                                         Wind          = Me%ExtOptions%LagrangianWindYes,&
                                         WaterQuality  = Me%ExtOptions%LagrangianWQMYes, &
                                         T90Variable   = Me%ExtOptions%LagrangianT90Yes, &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR35'
#endif         
        endif
#endif

#ifndef _WAVES_
        if(Me%ObjWaves    /= 0)Me%ExtOptions%WavesWindYes             = ON
#endif
        if(Me%ObjTurbGOTM /= 0)Me%ExtOptions%GOTMWindShearVelocityYes = ON

        
        if (Me%ExtOptions%HeatFluxYes) then

            call Search_Property(PropertyX, NonSolarFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR40'

            if(.not. PropertyX%ID%SolutionFromFile .and. .not. PropertyX%Constant) then

                call Search_Property(PropertyX, LatentHeat_, .true., STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR50'

                call Search_Property(PropertyX, SensibleHeat_, .true., STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR60'

                call Search_Property(PropertyX, NetLongWaveRadiation_, .true., STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR70'

                if(.not. PropertyX%ID%SolutionFromFile .and. .not. PropertyX%Constant) then

                    call Search_Property(PropertyX, UpwardLongWaveRadiation_, .true., STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR80'

                    call Search_Property(PropertyX, DownwardLongWaveRadiation_, .true., STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR90'
                endif

            endif

            call Search_Property(PropertyX, SurfaceRadiation_, .true., STAT = STAT_CALL)
            !if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR100'

            if(.not. PropertyX%ID%SolutionFromFile .and. .not. PropertyX%Constant) then

                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, SolarRadiation_))then
                    write(*,*) 'Missing SolarRadiation in Module Atmosphere '
                    stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR110'
                endif

            endif

        endif

        if (Me%ExtOptions%OxygenFluxYes) then

            call Search_Property(PropertyX, OxygenFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR120'
            
        endif

        if (Me%ExtOptions%CarbonDioxideFluxYes) then

            call Search_Property(PropertyX, CarbonDioxideFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR130'         
           
        endif




        if (Me%ExtOptions%SurfaceWaterFluxYes) then

            call Search_Property(PropertyX, SurfaceWaterFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR140'

            Me%ExtOptions%Precipitation = AtmospherePropertyExists(Me%ObjAtmosphere, Precipitation_)
            Me%ExtOptions%Irrigation    = AtmospherePropertyExists(Me%ObjAtmosphere, Irrigation_   )         
            

        endif



        if (Me%ExtOptions%HydrodynamicWindYes) then

            call Search_Property(PropertyX, WindStressX_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR150'


            call Search_Property(PropertyY, WindStressY_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR160'

            if ((      PropertyX%ID%SolutionFromFile .and. .not. PropertyY%ID%SolutionFromFile) .or. &
                (.not. PropertyX%ID%SolutionFromFile .and.       PropertyY%ID%SolutionFromFile)) then
                
                write (*,*) 'wind stress X must be given in the same way as wind stress Y'
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR0170'

            endif

        endif

        !Check if the Atmospheric Pressure property in the Atmosphere module exists
        if (Me%ExtOptions%HydrodynamicAtmPressureYes) then

            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AtmosphericPressure_)) then
                write(*,*) 'Missing AtmosphericPressure in Module Atmosphere'
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR180'
            endif

        endif

        !Check if the Mean Sea Level Pressure property in the Atmosphere module exists
        if (Me%ExtOptions%HydrodynamicMslpYes) then

            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, MeanSeaLevelPressure_)) then
                write(*,*) 'Missing mslp in Module Atmosphere'
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR181'
            endif

        endif

        if (Me%ExtOptions%OilYes) then

            !TO_DO Waves
            !Checks if one of the following Atmospheric properties exist:
            !Atmospheric Pressur or Mslp (Mean Sea Level Pressure)
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, AtmosphericPressure_)          &
                .and. .not. AtmospherePropertyExists(Me%ObjAtmosphere, MeanSeaLevelPressure_)) then

                write(*,*) 'Missing AtmosphericPressure or Mslp in Module Atmosphere '
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR190'

            endif

        endif

        if (Me%ExtOptions%LagrangianWindYes) then

            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_))then
                write(*,*) 'Missing WindVelocity in Module Atmosphere '
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR200'
            endif

            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_))then
                write(*,*) 'Missing WindVelocity in Module Atmosphere '
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR210'
            endif

        endif

        if (Me%ExtOptions%WQMYes           .or. Me%ExtOptions%T90VariableYes .or.       &
            Me%ExtOptions%LagrangianWQMYes .or. Me%ExtOptions%LagrangianT90Yes) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, SolarRadiation_))then
                write(*,*) 'Missing SolarRadiation in Module Atmosphere '
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR220'
            endif

        endif
        
        if (Me%ExtOptions%GOTMWindShearVelocityYes .and. Me%IntOptions%WindStress) then
            
            Me%IntOptions%WindShearVelocity = .true.
            
            allocate(Me%WindShearVelocity(Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
            
        else
        
            Me%IntOptions%WindShearVelocity = .false.
            
        endif


    end subroutine CheckOptionsWater
    
    !--------------------------------------------------------------------------
    !Module interface Water Air needs some atmosphere properties. This subroutine
    !checks if those properties are defined in atmosphere.
    !--------------------------------------------------------------------------
    
    subroutine CheckOptionsAir

        !External--------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB
        
        if (AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_))then
            
            Me%LocalAtm%WindVelocityX%Yes = ON
            nullify (Me%LocalAtm%WindVelocityX%Field)
            allocate(Me%LocalAtm%WindVelocityX%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%WindVelocityX%Field(:,:) = FillValueReal

        endif

        if (AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_))then
            
            Me%LocalAtm%WindVelocityY%Yes = ON
            nullify (Me%LocalAtm%WindVelocityY%Field)
            allocate(Me%LocalAtm%WindVelocityY%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%WindVelocityY%Field(:,:) = FillValueReal
        endif

        if (AtmospherePropertyExists(Me%ObjAtmosphere, WindDirection_))then
            
            Me%LocalAtm%WindDirection%Yes = ON
            nullify (Me%LocalAtm%WindDirection%Field)
            allocate(Me%LocalAtm%WindDirection%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%WindDirection%Field(:,:) = FillValueReal
        endif
        
        !Allocate field for the Atmospheric Pressure if present
        if (AtmospherePropertyExists(Me%ObjAtmosphere, AtmosphericPressure_)) then
           
            Me%LocalAtm%AtmosphericPressure%Yes = ON
            nullify (Me%LocalAtm%AtmosphericPressure%Field)
            allocate(Me%LocalAtm%AtmosphericPressure%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%AtmosphericPressure%Field(:,:) = FillValueReal
        endif

        !Allocate field for the Mean Sea Level Pressure if present
        if (AtmospherePropertyExists(Me%ObjAtmosphere, MeanSeaLevelPressure_)) then
           
            Me%LocalAtm%Mslp%Yes = ON
            nullify (Me%LocalAtm%Mslp%Field)
            allocate(Me%LocalAtm%Mslp%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%Mslp%Field(:,:) = FillValueReal
        endif

        call Search_Property(PropertyX, LatentHeat_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_)then
            if((.not. PropertyX%ID%SolutionFromFile) &
                .and. (.not. PropertyX%Constant)) then

                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_))then
                    write(*,*) 'Missing WindVelocity in Module Atmosphere '
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR10'
                endif

                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_))then
                    write(*,*) 'Missing WindVelocity in Module Atmosphere '
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR20'
                endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AirTemperature_)) then
                    write(*,*) 'Specify Atmosphere AirTemperature to calculate Interface'
                    write(*,*) 'LatentHeat'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR30'
                endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, RelativeHumidity_)) then
                    write(*,*) 'Specify Atmosphere RelativeHumidity to calculate Interface'
                    write(*,*) 'LatentHeat'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR40'
                endif
            endif

        endif

        call Search_Property(PropertyX, SensibleHeat_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_)then
            if ((.not. PropertyX%ID%SolutionFromFile) &
                 .and. (.not. PropertyX%Constant)) then

                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_))then
                    write(*,*) 'Missing WindVelocity X in Module Atmosphere '
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR50'
                endif

                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_))then
                    write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
                    stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR60'
                endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AirTemperature_)) then
                    write(*,*) 'Specify Atmosphere AirTemperature to calculate Interface'
                    write(*,*) 'SensibleHeat'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR70'
                endif
            endif
        endif

        call Search_Property(PropertyX, DownwardLongWaveRadiation_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if ((.not. PropertyX%ID%SolutionFromFile) &
                 .and. (.not. PropertyX%Constant)) then

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, RelativeHumidity_)) then
                    write(*,*) 'Specify Atmosphere RelativeHumidity to calculate Interface'
                    write(*,*) 'NetLongWaveRadiation'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR80'
                endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, CloudCover_)) then
                    write(*,*) 'Specify Atmosphere CloudCover to calculate Interface'
                    write(*,*) 'NetLongWaveRadiation'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR90'
                endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AirTemperature_)) then
                    write(*,*) 'Specify Atmosphere AirTemperature to calculate Interface'
                    write(*,*) 'SensibleHeat'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR100'
                endif
            endif

        endif


        call Search_Property(PropertyX, OxygenFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_))then
                write(*,*) 'Missing WindVelocity X in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR110'
            endif

            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_))then
                write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR120'
            endif

        endif
        
        
        call Search_Property(PropertyX, CarbonDioxideFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_))then
                write(*,*) 'Missing WindVelocity X in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR121'
            endif

            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_))then
                write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR122'
            endif

        endif


        call Search_Property(PropertyX, SpecificOxygenFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_))then
                write(*,*) 'Missing WindVelocity X in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR123'
            endif

            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_))then
                write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR124'
            endif

        endif
        
        
        call Search_Property(PropertyX, SpecificCarbonDioxideFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_))then
                write(*,*) 'Missing WindVelocity X in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR125'
            endif

            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_))then
                write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR126'
            endif

        endif


        call Search_Property(PropertyX, SurfaceRadiation_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, SolarRadiation_)) then
                write(*,*) 'Specify Atmosphere SolarRadiation to calculate Interface'
                write(*,*) 'SurfaceRadiation'
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR130'
            endif
        endif
        
        

    end subroutine CheckOptionsAir


    subroutine SetSubModulesConstructor 

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL 
        
        !Begin-----------------------------------------------------------------

        if(Me%ObjTurbGOTM /= 0)then

            if (.not. Me%Rugosity%ON)                                                   &
                stop 'SetSubModulesConstructor - ModuleInterfaceWaterAir - ERR10'

            call SetTurbGOTMSurfaceRugosity(TurbGOTMID      = Me%ObjTurbGOTM,           &
                                            SurfaceRugosity = Me%Rugosity%Field,        &
                                            STAT            = STAT_CALL)                
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'SetSubModulesConstructor - ModuleInterfaceWaterAir - ERR20'


        endif


    end subroutine SetSubModulesConstructor

    !----------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------
    subroutine Search_Property(PropertyX, PropertyXID, PrintWarning, STAT)


        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer             :: PropertyX
        integer                   , intent (IN)         :: PropertyXID
        logical,          optional, intent (IN)         :: PrintWarning
        integer,          optional, intent (OUT)        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX)) 
            if (PropertyX%ID%IDNumber == PropertyXID) then
                exit
            else
                PropertyX => PropertyX%Next                 
            endif
        end do

       !A PropertyX was found
       if (associated(PropertyX)) then
            STAT_ = SUCCESS_  
        else
            if (present(PrintWarning)) then
                if (PrintWarning) write (*,*)'Property Not Found in Module InterfaceWaterAir ',trim(GetPropertyName(PropertyXID))
            endif
            STAT_  = NOT_FOUND_ERR_  
        endif

        if (present(STAT))STAT = STAT_

    end subroutine Search_Property

    !--------------------------------------------------------------------------
    
    subroutine ReadLockExternalGlobal
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Now
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalGlobal - ModuleInterfaceWaterAir - ERR01'

        call GetGridCellArea (Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalGlobal - ModuleInterfaceWaterAir - ERR02'

    end subroutine ReadLockExternalGlobal

    
    !--------------------------------------------------------------------------

    
    subroutine ReadLockExternalWater
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !WaterPoints2D
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExtWater%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceWaterAir - ERR01'

        call GetDensity(WaterPropertiesID = Me%ObjWaterProperties,          &
                        Density           = Me%ExtWater%Density,            &
                        CurrentTime       = Me%ExternalVar%Now,             &
                        STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceWaterAir - ERR02'

        call GetGridData(Me%ObjGridData, Me%ExtWater%Bathymetry, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterfaceWaterAir - ERR04'

        !WaterPoints3D
        call GetWaterPoints3D(Me%ObjMap, Me%ExtWater%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterfaceWaterAir - ERR05'


    end subroutine ReadLockExternalWater


    !----------------------------------------------------------------------

    subroutine ReadUnlockExternalGlobal
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalGlobal - ModuleInterfaceWaterAir - ERR01'

        call null_time(Me%ExternalVar%Now)

    end subroutine ReadUnlockExternalGlobal

    
    !----------------------------------------------------------------------


    subroutine ReadUnlockExternalWater
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Ungets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap,  Me%ExtWater%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalGlobal - ModuleInterfaceWaterAir - ERR01'

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%Density, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceWaterAir - ERR02'

        call UnGetGridData(Me%ObjGridData, Me%ExtWater%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceWaterAir - ERR03'

        call UnGetMap(Me%ObjMap, Me%ExtWater%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceWaterAir - ERR04'

    end subroutine ReadUnlockExternalWater


    !----------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyInterfaceWaterAir(ObjInterfaceWaterAirID, LagrangianID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceWaterAirID, LagrangianID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleInterfaceWaterAir", "ModifyInterfaceWaterAir")

        STAT_ = UNKNOWN_

        call Ready(ObjInterfaceWaterAirID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if(Me%ObjLagrangian == 0 .and. LagrangianID /= 0)then
                Me%ObjLagrangian  = AssociateInstance(mLAGRANGIAN_, LagrangianID)
            endif

            call ReadLockExternalGlobal

            call ReadLockExternalWater

            call ModifyLocalAtmVariables

            call ModifyWaterAirFluxes

            call ModifyRugosity

            if(Me%Coupled%TimeSerie%Yes)            &
                call Output_TimeSeries

            if(Me%Coupled%BoxTimeSerie%Yes)         &
                call Output_BoxTimeSeries

            if(Me%OutPut%Yes)                       &
                call OutPut_Results_HDF

            call ReadUnlockExternalWater

            call ReadUnlockExternalGlobal

            call SetSubModulesModifier

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceWaterAir", "ModifyInterfaceWaterAir")

    end subroutine ModifyInterfaceWaterAir

    
    !--------------------------------------------------------------------------


    subroutine ModifyLocalAtmVariables

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-------------------------------------------------------------

        !Begin-------------------------------------------------------------
        
        if(Me%LocalAtm%WindVelocityX%Yes)then

            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
                                       Scalar       = Me%ExtAtm%WindVelocityX%Field,&
                                       ID           = WindVelocityX_,               &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR01'
            
            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%WindVelocityX%Field, Me%Size2D,         &
                                Me%ExtAtm%WindVelocityX%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%WindVelocityX%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR02'


        endif


        if(Me%LocalAtm%WindVelocityY%Yes)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
                                       Scalar       = Me%ExtAtm%WindVelocityY%Field,&
                                       ID           = WindVelocityY_,               &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR03'

            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%WindVelocityY%Field, Me%Size2D,         &
                                Me%ExtAtm%WindVelocityY%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%WindVelocityY%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR04'

        endif


        if(Me%LocalAtm%WindDirection%Yes)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
                                       Scalar       = Me%ExtAtm%WindDirection%Field,&
                                       ID           = WindDirection_,               &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR04a'

            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%WindDirection%Field, Me%Size2D,             &
                                Me%ExtAtm%WindDirection%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%WindDirection%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR04b'

        endif


        !Copy to local variable Atmospheric Pressure
        if(Me%LocalAtm%AtmosphericPressure%Yes)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
                                       Scalar       = Me%ExtAtm%AtmosphericPressure%Field,&
                                       ID           = AtmosphericPressure_,         &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR05'

            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%AtmosphericPressure%Field, Me%Size2D,   &
                                Me%ExtAtm%AtmosphericPressure%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AtmosphericPressure%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR06'

        endif

        !Copy to local variable Mean Sea Level Pressure
        if(Me%LocalAtm%Mslp%Yes)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
                                       Scalar       = Me%ExtAtm%Mslp%Field,&
                                       ID           = MeanSeaLevelPressure_,         &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR07'

            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%Mslp%Field, Me%Size2D,   &
                                Me%ExtAtm%Mslp%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%Mslp%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR08'

        endif

    end subroutine ModifyLocalAtmVariables

    !--------------------------------------------------------------------------

    subroutine ModifyWaterAirFluxes

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX, PropertyY

        !Begin-------------------------------------------------------------

        if(Me%IntOptions%LatentHeat)then
            call Search_Property(PropertyX, LatentHeat_,             STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyLatentHeat            (PropertyX)
        endif

        if(Me%IntOptions%SensibleHeat)then
            call Search_Property(PropertyX, SensibleHeat_,           STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifySensibleHeat          (PropertyX)
        endif

        if(Me%IntOptions%Evaporation)then
            call Search_Property(PropertyX, Evaporation_,            STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyEvaporation           (PropertyX)
        endif

        if(Me%IntOptions%NetLongWaveRadiation)then
            call Search_Property(PropertyX, NetLongWaveRadiation_,      STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyNetLongWaveRadiation     (PropertyX)
        endif

        if(Me%IntOptions%NonSolarFlux)then
            call Search_Property(PropertyX, NonSolarFlux_,           STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyNonSolarFlux          (PropertyX)
        endif

        if(Me%IntOptions%OxygenFlux)then
            call Search_Property(PropertyX, OxygenFlux_,             STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyOxygenFlux            (PropertyX)
        endif
           
        if(Me%IntOptions%CarbonDioxideFlux)then
            call Search_Property(PropertyX, CarbonDioxideFlux_,      STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyCarbonDioxideFlux  (PropertyX)
        endif

        if(Me%IntOptions%WindShearVelocity)then
            call ModifyWindShearVelocity     
        endif

        if(Me%IntOptions%TurbulentKineticEnergy)then
            call Search_Property(PropertyX, TurbulentKineticEnergy_, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyTurbulentKE           (PropertyX)
        endif


        if(Me%IntOptions%SurfaceRadiation)then
            call Search_Property(PropertyX, SurfaceRadiation_,       STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifySurfaceRadiation      (PropertyX)
        endif

        if(Me%IntOptions%SurfaceWaterFlux)then
            call Search_Property(PropertyX, SurfaceWaterFlux_,       STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifySurfaceWaterFlux      (PropertyX)
        endif
        
        if(Me%IntOptions%WindStress)then
            call Search_Property(PropertyX, WindStressX_,            STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                call Search_Property(PropertyY, WindStressY_,        STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_)  call ModifyWindStress    (PropertyX, PropertyY)
            endif
        endif

    end subroutine ModifyWaterAirFluxes

    !--------------------------------------------------------------------------

    subroutine ComputeLatentHeat (PropLatentHeat)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropLatentHeat

        !Local-----------------------------------------------------------------
        real,    dimension(:,:  ), pointer          :: UWIND, VWIND
        real,    dimension(:,:,:), pointer          :: WaterTemperature
        real,    dimension(:,:  ), pointer          :: AirTemperature, RelativeHumidity
        real                                        :: WindVelocity
        integer                                     :: ILB, IUB, JLB, JUB, KUB, i, j, STAT_CALL
        !Begin-----------------------------------------------------------------
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR01'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%RelativeHumidity%Field, &
                                   ID           = RelativeHumidity_,                &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR02'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AirTemperature%Field,   &
                                   ID           = AirTemperature_,                  &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR03'

            
        !Shorten Variables
        UWIND            => Me%LocalAtm%WindVelocityX%Field
        VWIND            => Me%LocalAtm%WindVelocityY%Field
        AirTemperature   => Me%ExtAtm%AirTemperature%Field
        WaterTemperature => Me%ExtWater%WaterTemperature
        RelativeHumidity => Me%ExtAtm%RelativeHumidity%Field

        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Calculates wind^2
                WindVelocity = sqrt(UWIND(i, j)**2. + VWIND(i, j)**2.)

                !We just consider the loss of heat (evaporation), once when water condenses
                !the heat gain is absorved by the air
                !Frank/Ramiro 26-10-2000

                !When the Water temperature > AirTemperature, the flux is negativ
!                PropLatentHeat%Field(i, j) = -1. * (19.0 + 0.95 * WindVelocity2)                  * &
!                                              (SaturatedVaporPressure(WaterTemperature(i, j,KUB)) - &
!                                              RelativeHumidity(i, j)                             * &
!                                              (SaturatedVaporPressure(AirTemperature  (i, j    ))))
                                              
                
!                if (PropLatentHeat%Field(i, j) > 0.) PropLatentHeat%Field(i, j) = 0.

                PropLatentHeat%Field(i, j) = LatentHeat    (ReferenceDensity, WaterTemperature(i, j, KUB),  &
                                                            AirTemperature  (i, j), RelativeHumidity(i, j), &
                                                            WindVelocity)
                   

            endif

        enddo
        enddo

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR04'

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%RelativeHumidity%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR05'

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AirTemperature%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR06'


        !Nullify auxiliar variables
        nullify(UWIND)
        nullify(VWIND)
        nullify(AirTemperature)
        nullify(WaterTemperature)
        nullify(RelativeHumidity)
    
    end subroutine ComputeLatentHeat


    subroutine ModifyLatentHeat(PropLatentHeat)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropLatentHeat

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropLatentHeat%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropLatentHeat%ID%ObjFillMatrix,  &
                                  Matrix2D          = PropLatentHeat%Field,             &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,        &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyLatentHeat - ModuleInterfaceWaterAir - ERR01'


        elseif(.not. PropLatentHeat%Constant)then

            call ComputeLatentHeat(PropLatentHeat)

        endif


    end subroutine ModifyLatentHeat

    
    !------------------------------------------------------------------------------

    
    subroutine ModifyNonSolarFlux(PropNonSolarFlux)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropNonSolarFLux

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropNonSolarFlux%ID%SolutionFromFile) then
            
            call ModifyFillMatrix(FillMatrixID      = PropNonSolarFlux%ID%ObjFillMatrix,    &
                                  Matrix2D          = PropNonSolarFlux%Field,               &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyNonSolarFlux - ModuleInterfaceWaterAir - ERR01'

        elseif (.not. PropNonSolarFlux%Constant) then

            call ComputeNonSolarFlux(PropNonSolarFlux)

        endif


    end subroutine ModifyNonSolarFlux
    
    
    !--------------------------------------------------------------------------

    subroutine ComputeNonSolarFlux(PropNonSolarFlux)

        !Arguments-------------------------------------------------------------
        type(T_Property),       pointer         :: PropNonSolarFlux
        
        !Local-----------------------------------------------------------------
        type(T_Property),       pointer         :: PropertyX
        real, dimension(:,:),   pointer         :: LatentHeat
        real, dimension(:,:),   pointer         :: SensibleHeat
        real, dimension(:,:),   pointer         :: NetLongWaveRadiation
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        call Search_Property(PropertyX, LatentHeat_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNonSolarFlux - ModuleInterfaceWaterAir - ERR01'
        LatentHeat => PropertyX%Field

        call Search_Property(PropertyX, SensibleHeat_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNonSolarFlux - ModuleInterfaceWaterAir - ERR02'
        SensibleHeat => PropertyX%Field


        call Search_Property(PropertyX, NetLongWaveRadiation_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNonSolarFlux - ModuleInterfaceWaterAir - ERR03'
        NetLongWaveRadiation => PropertyX%Field

        PropNonSolarFlux%Field = LatentHeat + SensibleHeat + NetLongWaveRadiation

        nullify(LatentHeat, SensibleHeat, NetLongWaveRadiation)

    end subroutine ComputeNonSolarFlux

    !--------------------------------------------------------------------------    

    
    subroutine ModifySurfaceRadiation (PropSurfaceRadiation)
        
        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropSurfaceRadiation

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropSurfaceRadiation%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropSurfaceRadiation%ID%ObjFillMatrix,  &
                                  Matrix2D          = PropSurfaceRadiation%Field,             &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,        &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifySurfaceRadiation - ModuleInterfaceWaterAir - ERR01a'


        elseif(.not. PropSurfaceRadiation%Constant)then

            call ComputeSurfaceRadiation(PropSurfaceRadiation)

        endif

    end subroutine ModifySurfaceRadiation


    !--------------------------------------------------------------------------    

    
    subroutine ComputeSurfaceRadiation (PropSurfaceRadiation)
        
        !Arguments--------------------------------------------------------------
        type(T_Property), pointer           :: PropSurfaceRadiation

        !External--------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, i, j
        real, dimension(:,:), pointer       :: SolarRadiation
        integer                             :: STAT_CALL
        integer                             :: CHUNK

        !Begin-----------------------------------------------------------------

        IUB = Me%WorkSize2D%IUB
        ILB = Me%WorkSize2D%ILB
        JUB = Me%WorkSize2D%JUB
        JLB = Me%WorkSize2D%JLB


        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,         &
                                   Scalar       = SolarRadiation,           &
                                   ID           = SolarRadiation_,          &
                                   STAT         = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_)                                        &
            stop 'ComputeSurfaceRadiation - ModuleInterfaceWaterAir - ERR02'

        if(Me%ReflectionCoef >= 0.)then

            if (MonitorPerformance) then
                call StartWatch ("ModuleInterfaceWaterAir", "ComputeSurfaceRadiation")
            endif

            CHUNK = CHUNK_J(JLB, JUB)
            !$OMP PARALLEL PRIVATE(i,j)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j=JLB, JUB
            do i=ILB, IUB
            
                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                    PropSurfaceRadiation%Field(i, j) = SolarRadiation(i, j) * (1 - Me%ReflectionCoef)
                                                                
                endif

            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            if (MonitorPerformance) then
                call StopWatch ("ModuleInterfaceWaterAir", "ComputeSurfaceRadiation")
            endif
            
        else

            call SetMatrixValue(PropSurfaceRadiation%Field, Me%Size2D, SolarRadiation)
        
        endif

        call UnGetAtmosphere(Me%ObjAtmosphere, SolarRadiation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeSurfaceRadiation - ModuleInterfaceWaterAir - ERR04'


    end subroutine ComputeSurfaceRadiation

    !--------------------------------------------------------------------------
        
    subroutine ComputeSensibleHeat(PropSensibleHeat)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropSensibleHeat

        !Local-----------------------------------------------------------------
        real,    dimension(:,:  ), pointer          :: UWIND, VWIND
        real,    dimension(:,:,:), pointer          :: WaterTemp
        real,    dimension(:,:  ), pointer          :: AirTemp
        real                                        :: WindVelocity
        integer                                     :: ILB, IUB, JLB, JUB, KUB
        integer                                     :: i, j, STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ComputeSensibleHeat - ModuleInterfaceWaterAir - ERR01'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AirTemperature%Field,   &
                                   ID           = AirTemperature_,                  &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeSensibleHeat - ModuleInterfaceWaterAir - ERR02'


        UWIND            => Me%LocalAtm%WindVelocityX%Field
        VWIND            => Me%LocalAtm%WindVelocityY%Field
        AirTemp          => Me%ExtAtm%AirTemperature%Field
        WaterTemp        => Me%ExtWater%WaterTemperature

        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Calculates the Modulus of the wind
                WindVelocity = sqrt(UWIND(i, j)**2. + VWIND(i, j)**2.)

                PropSensibleHeat%Field(i, j) = SensibleHeat (ReferenceDensity, WaterTemp(i, j, KUB), &
                                                             AirTemp(i, j), WindVelocity)

                !When the Water temperature > AirTemperature, the flux is negativ
!                PropSensibleHeat%Field(i, j) = -1. * BowenCoefficient               *   &
!                                                     (19.0 + 0.95 * WindVelocity2)  *   &
!                                                     (WaterTemp(i, j, KUB) - AirTemp(i, j))
            endif

        enddo
        enddo


        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeSensibleHeat - ModuleInterfaceWaterAir - ERR03'

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AirTemperature%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeSensibleHeat - ModuleInterfaceWaterAir - ERR04'


        !Nullify auxiliar variables
        nullify(UWIND)
        nullify(VWIND)
        nullify(AirTemp)
        nullify(WaterTemp)

    end subroutine ComputeSensibleHeat

    
    !--------------------------------------------------------------------------
    
    subroutine ModifySensibleHeat(PropSensibleHeat)

        !Arguments-------------------------------------------------------------

        type(T_Property), pointer :: PropSensibleHeat

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL


        !Begin-----------------------------------------------------------------

        if (PropSensibleHeat%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropSensibleHeat%ID%ObjFillMatrix,    &
                                  Matrix2D          = PropSensibleHeat%Field,               &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifySensibleHeat - ModuleInterfaceWaterAir - ERR01'

        elseif (.not. PropSensibleHeat%Constant) then

            call ComputeSensibleHeat(PropSensibleHeat)

        endif

    end subroutine ModifySensibleHeat
    
    !------------------------------------------------------------------------------

    subroutine ComputeEvaporation(PropEvaporation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropEvaporation

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: PropLatentHeat
        integer                                     :: ILB, IUB, JLB, JUB, i, j
        integer                                     :: STAT_CALL
        integer                                     :: CHUNK
         
        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB

        !Looks if the property also calculates latent heat. If so, use the evaporation
        !corresponding to the latent heat, if not check for other sources.
        call Search_Property (PropLatentHeat, LatentHeat_, STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then

            if (MonitorPerformance) then
                call StartWatch ("ModuleInterfaceWaterAir", "ComputeEvaporation")
            endif

            CHUNK = CHUNK_J(JLB, JUB)
            !$OMP PARALLEL PRIVATE(i,j)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB

                ![m3/s]                    = [J/m2/s] / [J/kg] / [kg/m3] * [m] * [m]
                PropEvaporation%Field(i, j) = PropLatentHeat%Field(i, j)         / &
                                              LatentHeatOfVaporization           / &
                                              ReferenceDensity                   * &
                                              Me%ExternalVar%GridCellArea(i, j)
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
        
            if (MonitorPerformance) then
                call StopWatch ("ModuleInterfaceWaterAir", "ComputeEvaporation")
            endif
            
        endif


    end subroutine ComputeEvaporation

    !--------------------------------------------------------------------------
    
    subroutine ModifyEvaporation(PropEvaporation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropEvaporation

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropEvaporation%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropEvaporation%ID%ObjFillMatrix,     &
                                  Matrix2D          = PropEvaporation%Field,                &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR01'

        elseif(.not. PropEvaporation%Constant)then

            call ComputeEvaporation(PropEvaporation)

        endif

    
    end subroutine ModifyEvaporation

    
    !--------------------------------------------------------------------------

    
    subroutine ModifySurfaceWaterFlux(PropSurfaceWaterFlux)

        !Arguments-------------------------------------------------------------
        type(T_Property),       pointer             :: PropSurfaceWaterFlux

        !Local-----------------------------------------------------------------
        type(T_Property),       pointer             :: PropEvaporation
        real, dimension(:,:),   pointer             :: Precipitation
        real, dimension(:,:),   pointer             :: Irrigation
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call SetMatrixValue(PropSurfaceWaterFlux%Field, Me%Size2D, 0.)

        if(Me%IntOptions%Evaporation)then
            !Gets Evaporation
            call Search_Property(PropEvaporation, Evaporation_, .false., STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_ .and. STAT_CALL .ne. NOT_FOUND_ERR_)then
                stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR01'
            elseif(STAT_CALL .eq. SUCCESS_)then
                PropSurfaceWaterFlux%Field = PropSurfaceWaterFlux%Field + PropEvaporation%Field
            endif
        endif


        if(Me%ExtOptions%Precipitation)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,         &
                                       Scalar       = Precipitation,            &
                                       ID           = Precipitation_,           &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                        &
                stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR02'
        
            PropSurfaceWaterFlux%Field = PropSurfaceWaterFlux%Field + Precipitation

            call UnGetAtmosphere(Me%ObjAtmosphere, Precipitation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR04'


        endif

        if(Me%ExtOptions%Irrigation)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,         &
                                       Scalar       = Irrigation,               &
                                       ID           = Irrigation_,              &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                        &
                stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR03'
            
        
            PropSurfaceWaterFlux%Field = PropSurfaceWaterFlux%Field + Irrigation
            
            call UnGetAtmosphere(Me%ObjAtmosphere, Irrigation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR05'

        endif

    end subroutine ModifySurfaceWaterFlux
    
    !--------------------------------------------------------------------------

    subroutine ComputeUpLongWaveRad(PropUpwardLongWaveRadiation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropUpwardLongWaveRadiation

        !Local-----------------------------------------------------------------
        real, dimension   (:,:,:), pointer          :: WaterTemp
        integer                                     :: ILB, IUB, JLB, JUB, KUB, i, j
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ComputeUpLongWaveRad - ModuleInterfaceWaterAir - ERR10'

       
        WaterTemp        => Me%ExtWater%WaterTemperature


        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                PropUpwardLongWaveRadiation%Field(i, j) = LongWaveUpward (WaterTemp(i, j, KUB))


            endif
        enddo
        enddo

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeUpLongWaveRad - ModuleInterfaceWaterAir - ERR20'

         nullify(WaterTemp)

    end subroutine ComputeUpLongWaveRad

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine ComputeDownLongWaveRad(PropDownwardLongWaveRadiation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropDownwardLongWaveRadiation

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, KUB, i, j
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB

       
        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%RelativeHumidity%Field, &
                                   ID           = RelativeHumidity_,                &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR10'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%CloudCover%Field,       &
                                   ID           = CloudCover_,                      &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR20'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AirTemperature%Field,   &
                                   ID           = AirTemperature_,                  &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR30'


        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                PropDownwardLongWaveRadiation%Field(i, j) = LongWaveDownward (Me%ExtAtm%CloudCover%Field    (i, j),      &
                                                                              Me%ExtAtm%AirTemperature%Field(i, j)) 

            endif
        enddo
        enddo

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%RelativeHumidity%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR40'
        
        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%CloudCover%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR50'

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AirTemperature%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR60'


    end subroutine ComputeDownLongWaveRad

    !--------------------------------------------------------------------------
    
    subroutine ModifyNetLongWaveRadiation(PropNetLongWaveRadiation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropNetLongWaveRadiation
        
        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: PropUpLongWaveRadiation
        type(T_Property), pointer                   :: PropDownLongWaveRadiation
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

i1:     if (PropNetLongWaveRadiation%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropNetLongWaveRadiation%ID%ObjFillMatrix,   &
                                  Matrix2D          = PropNetLongWaveRadiation%Field,              &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,                   &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyNetLongWaveRadiation - ModuleInterfaceWaterAir - ERR10'


        elseif (.not. PropNetLongWaveRadiation%Constant) then i1

            call Search_Property(PropUpLongWaveRadiation,       UpwardLongWaveRadiation_,      STAT = STAT_CALL)
            call Search_Property(PropDownLongWaveRadiation,     DownwardLongWaveRadiation_,    STAT = STAT_CALL)
            
i2:         if (PropUpLongWaveRadiation%ID%SolutionFromFile) then

                call ModifyFillMatrix(FillMatrixID      = PropUpLongWaveRadiation%ID%ObjFillMatrix,&
                                      Matrix2D          = PropUpLongWaveRadiation%Field,           &
                                      PointsToFill2D    = Me%ExtWater%WaterPoints2D,               &
                                      STAT              = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop 'ModifyNetLongWaveRadiation - ModuleInterfaceWaterAir - ERR20'

            elseif (.not. PropUpLongWaveRadiation%Constant) then i2

                call ComputeUpLongWaveRad(PropUpLongWaveRadiation)

            endif i2


i3:         if (PropDownLongWaveRadiation%ID%SolutionFromFile) then

                call ModifyFillMatrix(FillMatrixID      = PropDownLongWaveRadiation%ID%ObjFillMatrix,&
                                      Matrix2D          = PropDownLongWaveRadiation%Field,           &
                                      PointsToFill2D    = Me%ExtWater%WaterPoints2D,                     &
                                      STAT              = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop 'ModifyNetLongWaveRadiation - ModuleInterfaceWaterAir - ERR30'

            elseif (.not. PropDownLongWaveRadiation%Constant) then i3

                call ComputeDownLongWaveRad(PropDownLongWaveRadiation)

            endif i3

            PropNetLongWaveRadiation%Field(:,:) = PropUpLongWaveRadiation%Field(:,:) + &
                                                  PropDownLongWaveRadiation%Field(:,:)
        endif i1

    end subroutine ModifyNetLongWaveRadiation
    
    !--------------------------------------------------------------------------

    subroutine ModifyWindStress(PropWindStressX, PropWindStressY)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindStressX, PropWindStressY

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleInterfaceWaterAir", "ModifyWindStress")

        if (PropWindStressX%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropWindStressX%ID%ObjFillMatrix,     &
                                  Matrix2D          = PropWindStressX%Field,                &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyWindStress - ModuleInterfaceWaterAir - ERR01'

        endif

        if (PropWindStressY%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropWindStressY%ID%ObjFillMatrix,     &
                                  Matrix2D          = PropWindStressY%Field,                &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyWindStress - ModuleInterfaceWaterAir - ERR20'

        endif


        call RotateVectorFieldToGrid(HorizontalGridID  = Me%ObjHorizontalGrid,      &
                                     VectorInX         = PropWindStressX%Field,     &
                                     VectorInY         = PropWindStressY%Field,     &
                                     VectorOutX        = PropWindStressX%FieldGrid, &
                                     VectorOutY        = PropWindStressY%FieldGrid, &   
                                     WaterPoints2D     = Me%ExtWater%WaterPoints2D, &
                                     RotateX           = .true.,                    &
                                     RotateY           = .true.,                    &
                                     STAT              = STAT_CALL)

        if(STAT_CALL .ne. SUCCESS_) stop 'ModifyWindStress - ModuleInterfaceWaterAir - ERR30'


        if (.not. PropWindStressX%Constant            .and. &
            .not. PropWindStressY%Constant            .and. &
            .not. PropWindStressX%ID%SolutionFromFile .and. &
            .not. PropWindStressY%ID%SolutionFromFile ) then

            call ComputeTauWind (PropWindStressX, PropWindStressY)

            call RotateVectorGridToField(HorizontalGridID  = Me%ObjHorizontalGrid,      &
                                         VectorInX         = PropWindStressX%FieldGrid, &
                                         VectorInY         = PropWindStressY%FieldGrid, &
                                         VectorOutX        = PropWindStressX%Field    , &
                                         VectorOutY        = PropWindStressY%Field    , &   
                                         WaterPoints2D     = Me%ExtWater%WaterPoints2D, &
                                         RotateX           = .true.,                    &
                                         RotateY           = .true.,                    &
                                         STAT              = STAT_CALL)

        endif

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceWaterAir", "ModifyWindStress")

    end subroutine ModifyWindStress

    !--------------------------------------------------------------------------
    
    subroutine ComputeTauWind (PropWindStressX, PropWindStressY)

        !Arguments------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindStressX
        type(T_Property), pointer                   :: PropWindStressY

        !Local----------------------------------------------------------------
        real,    dimension(:,:), pointer            :: UWIND, VWIND
        integer                                     :: IUB, ILB, JUB, JLB, i, j
        real                                        :: VM, CDWIND
        real                                        :: Coef
        integer                                     :: CHUNK

        !Begin----------------------------------------------------------------

        IUB = Me%WorkSize2D%IUB
        ILB = Me%WorkSize2D%ILB
        JUB = Me%WorkSize2D%JUB
        JLB = Me%WorkSize2D%JLB

        UWIND   => Me%LocalAtm%WindVelocityX%Field
        VWIND   => Me%LocalAtm%WindVelocityY%Field

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceWaterAir", "ComputeTauWind")
        endif

        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,VM,Coef,CDWIND)
cd1:    if(Me%DefineCDWIND)then
        
            !The next line is done by each thread with OpenMP compilation.
            !This is ok because it is a private variable.
            CDWIND = Me%CDWIND
        
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j=JLB, JUB
            do i=ILB, IUB

                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                    
                    !Compute the square root 
                    VM = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)
                    
                    Coef = CDWIND * VM * Air_Density
                    
                    !Mellor, Introduction to Physical Oceanography, p52 (1996)---------
                    PropWindStressX%FieldGrid(I,J) = Coef * UWIND(I,J)
                    
                    !Mellor, Introduction to Physical Oceanography, p52 (1996)---------
                    PropWindStressY%FieldGrid(I,J) = Coef * VWIND(I,J)

                endif
                
            enddo
            enddo
            !$OMP END DO

        else cd1
        
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j=JLB, JUB
            do i=ILB, IUB

                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                    
                    !Compute the velocity modulus
                    VM = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)  
                    
                    !Compute the wind drag coefficient based on Large & Pond, 1981 
                    if (VM < 4.)then

                        !The Large & Pond, 1981 formulation is not valid for wind lower than 4 m/s
                        !A constant value is assumed   
                        CDWIND = 0.0012 
                    
                    elseif (VM >=  4  .and. VM <  11.)then

                        CDWIND = 0.0012

                    elseif (VM >= 11. .and. VM <= 26.)then

                        CDWIND = (0.00049 + 0.000065 * VM)                   

                    elseif (VM > 26.) then
                        
                        !The Large & Pond, 1981 formulation is not valid for wind higher than 26 m/s
                        !A constant value is assumed   
                        CDWIND = 0.00218

                    endif
                    
                    Coef = CDWIND * VM * Air_Density

                    !Mellor, Introduction to Physical Oceanography, p52 (1996)---------
                    PropWindStressX%FieldGrid(I,J) = Coef * UWIND(I,J)

                    !Mellor, Introduction to Physical Oceanography, p52 (1996)---------
                    PropWindStressY%FieldGrid(I,J) = Coef * VWIND(I,J)

                endif

            enddo
            enddo
            !$OMP END DO
            
        endif cd1
        !$OMP END PARALLEL

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceWaterAir", "ComputeTauWind")
        endif

        nullify(UWIND)
        nullify(VWIND)

    end subroutine ComputeTauWind        

    !--------------------------------------------------------------------------

    subroutine ModifyTurbulentKE(PropTurbulentKE)

        !Arguments-------------------------------------------------------------
        
        type(T_Property), pointer :: PropTurbulentKE

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropTurbulentKE%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropTurbulentKE%ID%ObjFillMatrix,   &
                                  Matrix2D          = PropTurbulentKE%Field,              &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,          &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyTurbulentKE - ModuleInterfaceWaterAir - ERR02'

        else if (.not. PropTurbulentKE%Constant) then

            call ComputeTKEWind (PropTurbulentKE) 

        endif

    end subroutine ModifyTurbulentKE

    !--------------------------------------------------------------------------

    subroutine ComputeTKEWind (PropTurbulentKE)

        !Arguments------------------------------------------------------------
        type(T_Property),       pointer             :: PropTurbulentKE

        !Local----------------------------------------------------------------
        real, dimension(:,:),   pointer             :: UWIND, VWIND
        integer                                     :: IUB, ILB, JUB, JLB, i, j
        real                                        :: VM
        integer                                     :: CHUNK

        !Begin----------------------------------------------------------------

        IUB = Me%WorkSize2D%IUB
        ILB = Me%WorkSize2D%ILB
        JUB = Me%WorkSize2D%JUB
        JLB = Me%WorkSize2D%JLB

        UWIND   => Me%LocalAtm%WindVelocityX%Field
        VWIND   => Me%LocalAtm%WindVelocityY%Field

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceWaterAir", "ComputeTKEWind")
        endif
    
        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,VM)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do3:    do j=JLB, JUB
do4:    do i=ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Compute the square root 
                VM     = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)       

                ![m/s2 * kg * m / m2 / s]         = [kg/m^3 * (m/s)^3]     
                PropTurbulentKE%Field(I,J) = (Air_Density * CDE * VM ** 3.0) 
            endif

        enddo do4
        enddo do3
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceWaterAir", "ComputeTKEWind")
        endif

        nullify(UWIND)
        nullify(VWIND)

    end subroutine ComputeTKEWind 
    
  !-------------------------------------------------------------------------- 

    subroutine ModifyOxygenFlux(PropertyO2)

        !Arguments--------------------------------------------------------------
        type (T_Property),        pointer       :: PropertyO2

        !Local------------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Begin------------------------------------------------------------------

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ModifyOxygenFlux - ModuleInterfaceWaterAir - ERR01'

        call ModifyAerationFlux(PropertyO2%Field)

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyOxygenFlux - ModuleInterfaceWaterAir - ERR02'

    end subroutine ModifyOxygenFlux

    !-------------------------------------------------------------------------- 

    subroutine ModifyAerationFlux(Flux)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: Flux

        !Local------------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: UWIND, VWIND
        real,    dimension(:,:,:),    pointer   :: WaterTemperature
        integer                                 :: i,j
        real                                    :: WindVelocity, WaterTemp
        integer                                 :: ILB, IUB, JLB, JUB, KUB
        
        !Begin------------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB 
        IUB = Me%WorkSize2D%IUB 
        JLB = Me%WorkSize2D%JLB 
        JUB = Me%WorkSize2D%JUB 
        KUB = Me%WorkSize3D%KUB

        !Surface temperature of the waterbody
        WaterTemperature => Me%ExtWater%WaterTemperature
        UWIND            => Me%LocalAtm%WindVelocityX%Field
        VWIND            => Me%LocalAtm%WindVelocityY%Field

        !The formulation is taken from CE-QUAL-W2 v3.1
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Module of the wind velocity 
                WindVelocity = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)
                
                WaterTemp = WaterTemperature(i, j, KUB)
                
                Flux(i,j) = AerationFlux(Me%AerationEquation, WindVelocity, WaterTemp)
               
            endif

        end do
        end do
         

        !Nullify auxiliar variables
        nullify(UWIND)
        nullify(VWIND)
        nullify(WaterTemperature)
        
    end subroutine ModifyAerationFlux
    

    !---------------------------------------------------------------------------
    
    subroutine ModifyCarbonDioxideFlux(PropCarbonDioxide)

        !Arguments--------------------------------------------------------------
        type (T_Property),         pointer :: PropCarbonDioxide

        !External---------------------------------------------------------------
        integer                            :: STAT_CALL

        !Local------------------------------------------------------------------
        integer                            :: ILB, IUB, JLB, JUB, KUB
        real,    dimension(:,:,:), pointer :: WaterTemperature
        integer                            :: i, j
        integer                            :: CHUNK
        
        !Begin------------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB 
        IUB = Me%WorkSize2D%IUB 
        JLB = Me%WorkSize2D%JLB 
        JUB = Me%WorkSize2D%JUB 
        KUB = Me%WorkSize3D%KUB

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ModifyCarbonDioxideFlux - ModuleInterfaceWaterAir - ERR10'


        If (Me%Coupled%CEQUALW2%Yes) then
        
            !Surface temperature of the waterbody
            WaterTemperature   => Me%ExtWater%WaterTemperature

            call ModifyCO2AerationFlux(PropCarbonDioxide%Field)

            if (MonitorPerformance) then
                call StartWatch ("ModuleInterfaceWaterAir", "ModifyCarbonDioxideFlux")
            endif

            CHUNK = CHUNK_J(JLB, JUB)
            !$OMP PARALLEL PRIVATE(i,j)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                    PropCarbonDioxide%Field(i,j)= PropCarbonDioxide%Field(i,j) * 0.923              * &
                                                  (0.286 * exp(-0.0314*(WaterTemperature(i,j,KUB))) * &
                                                   Me%AltitudeCorrection)

                endif

            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            if (MonitorPerformance) then
                call StopWatch ("ModuleInterfaceWaterAir", "ModifyCarbonDioxideFlux")
            endif         

            call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModuleInterfaceWaterAir - ERR040'


            !Nullify auxiliar variables
            nullify(WaterTemperature)
        
        else
        
            call ModifyCO2AerationFlux(PropCarbonDioxide%Field)

            call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModuleInterfaceWaterAir - ERR050'
 
        
        endif
        

    end subroutine ModifyCarbonDioxideFlux

   !-------------------------------------------------------------------------- 
    
    
    subroutine ModifyCO2AerationFlux(Flux)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: Flux

        !Local------------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: UWIND, VWIND
        real,    dimension(:,:,:),    pointer   :: WaterTemperature, WaterSalinity 
        real,    dimension(:,:,:),    pointer   :: Velocity_U, Velocity_V
        integer                                 :: i,j
        real                                    :: WindVelocity, WaterVel, WaterTemp, WaterSal
        integer                                 :: ILB, IUB, JLB, JUB, KUB, STAT_CALL
        
        !Begin------------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB 
        IUB = Me%WorkSize2D%IUB 
        JLB = Me%WorkSize2D%JLB 
        JUB = Me%WorkSize2D%JUB 
        KUB = Me%WorkSize3D%KUB

        !Surface temperature of the waterbody
        WaterTemperature => Me%ExtWater%WaterTemperature
        UWIND            => Me%LocalAtm%WindVelocityX%Field
        VWIND            => Me%LocalAtm%WindVelocityY%Field
        
        
        call GetHorizontalVelocity(HydrodynamicID      = Me%ObjHydrodynamic,        &
                                   Velocity_U          = Velocity_U,                &
                                   Velocity_V          = Velocity_V,                &
                                   STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR10'
        
        
        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterSalinity,        &
                              PropertyXIDNumber = Salinity_,                        &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR15'

        WaterSalinity    => Me%ExtWater%WaterSalinity

        
        !The formulation is taken from CE-QUAL-W2 v3.1
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Module of the wind velocity 
                WindVelocity = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)
                
                WaterTemp = WaterTemperature(i, j, KUB)
                
                WaterSal = WaterSalinity(i, j, KUB)
                
                WaterVel  = SQRT(Velocity_U(i, j, KUB)**2. + Velocity_V(i, j, KUB)**2.)
                
                Flux(i,j) = AerationFlux_CO2(Me%CO2AerationEquation, WindVelocity, WaterVel, WaterTemp, WaterSal)
               
            endif

        end do
        end do
 
             
        call UngetHydrodynamic(Me%ObjHydrodynamic, Velocity_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR020'    
        
        call UngetHydrodynamic(Me%ObjHydrodynamic, Velocity_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR030'     

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterSalinity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR040'

        !Nullify auxiliar variables
        nullify(UWIND)
        nullify(VWIND)
        nullify(Velocity_U)
        nullify(Velocity_V)
        nullify(WaterTemperature)
        nullify(WaterSalinity)

        
    end subroutine ModifyCO2AerationFlux
    
    !--------------------------------------------------------------------------  

    subroutine ModifyWindShearVelocity

        !Arguments--------------------------------------------------------------


        !External--------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, i, j
        real                                :: WindStressModule
        type(T_Property), pointer           :: WindStressX, WindStressY
        integer                             :: STAT_CALL
        integer                             :: CHUNK

        !Begin-----------------------------------------------------------------

        IUB = Me%WorkSize2D%IUB
        ILB = Me%WorkSize2D%ILB
        JUB = Me%WorkSize2D%JUB
        JLB = Me%WorkSize2D%JLB
        
        call Search_Property(WindStressX, WindStressX_, .true., STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ModifyWindShearVelocity - ModuleInterfaceWaterAir - ERR01'


        call Search_Property(WindStressY, WindStressY_, .true., STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ModifyWindShearVelocity - ModuleInterfaceWaterAir - ERR02'

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceWaterAir", "ModifyWindShearVelocity")
        endif

        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,WindStressModule)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do1:    do j=JLB, JUB
do2:    do i=ILB, IUB
            
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                WindStressModule                  = sqrt(WindStressX%Field(i, j)**2. + &
                                                         WindStressY%Field(i, j)**2.)

               Me%WindShearVelocity(i, j) = sqrt(WindStressModule / SigmaDensityReference)
                                                                
            endif

        enddo do2
        enddo do1
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceWaterAir", "ModifyWindShearVelocity")
        endif

    end subroutine ModifyWindShearVelocity

    !------------------------------------------------------------------------------

    subroutine ComputeWavesRugosity

        !Arguments-------------------------------------------------------------
#ifndef _WAVES_
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer               :: WaveHeight
        integer                                     :: ILB, IUB, JLB, JUB, i, j
        integer                                     :: STAT_CALL
        integer                                     :: CHUNK

        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB

        !Gets wave height
        call GetWaves (Me%ObjWaves, WaveHeight = WaveHeight, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ComputeWavesRugosity - InterfaceWaterAir - ERR10'

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceWaterAir", "ComputeWavesRugosity")
        endif

        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                Me%Rugosity%Field(i, j) = Me%Rugosity%WavesRelation * WaveHeight(i, j)
            endif
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceWaterAir", "ComputeWavesRugosity")
        endif
        
        call UnGetWaves (Me%ObjWaves, WaveHeight, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeWavesRugosity - InterfaceWaterAir - ERR20'
    
#endif
    end subroutine ComputeWavesRugosity


    !--------------------------------------------------------------------------
    
    subroutine ModifyRugosity

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (Me%Rugosity%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = Me%Rugosity%ID%ObjFillMatrix,     &
                                  Matrix2D          = Me%Rugosity%Field,                &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyRugosity - ModuleInterfaceWaterAir - ERR01'

        elseif(.not. Me%Rugosity%Constant)then

            if (Me%Rugosity%WavesFunction) then

                call ComputeWavesRugosity

            endif

        endif
    
    end subroutine ModifyRugosity

    !--------------------------------------------------------------------------
    
    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX

        !Begin-----------------------------------------------------------------
        
        PropertyX   => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                call WriteTimeSerie(Me%ObjTimeSerie,                            &
                                    Data2D  = PropertyX%Field,                  &
                                    STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'OutPut_TimeSeries - ModuleInterfaceWaterAir - ERR01'
            endif
            PropertyX=>PropertyX%Next
        enddo
    
    end subroutine OutPut_TimeSeries


    !--------------------------------------------------------------------------

    subroutine OutPut_BoxTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: CHUNK

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB 
        IUB = Me%Size2D%IUB 
        JLB = Me%Size2D%JLB 
        JUB = Me%Size2D%JUB 

        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
            if (PropertyX%BoxTimeSerie) then

                Me%Scalar2D(:,:) = 0.

                if (MonitorPerformance) then
                    call StartWatch ("ModuleInterfaceWaterAir", "OutPut_BoxTimeSeries")
                endif

                CHUNK = CHUNK_J(JLB, JUB)
                !$OMP PARALLEL PRIVATE(I,J)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do J = JLB, JUB
                do I = ILB, IUB

                    if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                        Me%Scalar2D(i,j) = PropertyX%Field(i, j) * Me%ExternalVar%GridCellArea(i, j)
                    
                    endif

                end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL
                
                if (MonitorPerformance) then
                    call StopWatch ("ModuleInterfaceWaterAir", "OutPut_BoxTimeSeries")
                endif
                               
                call BoxDif(Me%ObjBoxDif,                           &
                            Me%Scalar2D,                            &
                            "WaterAir "//trim(PropertyX%ID%name),   &
                            Me%ExtWater%WaterPoints2D,              &
                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  &
                    stop 'OutPut_BoxTimeSeries - ModuleInterfaceWaterAir - ERR01'

                Me%Scalar2D(:,:) = null_real

            endif
            PropertyX=>PropertyX%Next
        enddo

    end subroutine OutPut_BoxTimeSeries

    
    !--------------------------------------------------------------------------


    subroutine OutPut_Results_HDF
        
        !External--------------------------------------------------------------
        integer                            :: STAT_CALL
        real                               :: Year, Month, Day, Hour, Minute, Second
         
        !Local-----------------------------------------------------------------
        type (T_Property), pointer         :: PropertyX
        integer                            :: OutPutNumber
        integer, dimension(6    )          :: TimeAux
        real,    dimension(6    ), target  :: AuxTime
        real,    dimension(:    ), pointer :: TimePtr
        integer                            :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                            :: WorkKLB, WorkKUB

        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize3D%ILB 
        WorkIUB = Me%WorkSize3D%IUB 
        WorkJLB = Me%WorkSize3D%JLB 
        WorkJUB = Me%WorkSize3D%JUB
        WorkKLB = Me%WorkSize3D%KLB 
        WorkKUB = Me%WorkSize3D%KUB


        OutPutNumber = Me%OutPut%NextOutPut


TOut:   if (Me%ExternalVar%Now >= Me%OutPut%OutTime(OutPutNumber)) then

                call ExtractDate(Me%ExternalVar%Now,                         &
                                 Year = Year, Month  = Month,  Day    = Day, &
                                 Hour = Hour, Minute = Minute, Second = Second)

                TimeAux(1) = int(Year  )
                TimeAux(2) = int(Month )
                TimeAux(3) = int(Day   )
                TimeAux(4) = int(Hour  )
                TimeAux(5) = int(Minute)
                TimeAux(6) = int(Second)

                !Writes current time
                call ExtractDate   (Me%ExternalVar%Now,                  &
                                    AuxTime(1), AuxTime(2), AuxTime(3),  &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime
                
                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR01'

                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS", &
                                     Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR02'

                PropertyX => Me%FirstProperty

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjHDF5,                                        &
                                      Me%WorkSize2D%ILB,                                 &
                                      Me%WorkSize2D%IUB,                                 &
                                      Me%WorkSize2D%JLB,                                 &
                                      Me%WorkSize2D%JUB,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR03'


PropX:          do while (associated(PropertyX))

                    if(PropertyX%OutputHDF)then

                        if (PropertyX%ID%IDNumber /= SpecificOxygenFlux_ ) then
                            
                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//PropertyX%ID%Name, &
                                                 PropertyX%ID%Name, PropertyX%ID%Units,      &
                                                 Array2D      = PropertyX%Field,             &
                                                 OutputNumber = OutPutNumber,                &
                                                 STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                       &
                                stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR06a'
                           
                         endif
                           
                                                 
                        !Writes everything to disk
                        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR07'

                    endif

                PropertyX => PropertyX%Next

            enddo PropX

            Me%OutPut%NextOutPut = OutPutNumber + 1


        endif  TOut    

        nullify(PropertyX)


    end subroutine OutPut_Results_HDF


    !--------------------------------------------------------------------------
    
    
    subroutine SetSubModulesModifier

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_Property), pointer               :: PropertyX, PropertyY
        real, pointer, dimension(:,:)           :: AtmPressure, Precipitation

        !Begin-----------------------------------------------------------------

        
        if(Me%ExtOptions%HeatFluxYes)then

            call Search_Property(PropertyX, NonSolarFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR10'
            
            endif


            call Search_Property(PropertyX, SurfaceRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR20'
            


            endif

        endif

        if(Me%ExtOptions%OxygenFluxYes)then

            call Search_Property(PropertyX, OxygenFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR30'
            endif

        endif

        if(Me%ExtOptions%CarbonDioxideFluxYes)then

            call Search_Property(PropertyX, CarbonDioxideFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR40'
            endif

        endif


        if(Me%ExtOptions%WQMYes)then

            call Search_Property(PropertyX, SurfaceRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR50'
           
            endif

        endif

        if(Me%ExtOptions%T90VariableYes)then

            call Search_Property(PropertyX, SurfaceRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR60'
           
            endif

        endif

        if(Me%ExtOptions%SurfaceWaterFluxYes)then

            call Search_Property(PropertyX, SurfaceWaterFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceWaterFlux(Me%ObjHydrodynamic,                  &
                                         PropertyX%Field,                     &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR70'
                    
                if (Me%ExtOptions%Precipitation) then
                    
                    call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,         &
                                               Scalar       = Precipitation,            &
                                               ID           = Precipitation_,           &
                                               STAT         = STAT_CALL)
                    if (STAT_CALL .eq. SUCCESS_) then

                        call SetSurfaceFlux(Me%ObjWaterProperties,                      &
                                            Precipitation_,                             &
                                            Precipitation,                              &
                                            STAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_) &
                            stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR75'                
                    else

                        stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR77'                                            
                    
                    endif
                    
                   call UnGetAtmosphere(Me%ObjAtmosphere, Precipitation, &
                                        STAT         = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) &
                        stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR79'
                    
                
                endif                    
            endif

        endif

        if(Me%ExtOptions%HydrodynamicWindYes)then

            call Search_Property(PropertyX, WindStressX_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call Search_Property(PropertyY, WindStressY_, STAT = STAT_CALL) 

                if (STAT_CALL == SUCCESS_)then

                    call SetWindStress (Me%ObjHydrodynamic, PropertyX%FieldGrid, &
                                        PropertyY%FieldGrid, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR80'

                endif

            endif

        endif


        !Use the Atmospheric Pressure for Sea level pressure reference
        if(Me%ExtOptions%HydrodynamicAtmPressureYes)then

            call SetAtmosphericPressure(Me%ObjHydrodynamic,                         &
                                        Me%LocalAtm%AtmosphericPressure%Field,      &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR90'
        
        !Or Use the Mean Sea Level Pressure for Sea level pressure reference
        elseif(Me%ExtOptions%HydrodynamicMslpYes)then

            call SetAtmosphericPressure(Me%ObjHydrodynamic,                         &
                                        Me%LocalAtm%Mslp%Field,      &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR90'
        
        endif
        
        if (Me%IntOptions%WindShearVelocity)then
            
            call SetTurbGOTMWindShearVelocity(TurbGOTMID        = Me%ObjTurbGOTM,           &
                                              WindShearVelocity = Me%WindShearVelocity,     &
                                              STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR100'
            
        endif

#ifndef _WAVES_
        if(Me%ExtOptions%WavesWindYes)then

            call SetWavesWind(WavesID       = Me%ObjWaves,                                      &
                              WindX         = Me%LocalAtm%WindVelocityX%Field,                  &
                              WindY         = Me%LocalAtm%WindVelocityY%Field,                  &
                              STAT          = STAT_CALL)                                    
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR110'

        endif
#endif

i22:    if (Me%ObjLagrangian /= 0) then
#ifndef _LAGRANGIAN_
#ifdef  _LAGRANGIAN_GLOBAL_        
            !Checks Lagrangian 
            call GetLagrangianAirOptionsGlobal(LagrangianID  = Me%ObjLagrangian,                &
                                         Oil           = Me%ExtOptions%OilYes,                  &
                                         Wind          = Me%ExtOptions%LagrangianWindYes,       &
                                         WaterQuality  = Me%ExtOptions%LagrangianWQMYes,        &
                                         T90Variable   = Me%ExtOptions%LagrangianT90Yes,        &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR120'
#endif
#endif

#ifndef _LAGRANGIAN_
        if(Me%ExtOptions%LagrangianWindYes)then     
#ifdef  _LAGRANGIAN_GLOBAL_        
            !Checks Lagrangian 
            call SetLagrangianWindGlobal(LagrangianID = Me%ObjLagrangian,                       &
                                   ModelName    = Me%ModelName,                                 &
                                   WindX        = Me%LocalAtm%WindVelocityX%Field,              &
                                   WindY        = Me%LocalAtm%WindVelocityY%Field,              &
                                   STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR130'
#else
            call SetLagrangianWind(LagrangianID = Me%ObjLagrangian,                       &
                                   WindX        = Me%LocalAtm%WindVelocityX%Field,              &
                                   WindY        = Me%LocalAtm%WindVelocityY%Field,              &
                                   STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR140'
#endif                
        endif

        if(Me%ExtOptions%LagrangianWQMYes  .or. Me%ExtOptions%LagrangianT90Yes)then
            call Search_Property(PropertyX, SurfaceRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then        
#ifdef  _LAGRANGIAN_GLOBAL_ 
                call SetLagSolarRadiationGlobal( LagrangianID   = Me%ObjLagrangian,             &
                                                 ModelName      = Me%ModelName,                 &
                                                 SolarRadiation = PropertyX%Field,              &
                                                 STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                      &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR150'
#else
                call SetLagrangianSolarRadiation(LagrangianID   = Me%ObjLagrangian,             &
                                                 SolarRadiation = PropertyX%Field,              &
                                                 STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                      &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR160'
#endif
            endif
        
        endif

        if(Me%ExtOptions%OilYes)then

            !In this situation, if both Atmospheric Pressure and Mean Sea
            !Level Pressure are present, then MSLP will take precedence over
            !Atmospheric Pressure.
            
            nullify(AtmPressure)
                        
            !Checks if Mslp Field is present
            if(Me%LocalAtm%Mslp%Yes) then

                AtmPressure  => Me%LocalAtm%Mslp%Field

            !Checks if Atmospheric Pressure Field is present
            elseif(Me%LocalAtm%AtmosphericPressure%Yes) then

                AtmPressure  => Me%LocalAtm%AtmosphericPressure%Field

            endif

            if(Me%LocalAtm%Mslp%Yes .or. Me%LocalAtm%AtmosphericPressure%Yes) then
               
#ifdef  _LAGRANGIAN_GLOBAL_ 
                call SetLagrangianAtmPressureGlobal(LagrangianID = Me%ObjLagrangian,    &
                                              ModelName    = Me%ModelName,              &
                                              AtmPressure  = AtmPressure,               &
                                              STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR170'
#else
                call SetLagrangianAtmPressure(LagrangianID = Me%ObjLagrangian,          &
                                              AtmPressure  = AtmPressure,               &
                                              STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR180'
#endif                        
                nullify(AtmPressure)
            endif

        endif
#endif
        
        endif i22
        


    end subroutine SetSubModulesModifier

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillInterfaceWaterAir(ObjInterfaceWaterAirID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjInterfaceWaterAirID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, STAT_CALL, nUsers              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_           
        type (T_Property),    pointer       :: PropertyX

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterfaceWaterAirID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mINTERFACEWATERAIR_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%OutPut%Yes) then

                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR10'
                endif

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR20'
                
                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR30'
                
                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjGridData)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR40'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR50'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjGeometry)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR60'

                nUsers = DeassociateInstance(mMAP_,             Me%ObjMap)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR70'

                nUsers = DeassociateInstance(mHYDRODYNAMIC_,    Me%ObjHydrodynamic)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR80'

                nUsers = DeassociateInstance(mATMOSPHERE_,      Me%ObjAtmosphere)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR90'
                
                nUsers = DeassociateInstance(mWATERPROPERTIES_, Me%ObjWaterProperties)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR100'

#ifndef _LAGRANGIAN_
                if(Me%ObjLagrangian /= 0)then
                    nUsers = DeassociateInstance(mLAGRANGIAN_, Me%ObjLagrangian)
                    if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR110'
                endif
#endif
                if(Me%ObjTurbGOTM /= 0)then
                    nUsers = DeassociateInstance(mTURBGOTM_, Me%ObjTurbGOTM)
                    if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR120'
                endif

#ifndef _WAVES_
                if(Me%ObjWaves /= 0)then
                    nUsers = DeassociateInstance(mWAVES_, Me%ObjWaves)
                    if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR130'
                endif
#endif

                if(Me%Coupled%TimeSerie%Yes)then

                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR140'
                endif


                if(associated(Me%LocalAtm%WindVelocityX%Field      )) &
                   deallocate(Me%LocalAtm%WindVelocityX%Field      )
                if(associated(Me%LocalAtm%WindVelocityY%Field      )) &
                   deallocate(Me%LocalAtm%WindVelocityY%Field      )
                if(associated(Me%LocalAtm%WindDirection%Field      )) &
                   deallocate(Me%LocalAtm%WindDirection%Field      )
                if(associated(Me%LocalAtm%AtmosphericPressure%Field)) &
                   deallocate(Me%LocalAtm%AtmosphericPressure%Field)                
                if(associated(Me%LocalAtm%Mslp%Field)) &    
                   deallocate(Me%LocalAtm%Mslp%Field)   !Deallocate Mean Sea Level Pressure field

                PropertyX => Me%FirstProperty

                do while(associated(PropertyX))
                    
                    deallocate(PropertyX%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR150'
                    nullify(PropertyX%Field)

                    if (associated(PropertyX%FieldGrid)) then
                        deallocate(PropertyX%FieldGrid)
                        nullify   (PropertyX%FieldGrid)
                    endif

                    if(PropertyX%ID%SolutionFromFile)then

                        call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR160'
                    endif


                    PropertyX => PropertyX%Next

                end do

                if (Me%Rugosity%ON) then

                    deallocate(Me%Rugosity%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR170'
                    nullify(Me%Rugosity%Field)

                    if(Me%Rugosity%ID%SolutionFromFile) then

                        call KillFillMatrix(Me%Rugosity%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)&
                            stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR180'
                    endif

                endif

                if(Me%Coupled%BoxTimeSerie%Yes)then
                
                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceSedimentWater - ERR190'

                    deallocate(Me%Scalar2D)
                    nullify   (Me%Scalar2D)
                endif
                
                if (associated(Me%WindShearVelocity))then
                    deallocate(Me%WindShearVelocity)
                    nullify   (Me%WindShearVelocity)
                endif

                
                !Deallocates Instance
                call DeallocateInstance

                ObjInterfaceWaterAirID = 0
                STAT_      = SUCCESS_

            endif
        else 
            STAT_ = ready_
        endif cd1

        if (present(STAT)) STAT = STAT_
           

    end subroutine KillInterfaceWaterAir
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_InterfaceWaterAir), pointer          :: AuxObjInterfaceWaterAir
        type (T_InterfaceWaterAir), pointer          :: PreviousObjWaterSedInterface

        !Updates pointers
        if (Me%InstanceID == FirstObjInterfaceWaterAir%InstanceID) then
            FirstObjInterfaceWaterAir           => FirstObjInterfaceWaterAir%Next
        else
            PreviousObjWaterSedInterface        => FirstObjInterfaceWaterAir
            AuxObjInterfaceWaterAir             => FirstObjInterfaceWaterAir%Next
            do while (AuxObjInterfaceWaterAir%InstanceID /= Me%InstanceID)
                PreviousObjWaterSedInterface    => AuxObjInterfaceWaterAir
                AuxObjInterfaceWaterAir         => AuxObjInterfaceWaterAir%Next
            enddo

            !Now update linked list
            PreviousObjWaterSedInterface%Next   => AuxObjInterfaceWaterAir%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjInterfaceWaterAir_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceWaterAir_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjInterfaceWaterAir_ID > 0) then
            call LocateObjInterfaceWaterAir (ObjInterfaceWaterAir_ID)
            ready_ = VerifyReadLock (mINTERFACEWATERAIR_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        endif cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjInterfaceWaterAir (ObjInterfaceWaterAirID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceWaterAirID

        !Local-----------------------------------------------------------------

        Me => FirstObjInterfaceWaterAir
        do while (associated (Me))
            if (Me%InstanceID == ObjInterfaceWaterAirID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleInterfaceWaterAir - LocateObjInterfaceWaterAir - ERR01'

    end subroutine LocateObjInterfaceWaterAir

    !--------------------------------------------------------------------------

end module ModuleInterfaceWaterAir

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------






