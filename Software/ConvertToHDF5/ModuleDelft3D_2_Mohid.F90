!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Delft3D_2_Mohid
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Paulo Leitão - v4.0
! DESCRIPTION   : Module to serve as Delft3D_2_Mohid to create new modules
!
!------------------------------------------------------------------------------


Module ModuleDelft3D_2_Mohid

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleHDF5
    use ModuleStopWatch,      only : CreateWatchGroup, KillWatchGroup    
    use ModuleFunctions
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap
    use ModuleGeometry
    use ModuleMap
    use ModuleField4D

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertDelft3D_2_Mohid

    !Destructor
    private :: KillDelft3D_2_Mohid                                                     

    !Management
    private ::      Ready
    private ::          LocateObjDelft3D_2_Mohid 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetDelft3D_2_Mohid3D_I
    private :: UnGetDelft3D_2_Mohid3D_R8
    interface  UnGetDelft3D_2_Mohid
        module procedure UnGetDelft3D_2_Mohid3D_I
        module procedure UnGetDelft3D_2_Mohid3D_R8
    end interface  UnGetDelft3D_2_Mohid

    !Types---------------------------------------------------------------------
    
    private :: T_Delft3D_2_Mohid
    type       T_Delft3D_2_Mohid
        integer                                             :: InstanceID
        
        type (T_Size3D)                                     :: Size, WorkSize
        real(8), dimension(:, :, :),  pointer               :: Matrix
        
        type(T_Time)                                        :: BeginTime, EndTime, TimeRef
        real                                                :: DT

        logical                                             :: Atm_ON
        logical                                             :: Ocean_IN
        logical                                             :: Ocean_IN_Ini        
        logical                                             :: Ocean_IN_Bound        
        !testing
        logical                                             :: Ocean_Read_Hydro
        
                       
        integer                                             :: AtmNumber
        character(len=PathLength  ), dimension(:), pointer  :: AtmPropFile
        character(len=StringLength), dimension(:), pointer  :: AtmPropNameOut        
        type (T_PropertyID),         dimension(:), pointer  :: AtmProp
        character(len=StringLength), dimension(:), pointer  :: AtmOutFormat
        
        integer                                             :: AtmOutUnit        

        type (T_Size2D)                                     :: AtmWorkSize2D
        real,       dimension(:),   pointer                 :: AtmX                        
        real,       dimension(:),   pointer                 :: AtmY
        real,       dimension(:),   pointer                 :: AtmPropXY        
        real,       dimension(:,:), pointer                 :: AtmProp2D        
        logical,    dimension(:),   pointer                 :: AtmNoDataXY
        integer,    dimension(:),   pointer                 :: ObjField4D_Atm   

        real                                                :: AtmNODATA = 999.999
        real                                                :: AtmXorig
        real                                                :: AtmDx
        real                                                :: AtmYorig
        real                                                :: AtmDy

        character(len=PathLength)                           :: FileName_Atm_Grid_Out
        character(len=PathLength)                           :: FileName_Mod_Grid_Out
        character(len=PathLength)                           :: FileName_Mod_Geo_Out
        character(len=PathLength  )                         :: ModOutFile
        character(len=StringLength)                         :: ModIniFormat
        character(len=StringLength)                         :: ModBoundFormat
        integer                                             :: ModIniOutUnit
        logical                                             :: BoundAstro
                        
        integer                                             :: ClientNumber

        character(len=StringLength)                         :: SSHPropNameOut        
        type (T_PropertyID)                                 :: SSHProp
        real                                                :: SSHValueIni
        
        character(len=StringLength)                         :: VelXPropNameOut        
        type (T_PropertyID)                                 :: VelXProp
        real                                                :: VelXValueIni
        
        character(len=StringLength)                         :: VelYPropNameOut        
        type (T_PropertyID)                                 :: VelYProp
        real                                                :: VelYValueIni        

        character(len=StringLength)                         :: SalPropNameOut        
        type (T_PropertyID)                                 :: SalProp
        real                                                :: SalValueIni

        character(len=StringLength)                         :: TempPropNameOut        
        type (T_PropertyID)                                 :: TempProp
        real                                                :: TempValueIni        

        integer,    dimension(:,:,:), pointer               :: ModWaterPoints3D
        type (T_Size3D)                                     :: ModWorkSize3D, ModSize3D
        
        real,       dimension(:),     pointer               :: ModX3D        
        real,       dimension(:),     pointer               :: ModY3D        
        real,       dimension(:),     pointer               :: ModZ3D        
        real,       dimension(:),     pointer               :: ModPropXY3D                
        real,       dimension(:,:,:), pointer               :: ModProp3D   
        logical,    dimension(:),     pointer               :: ModNoDataXY3D 
        
        real,       dimension(:),     pointer               :: ModX2D        
        real,       dimension(:),     pointer               :: ModY2D        
        real,       dimension(:),     pointer               :: ModPropXY2D                
        real,       dimension(:,:),   pointer               :: ModProp2D   
        logical,    dimension(:),     pointer               :: ModNoDataXY2D   
        
        real                                                :: ModNODATA = 0.
        
        integer                                             :: ModNcells2D, ModNcells3D
        
        integer                                             :: BoundCellsNumber
        integer                                             :: BoundSectionsNumber    
        character(len=StringLength), dimension(:), pointer  :: BoundSectionsName
        integer                                             :: ModNLayers
        integer                                             :: ModNinstants     
        
        integer                                             :: BoundCellsUnit
        character(len=PathLength)                           :: BoundCellsFile        
        character(len=PathLength)                           :: BoundRiemannFile
        character(len=PathLength)                           :: BoundDensityFile
        real,       dimension(:),     pointer               :: BoundTime        
           

        integer,    dimension(:),     pointer               :: BoundCellsI
        integer,    dimension(:),     pointer               :: BoundCellsJ        
        real,       dimension(:),     pointer               :: BoundX2D        
        real,       dimension(:),     pointer               :: BoundY2D        
        real,       dimension(:),     pointer               :: BoundPropXY2D
        logical,    dimension(:),     pointer               :: BoundNoDataXY2D         
       
        !point XY, time
        real,       dimension(:,:),   pointer               :: BoundSSH_XYT2D
        real,       dimension(:,:),   pointer               :: BoundSSH_ASTRO_XYT2D        
        
       
        real,       dimension(:),     pointer               :: BoundX3D        
        real,       dimension(:),     pointer               :: BoundY3D        
        real,       dimension(:),     pointer               :: BoundZ3D        
        real,       dimension(:),     pointer               :: BoundPropXYZ3D
        logical,    dimension(:),     pointer               :: BoundNoDataXYZ3D         
        !point XY, layer, time
        real,       dimension(:,:,:), pointer               :: BoundSalXYZT3D   
        real,       dimension(:,:,:), pointer               :: BoundTempXYZT3D           
        real,       dimension(:,:,:), pointer               :: BoundVelx3D_XYZT3D 
        real,       dimension(:,:,:), pointer               :: BoundVely3D_XYZT3D    
        real,       dimension(:,:,:), pointer               :: BoundVelx3D_Astro_XYZT3D 
        real,       dimension(:,:,:), pointer               :: BoundVely3D_Astro_XYZT3D         
        real,       dimension(:,:,:), pointer               :: BoundRiemannXYZT3
        

        integer                                             :: ObjEnterData         = 0
        integer                                             :: ObjTime              = 0
        integer                                             :: ObjAtm_Grid_Out      = 0
        integer                                             :: ObjMod_Grid_Out      = 0
        integer                                             :: ObjMod_Bathym_Out    = 0        
        integer                                             :: ObjMod_Map2D_Out     = 0                
        integer                                             :: ObjMod_Geo_Out       = 0
        integer                                             :: ObjMod_Map3D_Out     = 0

        integer                                             :: ObjField4D_Ssh       = 0
        integer                                             :: ObjField4D_Ssh_Astro = 0       
        integer                                             :: ObjField4D_VelX      = 0        
        integer                                             :: ObjField4D_VelY      = 0
        integer                                             :: ObjField4D_VelX_Astro= 0        
        integer                                             :: ObjField4D_VelY_Astro= 0
        integer                                             :: ObjField4D_Sal       = 0         
        integer                                             :: ObjField4D_Temp      = 0  
        
        type(T_Delft3D_2_Mohid), pointer                    :: Next
        
    end type  T_Delft3D_2_Mohid

    !Global Module Variables
    type (T_Delft3D_2_Mohid), pointer                         :: FirstObjDelft3D_2_Mohid
    type (T_Delft3D_2_Mohid), pointer                         :: Me

    integer                                                   :: mDelft3D_2_Mohid_ = 0 !just to compile

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertDelft3D_2_Mohid(ObjEnterData, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjEnterData 
        integer                                         :: ClientNumber        
        integer, optional, intent(OUT)                  :: STAT     


        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, nUsers

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, ObjEnterData)

        Me%ClientNumber = ClientNumber
        
        call SetDate(Me%TimeRef, 2018, 1, 1, 0, 0, 0)        
            
        call ReadGlobalOptions
        
        call ConstructTime
        
        if (Me%Atm_ON) call ConstructAtmInOut
        
        if (Me%Ocean_IN) then            
        
            call ConstructModGrid        
            
        endif
        
        if (Me%Ocean_IN) then
            call ConstructModOceanIn     
        endif
        
        if (Me%Atm_ON  ) call WriteAtmInOut

        if (Me%Ocean_IN) call WriteModInOut
        
        if (Me%Ocean_Read_Hydro) then
            call ReadHydro           
        endif               
        
        nUsers = DeassociateInstance(mENTERDATA_, ObjEnterData)
        if (nUsers == 0) stop 'ConvertDelft3D_2_Mohid - ModuleDelft3D_2_Mohid - ERR10'        

        STAT_ = SUCCESS_
        
        if (present(STAT)) STAT = STAT_


        !----------------------------------------------------------------------

    end subroutine ConvertDelft3D_2_Mohid
 

    !--------------------------------------------------------------------------
    
    subroutine ConstructModOceanIn
    
        !Local-------------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
    
    
        !Begin-------------------------------------------------------------------    
        

        call ConstructModProp( block_begin  = '<<begin_ssh>>',                          &
                               block_end    = '<<end_ssh>>',                            &
                               ObjField4D   = Me%ObjField4D_Ssh,                        & 
                               PropID       = Me%SSHProp,                               &
                               ValueIni     = Me%SSHValueIni)

        call ConstructModProp( block_begin  = '<<begin_velx>>',                         &
                               block_end    = '<<end_velx>>',                           &
                               ObjField4D   = Me%ObjField4D_VelX,                       & 
                               PropID       = Me%VelXProp,                              &
                               ValueIni     = Me%VelxValueIni)
                              
        call ConstructModProp( block_begin  = '<<begin_vely>>',                         &
                               block_end    = '<<end_vely>>',                           &
                               ObjField4D   = Me%ObjField4D_VelY,                       & 
                               PropID       = Me%VelYProp,                              &
                               ValueIni     = Me%VelyValueIni)                                      

        call ConstructModProp( block_begin  = '<<begin_temp>>',                         &
                               block_end    = '<<end_temp>>',                           &
                               ObjField4D   = Me%ObjField4D_Temp,                       & 
                               PropID       = Me%TempProp,                              &
                               ValueIni     = Me%TempValueIni)     
                               
        call ConstructModProp( block_begin  = '<<begin_sal>>',                          &
                               block_end    = '<<end_sal>>',                            &
                               ObjField4D   = Me%ObjField4D_Sal,                        & 
                               PropID       = Me%SalProp,                               &
                               ValueIni     = Me%SalValueIni)       
                               

        
        if (Me%Ocean_IN_Ini) then
    
            call GetData(Me%ModIniFormat,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'MOD_INI_FORMAT',                               &
                         default      = '(12(f15.8,2X))',                               &
                         ClientModule = 'ModuleDelft3D_2_Mohid',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR10'
            endif   
            
        endif
        
        Me%BoundAstro = .false. 
        
        if (Me%Ocean_IN_Bound) then        
        
            call GetData(Me%BoundAstro,                                                 &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOUND_ASTRO',                                  &
                         ClientModule = 'ModuleDelft3D_2_Mohid',                        &
                         default      = .false.,                                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR30'
            endif   


            call GetData(Me%BoundCellsFile,                                             &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOUND_CELLS_FILE',                             &
                         ClientModule = 'ModuleDelft3D_2_Mohid',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR40'
            endif           
            
            if (iflag == 0) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR50'
            endif        

            call GetData(Me%BoundRiemannFile,                                           &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOUND_RIEMANN_FILE',                           &
                         ClientModule = 'ModuleDelft3D_2_Mohid',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR60'
            endif           
            
            if (iflag == 0) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR70'
            endif                

            call GetData(Me%BoundDensityFile,                                           &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOUND_DENSITY_FILE',                           &
                         ClientModule = 'ModuleDelft3D_2_Mohid',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR80'
            endif           
            
            if (iflag == 0) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR90'
            endif               
            
            call GetData(Me%ModBoundFormat,                                             &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'MOD_BOUND_FORMAT',                             &
                         default      = '(e14.5)',                                      & 
                         ClientModule = 'ModuleDelft3D_2_Mohid',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructModOceanIn - ModuleDelft3D_2_Mohid - ERR100'
            endif   
            


            if (Me%BoundAstro) then

                call ConstructModProp( block_begin  = '<<begin_ssh_astro>>',            &
                                       block_end    = '<<end_ssh_astro>>',              &
                                       ObjField4D   = Me%ObjField4D_Ssh_Astro,          & 
                                       PropID       = Me%SSHProp)

                call ConstructModProp( block_begin  = '<<begin_velx_astro>>',           &
                                       block_end    = '<<end_velx_astro>>',             &
                                       ObjField4D   = Me%ObjField4D_VelX_Astro,         & 
                                       PropID       = Me%VelXProp)
                                      
                call ConstructModProp( block_begin  = '<<begin_vely_astro>>',           &
                                       block_end    = '<<end_vely_astro>>',             &
                                       ObjField4D   = Me%ObjField4D_VelY_Astro,         & 
                                       PropID       = Me%VelYProp)                                             
            endif    

        endif            
            
    end subroutine ConstructModOceanIn                                                                  

    !--------------------------------------------------------------------------
                
    
    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%Atm_ON,                                                         &    
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'ATM_IN',                                           &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleDelft3D_2_Mohid',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadGlobalOptions - ModuleDelft3D_2_Mohid - ERR10'
        endif            
        
        call GetData(Me%Ocean_IN,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OCEAN_IN',                                         &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleDelft3D_2_Mohid',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadGlobalOptions - ModuleDelft3D_2_Mohid - ERR20'
        endif            
        
        if (Me%Ocean_IN) then

            call GetData(Me%Ocean_IN_Ini,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'OCEAN_IN_INI',                                 &
                         default      = .true.,                                         &
                         ClientModule = 'ModuleDelft3D_2_Mohid',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                stop 'ReadGlobalOptions - ModuleDelft3D_2_Mohid - ERR30'
            endif            
            
            call GetData(Me%Ocean_IN_Bound,                                             &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'OCEAN_IN_BOUND',                               &
                         default      = .true.,                                         &
                         ClientModule = 'ModuleDelft3D_2_Mohid',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                stop 'ReadGlobalOptions - ModuleDelft3D_2_Mohid - ERR40'
            endif                 
        
        endif

        call GetData(Me%Ocean_Read_Hydro,                                               &
                        Me%ObjEnterData, iflag,                                         &
                        SearchType   = FromBlock,                                       &
                        keyword      = 'OCEAN_READ_HYDRO',                              &
                        default      = .false.,                                         &
                        ClientModule = 'ModuleDelft3D_2_Mohid',                         &
                        STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadGlobalOptions - ModuleDelft3D_2_Mohid - ERR50'
        endif            
        
        
        
    end subroutine ReadGlobalOptions

!--------------------------------------------------------------------------

    subroutine ConstructTime

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: Dummy

        !Begin-----------------------------------------------------------------


        call ReadTimeKeyWords(ObjEnterData      = Me%ObjEnterData,                      &
                              ExtractTime       = FromBlock,                            &
                              BeginTime         = Me%BeginTime,                         &
                              EndTime           = Me%EndTime,                           &
                              DT                = Me%DT,                                &
                              VariableDT        = Dummy,                                &
                              ClientModule      = "ModuleDelft3D_2_Mohid")   


        call StartComputeTime(TimeID           = Me%ObjTime,                            &
                              InitialSystemTime= Me%BeginTime,                          &
                              BeginTime        = Me%BeginTime,                          &
                              EndTime          = Me%EndTime,                            &
                              DT               = Me%DT,                                 & 
                              VariableDT       = .false.,                               &
                              STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTime - ModuleDelft3D_2_Mohid - ERR10'
        

    end subroutine ConstructTime

    !--------------------------------------------------------------------------    
    
    subroutine ConstructAtmInOut

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag, i
        real                                        :: LatReference, LonReference
        logical                                     :: BlockInBlockFound
        character(len=PathLength)                   :: FilenameIn
        

        !Begin-----------------------------------------------------------------
        

        call GetData(Me%FileName_Atm_Grid_Out,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'ATM_GRID_OUT',                                     &
                     ClientModule = 'ModuleDelft3D_2_Mohid',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR10'
        endif                        
        
        if (iflag     == 0) then
            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR20'
        endif        
        
        call ConstructHorizontalGrid(HorizontalGridID = Me%ObjAtm_Grid_Out,             &
                                     DataFile         = Me%FileName_Atm_Grid_Out,       &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR30'
        endif 
        

        call GetLatitudeLongitude(Me%ObjAtm_Grid_Out, Latitude  = LatReference,         &
                                                      Longitude = LonReference,         &  
                                                      STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR40'
        endif         
        

        call GetData(Me%AtmNumber,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'ATM_PROP_NUMBER',                                  &
                     ClientModule = 'ModuleDelft3D_2_Mohid',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR50'
        endif
        
        if (iflag     == 0) then
            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR60'
        endif
                
        
        allocate (Me%AtmPropFile       (1:Me%AtmNumber))
        allocate (Me%AtmProp           (1:Me%AtmNumber))
        allocate (Me%AtmPropNameOut    (1:Me%AtmNumber))
        allocate (Me%AtmOutFormat      (1:Me%AtmNumber))
        allocate (Me%ObjField4D_Atm    (1:Me%AtmNumber))        
        
        Me%ObjField4D_Atm(:) = 0
        
d1:     do i=1, Me%AtmNumber
            call ExtractBlockFromBlock (Me%ObjEnterData,                                &
                                        ClientNumber        = Me%ClientNumber,          &
                                        block_begin         = "<<begin_atmosphere>>",   &
                                        block_end           = "<<end_atmosphere>>",     &
                                        BlockInBlockFound   = BlockInBlockFound,        &
                                        STAT                = STAT_CALL)
cd1 :       if (STAT_CALL .EQ. SUCCESS_) then    
cd2 :           if (BlockInBlockFound) then  

                    call GetData(FilenameIn,                                            &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'FILENAME_IN',                          &
                                 ClientModule = 'ModuleDelft3D_2_Mohid',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ConstructVelXInOut - ModuleDelft3D_2_Mohid - ERR70'
                    endif    

                    call ConstructPropertyID (Me%AtmProp(i), Me%ObjEnterData,           &
                                              ExtractType = FromBlockInBlock)                                                
                        
                    call ConstructField4D(Field4DID     = Me%ObjField4D_Atm(i),         &
                                          EnterDataID   = Me%ObjEnterData,              &
                                          ExtractType   = FromBlockInBlock,             &
                                          TimeID        = Me%ObjTime,                   &   
                                          FileName      = FilenameIn,                   &  
                                          LatReference  = LatReference,                 &
                                          LonReference  = LonReference,                 & 
                                          Extrapolate   = .true.,                       &        
                                          PropertyID    = Me%AtmProp(i),                &
                                          ClientID      = Me%ClientNumber,              &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) then                                          
                        stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR80'
                    endif                        
                    
                    call GetData(Me%AtmPropFile(i),                                     &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'ATM_PROP_FILE',                        &
                                 ClientModule = 'ModuleDelft3D_2_Mohid',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR90'
                    endif                    
                                   
                    call GetData(Me%AtmPropNameOut(i),                                  &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'ATM_PROP_NAME_OUT',                    &
                                 ClientModule = 'ModuleDelft3D_2_Mohid',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR100'
                    endif                                            
                                
                    call GetData(Me%AtmOutFormat(i),                                    &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'ATM_PROP_FORMAT_OUT',                  &
                                 ClientModule = 'ModuleDelft3D_2_Mohid',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR110'
                    endif                                   
                                
                else cd2        
                    
                    stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR120'                
                    
                endif cd2
            else cd1
                stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR130'
            endif cd1                

        enddo d1
        
        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR140'
        endif     

!        call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)         
!        if (STAT_CALL /= SUCCESS_) then
!            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR130'
!        endif        

    end subroutine ConstructAtmInOut


    !--------------------------------------------------------------------------    

   subroutine ConstructModGrid

        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer             :: SurfaceElevation
        integer                                     :: STAT_CALL
        integer                                     :: iflag
       

        !Begin-----------------------------------------------------------------
            
        call GetData(Me%FileName_Mod_Grid_Out,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MOD_GRID_OUT',                                     &
                     ClientModule = 'ModuleDelft3D_2_Mohid',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR10'
        endif                        
        
        
        call GetData(Me%FileName_Mod_Geo_Out,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MOD_GEOMETRY_OUT',                                 &
                     ClientModule = 'ModuleDelft3D_2_Mohid',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR40'
        endif             

        call GetData(Me%ModOutFile,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MOD_FILENAME_OUT',                                 &
                     ClientModule = 'ModuleDelft3D_2_Mohid',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR50'
        endif         
        
        call ConstructHorizontalGrid(HorizontalGridID = Me%ObjMod_Grid_Out,             &
                                     DataFile         = Me%FileName_Mod_Grid_Out,       &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR60'
        endif            

        call ConstructGridData      (GridDataID       = Me%ObjMod_Bathym_Out,           &
                                     HorizontalGridID = Me%ObjMod_Grid_Out,             &
                                     FileName         = Me%FileName_Mod_Grid_Out,       &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR70'
        endif   

        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjMod_Map2D_Out,            &
                                     GridDataID       = Me%ObjMod_Bathym_Out,           &
                                     HorizontalGridID = Me%ObjMod_Grid_Out,             &
                                     ActualTime       = Me%BeginTime,                   & 
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR80'
        endif   

        call ConstructGeometry      (GeometryID       = Me%ObjMod_Geo_Out,              &
                                     GridDataID       = Me%ObjMod_Bathym_Out,           &
                                     HorizontalGridID = Me%ObjMod_Grid_Out,             &
                                     HorizontalMapID  = Me%ObjMod_Map2D_Out,            &
                                     ActualTime       = Me%BeginTime,                   &
                                     NewDomain        = Me%FileName_Mod_Geo_Out,        &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR90'
        endif
        

        call GetGeometrySize(GeometryID     = Me%ObjMod_Geo_Out,                        &
                             Size           = Me%ModSize3D,                             &
                             WorkSize       = Me%ModWorkSize3D,                         &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR100'
        endif

        call ConstructMap ( Map_ID          = Me%ObjMod_Map3D_Out,                      &
                            GeometryID      = Me%ObjMod_Geo_Out,                        &
                            HorizontalMapID = Me%ObjMod_Map2D_Out,                      &
                            TimeID          = Me%ObjTime,                               &
                            STAT            = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR110'
        endif

        allocate(SurfaceElevation (Me%ModSize3D%ILB:Me%ModSize3D%IUB,                   &
                                   Me%ModSize3D%JLB:Me%ModSize3D%JUB))
        SurfaceElevation(:,:) = 0

        call GetWaterPoints3D(Me%ObjMod_Map3D_Out, Me%ModWaterPoints3D, STAT = STAT_CALL) 
        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR120'
        endif

        call ComputeInitialGeometry(GeometryID      = Me%ObjMod_Geo_Out,                &
                                    WaterPoints3D   = Me%ModWaterPoints3D,              &
                                    SurfaceElevation= SurfaceElevation,                 &
                                    ActualTime      = Me%BeginTime,                     &
                                    STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModGrid - ModuleDelft3D_2_Mohid - ERR130'
        endif

        deallocate(SurfaceElevation)
        

    end subroutine ConstructModGrid        
    
    !--------------------------------------------------------------------------        
    
    subroutine ConstructModProp(block_begin, block_end, ObjField4D, PropID, ValueIni)
        
        !Arguments--------------------------------------------------------------
        character(len=*)                            :: block_begin
        character(len=*)                            :: block_end
        integer                                     :: ObjField4D
        type (T_PropertyID)                         :: PropID        
        real, optional                              :: ValueIni

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        real                                        :: LatReference, LonReference, ValueIni_
        logical                                     :: BlockInBlockFound
        character(len=PathLength)                   :: FilenameIn
        real, dimension(1:2,1:2)                    :: WindowLimitsXY        
        real                                        :: West, East, South, North         

        !Begin-----------------------------------------------------------------

        call GetLatitudeLongitude(Me%ObjMod_Grid_Out, Latitude  = LatReference,         &
                                                      Longitude = LonReference,         &  
                                                      STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModProp - ModuleDelft3D_2_Mohid - ERR10'
        endif         
        
        call GetGridBorderLimits(Me%ObjMod_Grid_Out, West, East, South, North, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModProp - ModuleDelft3D_2_Mohid - ERR20'
        endif            
        
        WindowLimitsXY(2,1) = South
        WindowLimitsXY(2,2) = North
        WindowLimitsXY(1,1) = West
        WindowLimitsXY(1,2) = East        

        call ExtractBlockFromBlock (Me%ObjEnterData,                                    &
                                    ClientNumber        = Me%ClientNumber,              &
                                    block_begin         = trim(block_begin),            &
                                    block_end           = trim(block_end),              &
                                    BlockInBlockFound   = BlockInBlockFound,            &
                                    STAT                = STAT_CALL)
cd1 :   if (STAT_CALL .EQ. SUCCESS_) then    
cd2 :       if (BlockInBlockFound) then  

                call GetData(FilenameIn,                                                &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'FILENAME_IN',                              &
                             ClientModule = 'ModuleDelft3D_2_Mohid',                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) then
                    stop 'ConstructModProp - ModuleDelft3D_2_Mohid - ERR20'
                endif    
                
                if (iflag == 0) then
                    stop 'ConstructModProp - ModuleDelft3D_2_Mohid - ERR30'
                endif
                

                call ConstructPropertyID (PropID, Me%ObjEnterData,                      &
                                          ExtractType = FromBlockInBlock)                                                    
                
                call ConstructField4D(Field4DID     = ObjField4D,                       &
                                      EnterDataID   = Me%ObjEnterData,                  &
                                      ExtractType   = FromBlockInBlock,                 &
                                      TimeID        = Me%ObjTime,                       &   
                                      FileName      = FilenameIn,                       &
                                      LatReference  = LatReference,                     &
                                      LonReference  = LonReference,                     &
                                      WindowLimitsXY= WindowLimitsXY,                   &                                      
                                      Extrapolate   = .true.,                           &        
                                      PropertyID    = PropID,                           &
                                      ClientID      = Me%ClientNumber,                  &
                                      STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then                                          
                    stop 'ConstructModProp - ModuleDelft3D_2_Mohid - ERR40'
                endif                        
                
                call GetData(ValueIni_,                                             &
                             Me%ObjEnterData, iflag,                                &
                             SearchType   = FromBlockInBlock,                       &
                             keyword      = 'VALUE_INI',                            &
                             default      =  FillValueReal,                         &
                             ClientModule = 'ModuleDelft3D_2_Mohid',                &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) then
                    stop 'ConstructVelXInOut - ModuleDelft3D_2_Mohid - ERR50'
                endif    
                
                if (present(ValueIni)) ValueIni = ValueIni_
                
            else cd2        
                
                stop 'ConstructModProp - ModuleDelft3D_2_Mohid - ERR70'
                
            endif cd2
        else cd1
            stop 'ConstructModProp - ModuleDelft3D_2_Mohid - ERR120'
        endif cd1           
        
        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructModProp - ModuleDelft3D_2_Mohid - ERR130'
        endif            

        
!        call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)         
!        if (STAT_CALL /= SUCCESS_) then
!            stop 'ConstructAtmInOut - ModuleDelft3D_2_Mohid - ERR130'
!        endif        

    end subroutine ConstructModProp


    !--------------------------------------------------------------------------        

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine GetDelft3D_2_MohidPointer (ObjDelft3D_2_MohidID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjDelft3D_2_MohidID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDelft3D_2_MohidID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mDelft3D_2_Mohid_, Me%InstanceID)

            Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetDelft3D_2_MohidPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetDelft3D_2_MohidInteger (ObjDelft3D_2_MohidID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjDelft3D_2_MohidID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDelft3D_2_MohidID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetDelft3D_2_MohidInteger

    !--------------------------------------------------------------------------

    subroutine UnGetDelft3D_2_Mohid3D_I(ObjDelft3D_2_MohidID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjDelft3D_2_MohidID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDelft3D_2_MohidID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mDelft3D_2_Mohid_, Me%InstanceID, "UnGetDelft3D_2_Mohid3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetDelft3D_2_Mohid3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetDelft3D_2_Mohid3D_R8(ObjDelft3D_2_MohidID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjDelft3D_2_MohidID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDelft3D_2_MohidID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mDelft3D_2_Mohid_, Me%InstanceID,  "UnGetDelft3D_2_Mohid3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetDelft3D_2_Mohid3D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine WriteAtmInOut


        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, STAT_CALL
        type (T_Time)                               :: CurrentTime

        !----------------------------------------------------------------------
        
        call ReadAtmGridXYZ
        
        do i = 1, Me%AtmNumber        
            
            call OpenAtmOutFile(i)
            
            CurrentTime = Me%BeginTime
            
            do while (CurrentTime<=Me%EndTime)
            
                call ReadProp2D (CurrentTime, i)
            
                call WriteProp2D(CurrentTime, i)
                
                CurrentTime = CurrentTime + Me%DT
            
            enddo
            
            call UnitsManager(Me%AtmOutUnit, CLOSE_FILE, STAT = STAT_CALL) 
                                       
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteAtmInOut - ModuleDelft3D_2_Mohid - ERR10'
            endif              

            call KillField4D(Field4DID = Me%ObjField4D_Atm(i),                               &
                             STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteAtmInOut - ModuleDelft3D_2_Mohid - ERR10'
            endif              

            
        enddo
        
        
        
        deallocate(Me%AtmPropFile   ) 
        deallocate(Me%AtmPropNameOut)        
        deallocate(Me%AtmProp       )
        deallocate(Me%AtmOutFormat  )
        deallocate(Me%AtmX          )              
        deallocate(Me%AtmY          )
        deallocate(Me%AtmPropXY     )        
        deallocate(Me%AtmProp2D     )        
        deallocate(Me%AtmNoDataXY   )
        
        


    end subroutine WriteAtmInOut
    
    !---------------------------------------------------------------------------
    
    
    subroutine ReadAtmGridXYZ

        !Arguments-------------------------------------------------------------
                   
        !Local-----------------------------------------------------------------
        real, dimension(:,:),  pointer              :: CoordX, CoordY        
        type (T_Size2D)                             :: Size2D
        integer                                     :: iSizeVector, icount, i, j, STAT_CALL
        
        !----------------------------------------------------------------------

        call GetHorizontalGridSize(HorizontalGridID = Me%ObjAtm_Grid_Out,               &
                                   Size             = Size2D,                           &
                                   WorkSize         = Me%AtmWorkSize2D,                 &
                                   STAT             = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR10' 
        endif
        
        
        iSizeVector = (Me%AtmWorkSize2D%JUB - Me%AtmWorkSize2D%JLB + 1) *               &
                      (Me%AtmWorkSize2D%IUB - Me%AtmWorkSize2D%ILB + 1)
        
        
        allocate(Me%AtmX           (1:iSizeVector))
        allocate(Me%AtmY           (1:iSizeVector))
        allocate(Me%AtmPropXY      (1:iSizeVector))
        allocate(Me%AtmNoDataXY    (1:iSizeVector))
        
        
        allocate(Me%AtmProp2D(Size2D%ILB:Size2D%IUB,Size2D%JLB:Size2D%JUB))        
        

        call GetZCoordinates(Me%ObjAtm_Grid_Out, CoordX, CoordY, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR20' 
        endif
        
        icount = 0
        
        do j = Me%AtmWorkSize2D%JLB, Me%AtmWorkSize2D%JUB
        do i = Me%AtmWorkSize2D%ILB, Me%AtmWorkSize2D%IUB
        
            icount          = icount + 1
            Me%AtmX(icount) = CoordX(i, j)
            Me%AtmY(icount) = CoordY(i, j)
        
        enddo
        enddo
        
        Me%AtmDx = CoordX(1, 2) - CoordX(1, 1)
        Me%AtmDy = CoordY(2, 1) - CoordY(1, 1)        
        

        call GetGridOrigin(Me%ObjAtm_Grid_Out, Me%AtmXorig, Me%AtmYorig, STAT  = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR30' 
        endif
        
        
        call UnGetHorizontalGrid(Me%ObjAtm_Grid_Out, CoordX, STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR40' 
        endif        

        call UnGetHorizontalGrid(Me%ObjAtm_Grid_Out, CoordY, STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR50' 
        endif        

    
    end subroutine ReadAtmGridXYZ    

    !---------------------------------------------------------------------------
    
    subroutine OpenAtmOutFile(i)

        !Arguments-------------------------------------------------------------
        integer                                     :: i                   
        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: CharOut 
        integer                                     :: STAT_CALL, n_cols
        
        !----------------------------------------------------------------------

        call UnitsManager(Me%AtmOutUnit, OPEN_FILE, STAT = STAT_CALL) 
                                   
        if (STAT_CALL /= SUCCESS_) then
            stop 'OpenAtmOutFile - ModuleDelft3D_2_Mohid - ERR10'
        endif             
        
        open(UNIT   = Me%AtmOutUnit,                                                    &
             FILE   = trim(adjustl(Me%AtmPropFile(i))),                                 &
             FORM   = "FORMATTED",                                                      &   
             IOSTAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop 'OpenAtmOutFile - ModuleDelft3D_2_Mohid - ERR10'
        endif      
        
        n_cols = Me%AtmWorkSize2D%JUB - Me%AtmWorkSize2D%JLB + 1
        
        write(CharOut,'(I7)') n_cols
        Me%AtmOutFormat(i) = "("//trim(adjustl(CharOut))//trim(Me%AtmOutFormat(i))//")"
        
        write(CharOut, '(A23)'      ) 'FileVersion      = 1.03'
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))
        write(CharOut, '(A44)'      ) 'Filetype         = meteo_on_equidistant_grid'
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A19,I7)'   ) 'n_cols           = ', n_cols
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A19,I7)'   ) 'n_rows           = ', Me%AtmWorkSize2D%IUB
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A25)'      ) 'grid_unit        = degree'
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A19,f12.6)') 'x_llcenter       = ', Me%AtmXorig
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A19,f12.6)') 'dx               = ', Me%AtmDx
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A19,f12.6)') 'y_llcenter       = ', Me%AtmYorig
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A19,f12.6)') 'dy               = ', Me%AtmDy
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A19,f12.6)') 'NODATA_value     = ', Me%AtmNoData
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A20)'      ) 'n_quantity       = 1'
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A)'        ) 'quantity1        = '// trim(adjustl(Me%AtmPropNameOut(i)))
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        
        write(CharOut, '(A)'        ) 'unit1            = '// trim(adjustl(Me%AtmProp(i)%Units))  
        write(Me%AtmOutUnit,'(A)'   ) adjustl(trim(CharOut))        

    
    end subroutine OpenAtmOutFile    
    

    !---------------------------------------------------------------------------
    
    !---------------------------------------------------------------------------
    
    
    subroutine ReadModGridXYZ

        !Arguments-------------------------------------------------------------
                   
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer             :: ZCellCenter 
        real, dimension(:,:),  pointer              :: CoordX, CoordY        
        integer                                     :: Ncells2D, Ncells3D, icount
        integer                                     :: i, j, k, STAT_CALL
        
        
        !----------------------------------------------------------------------

        
        icount = (Me%ModWorkSize3D%IUB - Me%ModWorkSize3D%ILB+1) *                    &
                 (Me%ModWorkSize3D%JUB - Me%ModWorkSize3D%JLB+1)
        
        Ncells2D       = icount
        Me%ModNcells2D = Ncells2D
                            
        
        allocate(Me%ModX2D         (1:Ncells2D))
        allocate(Me%ModY2D         (1:Ncells2D))
        allocate(Me%ModPropXY2D    (1:Ncells2D))
        allocate(Me%ModNoDataXY2D  (1:Ncells2D))
        
        Me%ModNoDataXY2D(:) = .true.
        
        
        allocate(Me%ModProp2D(Me%ModSize3D%ILB:Me%ModSize3D%IUB,Me%ModSize3D%JLB:Me%ModSize3D%JUB))        
        
        Me%ModProp2D(:,:) = Me%ModNODATA

        call GetZCoordinates(Me%ObjMod_Grid_Out, CoordX, CoordY, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadModGridXYZ - ModuleDelft3D_2_Mohid - ERR20' 
        endif
        
        icount = 0
        do j = Me%ModWorkSize3D%JLB, Me%ModWorkSize3D%JUB
        do i = Me%ModWorkSize3D%ILB, Me%ModWorkSize3D%IUB
            icount = icount + 1
            Me%ModX2D(icount) = CoordX(i, j)
            Me%ModY2D(icount) = CoordY(i, j)
        enddo
        enddo        
        
       
        Ncells3D       = Ncells2D * (Me%ModWorkSize3D%KUB - Me%ModWorkSize3D%KLB + 1)
        Me%ModNcells3D = Ncells3D
                                                 
    
        allocate(Me%ModX3D         (1:Ncells3D))
        allocate(Me%ModY3D         (1:Ncells3D))
        allocate(Me%ModZ3D         (1:Ncells3D))
        allocate(Me%ModPropXY3D    (1:Ncells3D))
        allocate(Me%ModNoDataXY3D  (1:Ncells3D))
        
        Me%ModNoDataXY3D(:) = .true.
        
        allocate(Me%ModProp3D(Me%ModSize3D%ILB:Me%ModSize3D%IUB,                        &
                              Me%ModSize3D%JLB:Me%ModSize3D%JUB,                        &
                              Me%ModSize3D%KLB:Me%ModSize3D%KUB))            

        icount = 0
        
        do k = Me%ModWorkSize3D%KLB,   Me%ModWorkSize3D%KUB
        do j = Me%ModWorkSize3D%JLB,   Me%ModWorkSize3D%JUB
        do i = Me%ModWorkSize3D%ILB,   Me%ModWorkSize3D%IUB
            icount = icount + 1
            Me%ModX3D(icount) = CoordX(i, j)
            Me%ModY3D(icount) = CoordY(i, j)
        enddo
        enddo
        enddo  
          
        Me%ModProp3D(:,:,:) = Me%ModNODATA
        
        call UnGetHorizontalGrid(Me%ObjMod_Grid_Out, CoordX, STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR40' 
        endif        

        call UnGetHorizontalGrid(Me%ObjMod_Grid_Out, CoordY, STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR50' 
        endif               
    
        call GetGeometryDistances(Me%ObjMod_Geo_Out, ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR60' 
        endif         
    
        !
        !FillGridMethod = 2
        call FillMatrix3D(Me%ModWorkSize3D%ILB,                                         &
                          Me%ModWorkSize3D%IUB,                                         &
                          Me%ModWorkSize3D%JLB,                                         &
                          Me%ModWorkSize3D%JUB,                                         &
                          Me%ModWorkSize3D%KLB,                                         &
                          Me%ModWorkSize3D%KUB,                                         &
                          Me%ModWaterPoints3D,                                          &
                          ZCellCenter,                                                  &
                          FillGridMethod = 2)        
    
        icount = 0
        
        do k = Me%ModWorkSize3D%KLB, Me%ModWorkSize3D%KUB
        do j = Me%ModWorkSize3D%JLB, Me%ModWorkSize3D%JUB
        do i = Me%ModWorkSize3D%ILB, Me%ModWorkSize3D%IUB        
                icount            = icount + 1
                Me%ModZ3D(icount) = ZCellCenter(i, j, k)
        enddo
        enddo   
        enddo 
        
        call UnGetGeometry(Me%ObjMod_Geo_Out, ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadAtmGridXYZ - ModuleDelft3D_2_Mohid - ERR70' 
        endif         
        
    
    end subroutine ReadModGridXYZ    

    !---------------------------------------------------------------------------    
    

    subroutine ReadBoundGridXYZ

        !Arguments-------------------------------------------------------------
                   
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer             :: ZCellCenter 
        real, dimension(:,:),  pointer              :: CoordX, CoordY        
        integer                                     :: Ncells2D, Ncells3D, icount
        integer                                     :: i, j, k, ip, STAT_CALL
        integer                                     :: Ninstants, Nlayers
        

        !----------------------------------------------------------------------

        Ninstants = int((Me%EndTime - Me%BeginTime) / Me%DT) + 1
        
        allocate (Me%BoundTime(1:Ninstants))
        
        Me%BoundTime(:) = 0.

        call ReadBoundcells
        
        Ncells2D = Me%BoundCellsNumber
        
        allocate(Me%BoundX2D         (1:Ncells2D))
        allocate(Me%BoundY2D         (1:Ncells2D))
        allocate(Me%BoundPropXY2D    (1:Ncells2D))
        allocate(Me%BoundNoDataXY2D  (1:Ncells2D))
        
      
        !point XY, time
        allocate(Me%BoundSSH_XYT2D      (1:Ncells2D,1:Ninstants))
        allocate(Me%BoundSSH_ASTRO_XYT2D(1:Ncells2D,1:Ninstants))
        
        NLayers  =  (Me%ModWorkSize3D%KUB - Me%ModWorkSize3D%KLB + 1)
        
        Ncells3D = Ncells2D * NLayers
        
        Me%ModNLayers   = NLayers
        Me%ModNinstants = Ninstants
        
        
        allocate(Me%BoundX3D         (1:Ncells3D))
        allocate(Me%BoundY3D         (1:Ncells3D))
        allocate(Me%BoundZ3D         (1:Ncells3D))        
        allocate(Me%BoundPropXYZ3D   (1:Ncells3D))
        allocate(Me%BoundNoDataXYZ3D (1:Ncells3D))        
        
       
        !point XY, layer, time
        allocate(Me%BoundSalXYZT3D             (1:Ncells2D,1:NLayers,1:Ninstants))
        allocate(Me%BoundTempXYZT3D            (1:Ncells2D,1:NLayers,1:Ninstants))
        allocate(Me%BoundVelx3D_XYZT3D         (1:Ncells2D,1:NLayers,1:Ninstants))
        allocate(Me%BoundVely3D_XYZT3D         (1:Ncells2D,1:NLayers,1:Ninstants))
        allocate(Me%BoundVelx3D_Astro_XYZT3D   (1:Ncells2D,1:NLayers,1:Ninstants))
        allocate(Me%BoundVely3D_Astro_XYZT3D   (1:Ncells2D,1:NLayers,1:Ninstants))         
        allocate(Me%BoundRiemannXYZT3          (1:Ncells2D,1:NLayers,1:Ninstants))
        
        
        call GetZCoordinates(Me%ObjMod_Grid_Out, CoordX, CoordY, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadBoundGridXYZ - ModuleDelft3D_2_Mohid - ERR20' 
        endif                

        do ip=1, Ncells2D        
            i = Me%BoundCellsI(ip)
            j = Me%BoundCellsJ(ip)
            Me%BoundX2D(ip) = CoordX(i, j)
            Me%BoundY2D(ip) = CoordY(i, j)
        enddo        
        
        Me%BoundNoDataXY2D(:) = .true.
        
        icount = 0
        do ip=1, Ncells2D        
            i = Me%BoundCellsI(ip)
            j = Me%BoundCellsJ(ip)
            do k=1, Nlayers                        
                icount = icount + 1
                Me%BoundX3D(icount) = CoordX(i, j)
                Me%BoundY3D(icount) = CoordY(i, j)
            enddo                
        enddo    
        
        Me%BoundNoDataXYZ3D(:) = .true.            
        
        
        call UnGetHorizontalGrid(Me%ObjMod_Grid_Out, CoordX, STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadBoundGridXYZ - ModuleDelft3D_2_Mohid - ERR40' 
        endif        

        call UnGetHorizontalGrid(Me%ObjMod_Grid_Out, CoordY, STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadBoundGridXYZ - ModuleDelft3D_2_Mohid - ERR50' 
        endif               
    
        call GetGeometryDistances(Me%ObjMod_Geo_Out, ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadBoundGridXYZ - ModuleDelft3D_2_Mohid - ERR60' 
        endif         

        
        icount = 0
        do ip=1, Ncells2D        
            i = Me%BoundCellsI(ip)
            j = Me%BoundCellsJ(ip)
            do k= Me%ModWorkSize3D%KLB, Me%ModWorkSize3D%KUB
                icount = icount + 1
                Me%BoundZ3D(icount) = ZCellCenter(i, j, k)
            enddo                
        enddo           
    
        call UnGetGeometry(Me%ObjMod_Geo_Out, ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadBoundGridXYZ - ModuleDelft3D_2_Mohid - ERR70' 
        endif         
        
    
    end subroutine ReadBoundGridXYZ    

    !---------------------------------------------------------------------------    
    
        
    
    subroutine ReadBoundCells

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: CharRead
        integer                                     :: STAT_CALL, nlines, i
        logical                                     :: StopRun
        
        !----------------------------------------------------------------------

        call UnitsManager(Me%BoundCellsUnit, OPEN_FILE, STAT = STAT_CALL) 
                                   
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadBoundCells - ModuleDelft3D_2_Mohid - ERR10'
        endif             
        
        open(UNIT   = Me%BoundCellsUnit,                                                &
             FILE   = trim(adjustl(Me%BoundCellsFile)),                                 &
             FORM   = "FORMATTED",                                                      &   
             IOSTAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadBoundCells - ModuleDelft3D_2_Mohid - ERR20'
        endif      
        
        nlines = 0
        do  
            read(Me%BoundCellsUnit,*,iostat=STAT_CALL)
            if (STAT_CALL/=0) exit
            nlines = nlines + 1
        enddo 

        rewind(Me%BoundCellsUnit)
        
        Me%BoundCellsNumber    = nlines * 2
        Me%BoundSectionsNumber = nlines
        
        allocate(Me%BoundSectionsName(1:Me%BoundSectionsNumber))
        allocate(Me%BoundCellsI(1:Me%BoundCellsNumber))
        allocate(Me%BoundCellsJ(1:Me%BoundCellsNumber))        

        do i=1, nlines
            read(Me%BoundCellsUnit,'(A81)') CharRead
            read(CharRead(26:49),*) Me%BoundCellsJ(2*i-1), Me%BoundCellsI(2*i-1), Me%BoundCellsJ(2*i), Me%BoundCellsI(2*i)
            Me%BoundSectionsName(i) =   CharRead(1:15)          
        enddo
        
        StopRun = .false. 
        
        do i=1, Me%BoundCellsNumber
            if (Me%BoundCellsJ(i) > Me%ModWorkSize3D%JUB .or. Me%BoundCellsJ(i) < Me%ModWorkSize3D%JLB) then
                write(*,*) 'Bound cell number =', i, 'is not define correctly. Not valid column =',Me%BoundCellsJ(i)
                StopRun = .true.
            endif
            if (Me%BoundCellsI(i) > Me%ModWorkSize3D%IUB .or. Me%BoundCellsI(i) < Me%ModWorkSize3D%ILB) then
                write(*,*) 'Bound cell number =', i, 'is not define correctly. Not valid line =',Me%BoundCellsI(i)
                StopRun = .true.
            endif            
        enddo
        
        if (StopRun) then
            stop 'ReadBoundCells - ModuleDelft3D_2_Mohid - ERR30'
        endif            
        

        call UnitsManager(Me%BoundCellsUnit, CLOSE_FILE, STAT = STAT_CALL) 
                                   
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadBoundCells - ModuleDelft3D_2_Mohid - ERR40'
        endif          
    
    end subroutine ReadBoundCells    
    

    !---------------------------------------------------------------------------
    
        
    
    subroutine WriteBoundRiemann

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: CharOut         
        integer                                     :: Ninst, nlayers, FileUnit
        integer                                     :: STAT_CALL, i
        
        !----------------------------------------------------------------------

        call UnitsManager(FileUnit, OPEN_FILE, STAT = STAT_CALL) 
                                   
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBoundRiemann - ModuleDelft3D_2_Mohid - ERR10'
        endif             
        
        open(UNIT   = FileUnit,                                                         &
             FILE   = trim(adjustl(Me%BoundRiemannFile)),                               &
             FORM   = "FORMATTED",                                                      &   
             IOSTAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBoundRiemann - ModuleDelft3D_2_Mohid - ERR20'
        endif   
        
        Ninst   = Me%ModNinstants   
        nlayers = Me%ModNlayers
                
        write(CharOut,'(I7)') nlayers * 2  
        
        do i = 1, Me%BoundSectionsNumber

            call WriteBoundTable(BoundXYZT3     = Me%BoundRiemannXYZT3,                 & 
                                 PropOutName    = "Riemann         (R)  End",           &
                                 UnitOutName    = "[m/s]",                              &
                                 FileUnit       = FileUnit,                             &
                                 Section        = i,                                    &
                                 layer          = "layer:")      
                          
        enddo 


        call UnitsManager(FileUnit, CLOSE_FILE, STAT = STAT_CALL) 
                                   
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBoundRiemann - ModuleDelft3D_2_Mohid - ERR300'
        endif          
    
    end subroutine WriteBoundRiemann    
    

    !---------------------------------------------------------------------------
    
    !---------------------------------------------------------------------------
    
        
    
    subroutine WriteBoundDensity

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: CharOut         
        integer                                     :: Ninst, nlayers, FileUnit
        integer                                     :: STAT_CALL, i
        
        !----------------------------------------------------------------------

        call UnitsManager(FileUnit, OPEN_FILE, STAT = STAT_CALL) 
                                   
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBoundDensity - ModuleDelft3D_2_Mohid - ERR10'
        endif             
        
        open(UNIT   = FileUnit,                                                         &
             FILE   = trim(adjustl(Me%BoundDensityFile)),                               &
             FORM   = "FORMATTED",                                                      &   
             IOSTAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBoundDensity - ModuleDelft3D_2_Mohid - ERR20'
        endif   
        
        Ninst   = Me%ModNinstants   
        nlayers = Me%ModNlayers
                
        write(CharOut,'(I7)') nlayers * 2  
        
        do i = 1, Me%BoundSectionsNumber

            call WriteBoundTable(BoundXYZT3     = Me%BoundSalXYZT3D,                    & 
                                 PropOutName    = "Salinity             end",           &
                                 UnitOutName    = "[ppt]",                              &
                                 FileUnit       = FileUnit,                             &
                                 Section        = i,                                    &
                                 layer          = "layer")      
                          
            call WriteBoundTable(BoundXYZT3     = Me%BoundTempXYZT3D,                   & 
                                 PropOutName    = "Temperature          end",           &
                                 UnitOutName    = "[C]",                                &
                                 FileUnit       = FileUnit,                             &
                                 Section        = i,                                    &
                                 layer          = "layer")      

        enddo 


        call UnitsManager(FileUnit, CLOSE_FILE, STAT = STAT_CALL) 
                                   
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBoundDensity - ModuleDelft3D_2_Mohid - ERR300'
        endif          
    
    end subroutine WriteBoundDensity    
    

    !---------------------------------------------------------------------------    
    
    subroutine WriteBoundTable(BoundXYZT3, PropOutName, UnitOutName, FileUnit, Section, layer)
    
        !Arguments-------------------------------------------------------------
        real,   dimension(:,:,:), pointer           :: BoundXYZT3
        character(len=*)                            :: PropOutName, UnitOutName, layer     
        integer                                     :: FileUnit, Section
        
        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: CharOut, FormatBound         
        integer                                     :: Ninst, nlayers
        integer                                     :: STAT_CALL, i, in, k
        
        !---------------------------------------------------------------------    

        Ninst   = Me%ModNinstants   
        nlayers = Me%ModNlayers  

        i       = Section
        
        write(CharOut,'(I6)') i
        write(FileUnit,'(A)')      "table-name          'Boundary Section : "//trim(adjustl(CharOut))//"'"
        write(FileUnit,'(A)')      "contents            '3d-profile'"
        write(FileUnit,'(A)')      "location            '"//trim(adjustl(Me%BoundSectionsName(i)))//"'"
        write(FileUnit,'(A)')      "time-function       'non-equidistant'"
        write(FileUnit,'(A)')      "reference-time       "//trim(adjustl(ConvertDateToString(Me%TimeRef)))
        write(FileUnit,'(A)')      "time-unit           'minutes'"
        write(FileUnit,'(A)')      "interpolation       'linear'"
        write(FileUnit,'(A)')      "parameter           'time' unit '[min]'"        

        write(CharOut,'(I7)') nlayers * 2  
        FormatBound = '(f15.2,'//trim(adjustl(CharOut))//trim(adjustl(Me%ModBoundFormat))//')'
       
   
        do k=1, nlayers 
            write(CharOut,'(I3)') k
            write(FileUnit,'(A)') "parameter           '"//trim(adjustl(PropOutName))//" A "//trim(adjustl(layer))//" "//&
                                  trim(adjustl(CharOut))//"' unit '"//trim(adjustl(UnitOutName))//"'"  
        enddo
        
        do k=1, nlayers 
            write(CharOut,'(I3)') k
            write(FileUnit,'(A)') "parameter           '"//trim(adjustl(PropOutName))//" B "//trim(adjustl(layer))//" "//&
                                  trim(adjustl(CharOut))//"' unit '"//trim(adjustl(UnitOutName))//"'"  
        enddo            
        
        write(FileUnit,*) 'records-in-table', Ninst            
        do in = 1, Ninst
            write(FileUnit,FormatBound,iostat=STAT_CALL) Me%BoundTime(in),          &
                  (BoundXYZT3(i*2-1,   k, in),k=1,nlayers),               &
                  (BoundXYZT3(i*2,     k, in),k=1,nlayers)
        enddo            
        
    end subroutine WriteBoundTable        
        
    !---------------------------------------------------------------------------        
    
    subroutine OpenModOutFile_Ini

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------

        call UnitsManager(Me%ModIniOutUnit, OPEN_FILE, STAT = STAT_CALL) 
                                   
        if (STAT_CALL /= SUCCESS_) then
            stop 'OpenModOutFile_Ini - ModuleDelft3D_2_Mohid - ERR10'
        endif             
        
        open(UNIT   = Me%ModIniOutUnit,                                                 &
             FILE   = trim(adjustl(Me%ModOutFile)),                                     &
             FORM   = "FORMATTED",                                                      &   
             IOSTAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop 'OpenModOutFile_Ini - ModuleDelft3D_2_Mohid - ERR20'
        endif      
        
    
    end subroutine OpenModOutFile_Ini    
    

    !---------------------------------------------------------------------------
    
    subroutine ReadProp2D(CurrentTime, i)

        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: i                        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, icount, ii, jj
        
        !----------------------------------------------------------------------
        
        Me%AtmNoDataXY(:) = .true.
        Me%AtmPropXY  (:) = FillValueReal

        call ModifyField4DXYZ(Field4DID             = Me%ObjField4D_Atm(i),             &
                              PropertyIDNumber      = Me%AtmProp(i)%IDNumber,           &
                              CurrentTime           = CurrentTime,                      &
                              X                     = Me%AtmX,                          &    
                              Y                     = Me%AtmY,                          &
                              Field                 = Me%AtmPropXY,                     &
                              NoData                = Me%AtmNoDataXY,                   &
                              STAT                  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadProp2D - ModuleDelft3D_2_Mohid - ERR10'
        endif
        
        icount = 0
        
        do jj = Me%AtmWorkSize2D%JLB, Me%AtmWorkSize2D%JUB
        do ii = Me%AtmWorkSize2D%ILB, Me%AtmWorkSize2D%IUB        
            
            icount           = icount + 1
            if (Me%AtmNoDataXY(icount)) then
                Me%AtmProp2D(ii, jj) = Me%AtmNODATA
            else                        
                Me%AtmProp2D(ii, jj) = Me%AtmPropXY(icount)
            endif
            
        enddo   
        enddo      
    
    end subroutine ReadProp2D
    
    !---------------------------------------------------------------------------
    
    !---------------------------------------------------------------------------
    
    subroutine ReadMod2D(CurrentTime, ObjField4D, PropIDNumber, ValueIni)

        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: ObjField4D
        integer                                     :: PropIDNumber     
        real, optional                              :: ValueIni

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, icount, ii, jj
        
        !----------------------------------------------------------------------
    
        Me%ModNoDataXY2D(:) = .true.
        Me%ModPropXY2D  (:) = FillValueReal

        call ModifyField4DXYZ(Field4DID             = ObjField4D,                       &
                              PropertyIDNumber      = PropIDNumber,                     &
                              CurrentTime           = CurrentTime,                      &
                              X                     = Me%ModX2D,                        &    
                              Y                     = Me%ModY2D,                        &
                              Field                 = Me%ModPropXY2D,                   &
                              NoData                = Me%ModNoDataXY2D,                 &
                              STAT                  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadMod2D - ModuleDelft3D_2_Mohid - ERR10'
        endif
        
        icount = 0
        
        do jj = Me%ModWorkSize3D%JLB, Me%ModWorkSize3D%JUB
        do ii = Me%ModWorkSize3D%ILB, Me%ModWorkSize3D%IUB        
            icount = icount + 1
            if (Me%ModNoDataXY2D(icount)) then
                !Me%ModProp2D(ii, jj) = Me%ModNODATA
                write(*,*) 'i,j=', ii, jj
                stop 'ReadMod2D - ModuleDelft3D_2_Mohid - ERR20' 
            else                        
                Me%ModProp2D(ii, jj) = Me%ModPropXY2D(icount)
            endif
        enddo   
        enddo      
    
        if (present(ValueIni)) then    
            if (CurrentTime==Me%BeginTime) then
                call WriteModIn2D(ValueIni = ValueIni)
            endif        
        endif    
    end subroutine ReadMod2D
    
    !---------------------------------------------------------------------------    
    
    !---------------------------------------------------------------------------
    
    subroutine ReadMod3D(CurrentTime, ObjField4D, PropIDNumber, ValueIni)

        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: ObjField4D
        integer                                     :: PropIDNumber
        real, optional                              :: ValueIni
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, icount, ii, jj, k
        
        !----------------------------------------------------------------------
        
    
        Me%ModNoDataXY3D(:) = .true.
        Me%ModPropXY3D  (:) = FillValueReal

        call ModifyField4DXYZ(Field4DID             = ObjField4D,                       &
                              PropertyIDNumber      = PropIDNumber,                     &
                              CurrentTime           = CurrentTime,                      &
                              X                     = Me%ModX3D,                        &    
                              Y                     = Me%ModY3D,                        &
                              Z                     = Me%ModZ3D,                        &                              
                              Field                 = Me%ModPropXY3D,                   &
                              NoData                = Me%ModNoDataXY3D,                 &
                              STAT                  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadMod3D - ModuleDelft3D_2_Mohid - ERR10'
        endif
        
        icount = 0

        do  k = Me%ModWorkSize3D%KLB, Me%ModWorkSize3D%KUB        
        do jj = Me%ModWorkSize3D%JLB, Me%ModWorkSize3D%JUB
        do ii = Me%ModWorkSize3D%ILB, Me%ModWorkSize3D%IUB        
            icount = icount + 1                   
            if (Me%ModNoDataXY3D(icount)) then
                !Me%ModProp3D(ii, jj,k) = Me%ModNODATA
                write(*,*) 'i,j,k=', ii, jj, k
                stop 'ReadMod3D - ModuleDelft3D_2_Mohid - ERR20' 
            else                        
                Me%ModProp3D(ii, jj, k) = Me%ModPropXY3D(icount)
            endif
        enddo   
        enddo      
        enddo
        
    
        if (present(ValueIni)) then    
            if (CurrentTime==Me%BeginTime) then
                call WriteModIn3D(ValueIni = ValueIni)
            endif        
        endif      
        
        
    end subroutine ReadMod3D
    
    !---------------------------------------------------------------------------  
    
    subroutine ReadBound2D(CurrentTime, Instant, ObjField4D, PropIDNumber, OutProp2D)

        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: Instant
        integer                                     :: ObjField4D
        integer                                     :: PropIDNumber     
        real,       dimension(:,:), pointer         :: OutProp2D

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ip
        
        !----------------------------------------------------------------------
    
        Me%BoundNoDataXY2D(:) = .true.
        Me%BoundPropXY2D  (:) = FillValueReal

        call ModifyField4DXYZ(Field4DID             = ObjField4D,                       &
                              PropertyIDNumber      = PropIDNumber,                     &
                              CurrentTime           = CurrentTime,                      &
                              X                     = Me%BoundX2D,                      &    
                              Y                     = Me%BoundY2D,                      &
                              Field                 = Me%BoundPropXY2D,                 &
                              NoData                = Me%BoundNoDataXY2D,               &
                              STAT                  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadBound2D - ModuleDelft3D_2_Mohid - ERR10'
        endif
        
        do ip= 1, Me%BoundCellsNumber
            if (Me%BoundNoDataXY2D(ip)) then
                write(*,*) "Boundary cell =", ip
                write(*,*) "Instant       =", Instant              
                stop 'ReadBound2D - ModuleDelft3D_2_Mohid - ERR20'
            else
                OutProp2D(ip, Instant) = Me%BoundPropXY2D(iP)
            endif
        enddo

    end subroutine ReadBound2D
    
    !---------------------------------------------------------------------------  
    
    subroutine ReadBound3D(CurrentTime, Instant, ObjField4D, PropIDNumber, OutProp3D)

        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: Instant
        integer                                     :: ObjField4D
        integer                                     :: PropIDNumber     
        real,       dimension(:,:,:), pointer       :: OutProp3D

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ip, k, icount
        
        !----------------------------------------------------------------------
    
        Me%BoundNoDataXYZ3D(:) = .true.
        Me%BoundPropXYZ3D  (:) = FillValueReal        

        call ModifyField4DXYZ(Field4DID             = ObjField4D,                       &
                              PropertyIDNumber      = PropIDNumber,                     &
                              CurrentTime           = CurrentTime,                      &
                              X                     = Me%BoundX3D,                      &    
                              Y                     = Me%BoundY3D,                      &
                              Z                     = Me%BoundZ3D,                      &
                              Field                 = Me%BoundPropXYZ3D,                &
                              NoData                = Me%BoundNoDataXYZ3D,              &
                              STAT                  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadBound3D - ModuleDelft3D_2_Mohid - ERR10'
        endif
        
        icount = 0
        do ip= 1, Me%BoundCellsNumber
            do k=1, Me%ModNlayers
                icount = icount + 1
                if (Me%BoundNoDataXYZ3D(icount)) then
                    write(*,*) "Boundary cell =", ip
                    write(*,*) "Layer         =", k
                    write(*,*) "Instant       =", Instant                    
                    stop 'ReadBound3D - ModuleDelft3D_2_Mohid - ERR20'
                else
                    OutProp3D(ip, k, Instant) = Me%BoundPropXYZ3D(icount)
                endif
            enddo                                
        enddo

    end subroutine ReadBound3D
    
    !---------------------------------------------------------------------------

    subroutine ComputeRiemann(Instant)

        !Arguments-------------------------------------------------------------
        integer                                     :: Instant

        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer             :: Bathymetry
        real                                        :: CH, Riemann
        integer                                     :: STAT_CALL, ip, i, j, k, icount
        
        !----------------------------------------------------------------------
    
        !Bathymetry
        call GetGridData( Me%ObjMod_Bathym_Out, Bathymetry, STAT_CALL)     
        if (STAT_CALL  /= SUCCESS_) stop 'ComputeRiemann - ModuleDelft3D_2_Mohid - ERR10'

    
        
        icount = 0
        do ip= 1, Me%BoundCellsNumber
            i        = Me%BoundCellsI(ip)
            j        = Me%BoundCellsJ(ip)
            if (Bathymetry(i, j) > 1.) then
                !s-1
                CH = sqrt(Gravity * Bathymetry(i, j)) / Bathymetry(i, j)
            else
                CH = sqrt(Gravity * 1.) / 1.
            endif                
            do k=1, Me%ModNlayers
                !Initialization
                Riemann = 0.

                Riemann = Riemann + Me%BoundVely3D_XYZT3D(ip,k,instant)
                    
                if (Me%BoundAstro) then
                    Riemann = Riemann + Me%BoundVely3D_Astro_XYZT3D(ip,k,instant)                    
                endif                        
                
                !south boundary
                if      (i == Me%ModWorkSize3d%ILB) then

                    !m/s                                   s-1 * m 
                    Riemann = Riemann + CH * (Me%BoundSSH_XYT2D(ip, Instant))
                
                    if (Me%BoundAstro) then
                        Riemann = Riemann + CH * (Me%BoundSSH_ASTRO_XYT2D(ip, Instant))
                    endif                                        
                
                
                !north boundary
                elseif  (i == Me%ModWorkSize3d%IUB) then

                    !m/s                                   s-1 * m 
                    Riemann = Riemann - CH * (Me%BoundSSH_XYT2D(ip, Instant))
                
                    if (Me%BoundAstro) then
                        Riemann = Riemann - CH * (Me%BoundSSH_ASTRO_XYT2D(ip, Instant))
                    endif                                        
                
                !west boundary
                elseif  (j == Me%ModWorkSize3d%JLB) then

                    !m/s                                   s-1 * m 
                    Riemann = Riemann + CH * (Me%BoundSSH_XYT2D(ip, Instant))
                
                    if (Me%BoundAstro) then
                        Riemann = Riemann + CH * (Me%BoundSSH_ASTRO_XYT2D(ip, Instant))
                    endif                                        
                

                !east boundary
                elseif  (j == Me%ModWorkSize3d%JUB) then
                    !m/s                                   s-1 * m 
                    Riemann = Riemann - CH * (Me%BoundSSH_XYT2D(ip, Instant))
                
                    if (Me%BoundAstro) then
                        Riemann = Riemann - CH * (Me%BoundSSH_ASTRO_XYT2D(ip, Instant))
                    endif                  

                endif
                
                Me%BoundRiemannXYZT3(ip, k, Instant) = Riemann
                
            enddo                                
        enddo

        call UnGetGridData( Me%ObjMod_Bathym_Out, Bathymetry, STAT_CALL)     
        if (STAT_CALL  /= SUCCESS_) stop 'ComputeRiemann - ModuleDelft3D_2_Mohid - ERR20'


    end subroutine ComputeRiemann
    


    !---------------------------------------------------------------------------      
    subroutine WriteProp2D(CurrentTime, i)

        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: i
        !Local-----------------------------------------------------------------
        character(len=6)                            :: a1
        character(len=11)                           :: a2
        character(len=39)                           :: a3
        character(len=56)                           :: at
        integer                                     :: Hours, il, jl
        
        !----------------------------------------------------------------------
        
        
        Hours = int((CurrentTime - Me%TimeRef) / 3600.)
        
        a1 = "TIME ="
        write(a2,'(I11)') Hours
        a3 = " hours since 2018-01-01 00:00:00 +00:00"        
        at = a1//a2//a3
        write(Me%AtmOutUnit,'(A56)') at
        do il = Me%AtmWorkSize2D%ILB, Me%AtmWorkSize2D%IUB        
            write(Me%AtmOutUnit,trim(Me%AtmOutFormat(i))) (Me%AtmProp2D(il, jl),jl = Me%AtmWorkSize2D%JLB, Me%AtmWorkSize2D%JUB) 
        enddo   
    
    end subroutine WriteProp2D    
  
      !---------------------------------------------------------------------------  
    
    subroutine WriteModIn2D(ValueIni)

        !Arguments-------------------------------------------------------------
        real                                :: ValueIni
        !Local-----------------------------------------------------------------
        integer                             :: i 
        !----------------------------------------------------------------------

        if (ValueIni > FillValueReal) then        
            Me%ModProp2D(:,:) = ValueIni
        endif
          
        do i = Me%ModWorkSize3D%ILB-1, Me%ModWorkSize3D%IUB+1        
            write(Me%ModIniOutUnit,trim(Me%ModIniFormat)) Me%ModProp2D(i, Me%ModWorkSize3D%JLB-1:Me%ModWorkSize3D%JUB+1) 
        enddo   
  
    end subroutine WriteModIn2D    
  
    !---------------------------------------------------------------------------          
        
        
    subroutine WriteModIn3D(ValueIni)

        !Arguments-------------------------------------------------------------
        real                                :: ValueIni
        !Local-----------------------------------------------------------------
        integer                             :: i, k 
        !----------------------------------------------------------------------

        if (ValueIni > FillValueReal) then        
            Me%ModProp3D(:,:,:) = ValueIni
        endif

        do k = Me%ModWorkSize3D%KLB, Me%ModWorkSize3D%KUB                  
        do i = Me%ModWorkSize3D%ILB-1, Me%ModWorkSize3D%IUB+1        
            write(Me%ModIniOutUnit,trim(Me%ModIniFormat)) Me%ModProp3D(i, Me%ModWorkSize3D%JLB-1:Me%ModWorkSize3D%JUB+1, k) 
        enddo 
        enddo          

    end subroutine WriteModIn3D    
  
  !---------------------------------------------------------------------------    
   
  
    subroutine WriteModInOut


        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------
        

        if (Me%Ocean_IN_Ini) then

            call ModIni(CurrentTime = Me%BeginTime )    
            
        endif            
        
        if (Me%Ocean_IN_Bound) then
            call ModBound           
        endif            
        
     
        
        call KillField4D(Field4DID = Me%ObjField4D_SSH,                                 &
                         STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteModInOut - ModuleDelft3D_2_Mohid - ERR10'
        endif              

        call KillField4D(Field4DID = Me%ObjField4D_VelX,                                &
                         STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteModInOut - ModuleDelft3D_2_Mohid - ERR20'
        endif                

        call KillField4D(Field4DID = Me%ObjField4D_VelY,                                &
                         STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteModInOut - ModuleDelft3D_2_Mohid - ERR30'
        endif                

        call KillField4D(Field4DID = Me%ObjField4D_Sal,                                &
                         STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteModInOut - ModuleDelft3D_2_Mohid - ERR40'
        endif                
        
        call KillField4D(Field4DID = Me%ObjField4D_Temp,                                &
                         STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteModInOut - ModuleDelft3D_2_Mohid - ERR50'
        endif                        
        
        
        if (Me%BoundAstro) then
            call KillField4D(Field4DID = Me%ObjField4D_SSH_Astro,                       &
                             STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteModInOut - ModuleDelft3D_2_Mohid - ERR60'
            endif              

            call KillField4D(Field4DID = Me%ObjField4D_VelX_Astro,                      &
                             STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteModInOut - ModuleDelft3D_2_Mohid - ERR70'
            endif                

            call KillField4D(Field4DID = Me%ObjField4D_VelY_Astro,                      &
                             STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteModInOut - ModuleDelft3D_2_Mohid - ERR80'
            endif                      
        endif
        
    end subroutine WriteModInOut
    
  !---------------------------------------------------------------------------  

    subroutine ModBound    


        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: it
        type (T_Time)                               :: CurrentTime

        !----------------------------------------------------------------------
        
        CurrentTime = Me%BeginTime       
        
        call ReadBoundGridXYZ                 
        
        it = 0
        
        do while (CurrentTime<=Me%EndTime)
        
             it = it + 1
        
            !SSH        
            call ReadBound2D(   CurrentTime     = CurrentTime,                          &
                                Instant         = it,                                   &
                                ObjField4D      = Me%ObjField4D_SSH,                    &     
                                PropIDNumber    = Me%SSHProp%IDNumber,                  &    
                                OutProp2D       = Me%BoundSSH_XYT2D)    
                                
            !velocity X
            call ReadBound3D(   CurrentTime     = CurrentTime,                          &
                                Instant         = it,                                   &
                                ObjField4D      = Me%ObjField4D_VelX,                   &     
                                PropIDNumber    = Me%VelXProp%IDNumber,                 &
                                OutProp3D       = Me%BoundVelx3D_XYZT3D)               

                            
            !velocity Y
            call ReadBound3D(   CurrentTime     = CurrentTime,                          &
                                Instant         = it,                                   &
                                ObjField4D      = Me%ObjField4D_VelY,                   &     
                                PropIDNumber    = Me%VelYProp%IDNumber,                 &
                                OutProp3D       = Me%BoundVely3D_XYZT3D)                           

            if (Me%BoundAstro) then
                !SSH Astro       
                call ReadBound2D(   CurrentTime     = CurrentTime,                      &
                                    Instant         = it,                               &
                                    ObjField4D      = Me%ObjField4D_SSH_Astro,          &     
                                    PropIDNumber    = Me%SSHProp%IDNumber,              &
                                    OutProp2D       = Me%BoundSSH_ASTRO_XYT2D)   
                                    

                !velocity X astro
                call ReadBound3D(   CurrentTime     = CurrentTime,                      &
                                    Instant         = it,                               &
                                    ObjField4D      = Me%ObjField4D_VelX_Astro,         &     
                                    PropIDNumber    = Me%VelXProp%IDNumber,             &
                                    OutProp3D       = Me%BoundVelx3D_Astro_XYZT3D)           

                                
                !velocity Y astro
                call ReadBound3D(   CurrentTime     = CurrentTime,                      &
                                    Instant         = it,                               &
                                    ObjField4D      = Me%ObjField4D_VelY_Astro,         &     
                                    PropIDNumber    = Me%VelYProp%IDNumber,             &
                                    OutProp3D       = Me%BoundVely3D_Astro_XYZT3D)
            endif

                                
            call ComputeRiemann(it)
            
            !salinity
            call ReadBound3D(   CurrentTime     = CurrentTime,                          &
                                Instant         = it,                                   &
                                ObjField4D      = Me%ObjField4D_Sal,                    &     
                                PropIDNumber    = Me%SalProp%IDNumber,                  &
                                OutProp3D       = Me%BoundSalXYZT3D)   

            !temperature
            call ReadBound3D(   CurrentTime     = CurrentTime,                          &
                                Instant         = it,                                   &
                                ObjField4D      = Me%ObjField4D_Temp,                   &     
                                PropIDNumber    = Me%TempProp%IDNumber,                 &
                                OutProp3D       = Me%BoundTempXYZT3D)  

            !Minutes
            Me%BoundTime(it) = (CurrentTime - Me%TimeRef) / 60.


            CurrentTime = CurrentTime + Me%DT
            
            
        
        enddo
        


        call WriteBoundRiemann
        
        call WriteBoundDensity
                                   
      


    end subroutine ModBound
    
    !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------  

    subroutine ReadHydro    


        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: it
        type (T_Time)                               :: CurrentTime
        real, dimension(:,:),  pointer              :: CoordX, CoordY 
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------
        
        if (Me%ObjMod_Grid_Out == 0) call ConstructModGrid       
        
        call ConstructModProp( block_begin  = '<<begin_ssh>>',                          &
                               block_end    = '<<end_ssh>>',                            &
                               ObjField4D   = Me%ObjField4D_Ssh,                        & 
                               PropID       = Me%SSHProp,                               &
                               ValueIni     = Me%SSHValueIni)

        call ConstructModProp( block_begin  = '<<begin_velx>>',                         &
                               block_end    = '<<end_velx>>',                           &
                               ObjField4D   = Me%ObjField4D_VelX,                       & 
                               PropID       = Me%VelXProp,                              &
                               ValueIni     = Me%VelxValueIni)
                              
        call ConstructModProp( block_begin  = '<<begin_vely>>',                         &
                               block_end    = '<<end_vely>>',                           &
                               ObjField4D   = Me%ObjField4D_VelY,                       & 
                               PropID       = Me%VelYProp,                              &
                               ValueIni     = Me%VelyValueIni)            
        
        CurrentTime = Me%BeginTime       
        
        call GetZCoordinates(Me%ObjMod_Grid_Out, CoordX, CoordY, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadHydro - ModuleDelft3D_2_Mohid - ERR10' 
        endif
        
  
        
        allocate(Me%BoundX2D            (1))
        allocate(Me%BoundY2D            (1))       
        allocate(Me%BoundPropXY2D       (1)) 
        allocate(Me%BoundNoDataXY2D     (1))        
         
        allocate(Me%BoundX3D            (1))
        allocate(Me%BoundY3D            (1))
        allocate(Me%BoundZ3D            (1))
        allocate(Me%BoundPropXYZ3D      (1))
        allocate(Me%BoundNoDataXYZ3D    (1))    

        Me%BoundX2D         (1) = CoordX(1, 1)
        Me%BoundY2D         (1) = CoordY(1, 1)
         
        Me%BoundX3D         (1) = CoordX(1, 1)       
        Me%BoundY3D         (1) = CoordY(1, 1)       
        Me%BoundZ3D         (1) = 0.
        
        it = 0
        
        do while (CurrentTime<=Me%EndTime)
        
             it = it + 1
        
            !SSH        
            call ReadBound2D(   CurrentTime     = CurrentTime,                          &
                                Instant         = it,                                   &
                                ObjField4D      = Me%ObjField4D_SSH,                    &     
                                PropIDNumber    = Me%SSHProp%IDNumber,                  &    
                                OutProp2D       = Me%BoundSSH_XYT2D)    
                                
            !velocity X
            call ReadBound3D(   CurrentTime     = CurrentTime,                          &
                                Instant         = it,                                   &
                                ObjField4D      = Me%ObjField4D_VelX,                   &     
                                PropIDNumber    = Me%VelXProp%IDNumber,                 &
                                OutProp3D       = Me%BoundVelx3D_XYZT3D)               

                            
            !velocity Y
            call ReadBound3D(   CurrentTime     = CurrentTime,                          &
                                Instant         = it,                                   &
                                ObjField4D      = Me%ObjField4D_VelY,                   &     
                                PropIDNumber    = Me%VelYProp%IDNumber,                 &
                                OutProp3D       = Me%BoundVely3D_XYZT3D)                           


            CurrentTime = CurrentTime + Me%DT
            
            
        
        enddo
        
        deallocate(Me%BoundX2D        )
        deallocate(Me%BoundY2D        )       
        deallocate(Me%BoundPropXY2D   ) 
        deallocate(Me%BoundNoDataXY2D )        
         
        deallocate(Me%BoundX3D        )
        deallocate(Me%BoundY3D        )
        deallocate(Me%BoundZ3D        )
        deallocate(Me%BoundPropXYZ3D  )
        deallocate(Me%BoundNoDataXYZ3D) 

        call UnGetHorizontalGrid(Me%ObjMod_Grid_Out, CoordX, STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadHydro - ModuleDelft3D_2_Mohid - ERR20' 
        endif        

        call UnGetHorizontalGrid(Me%ObjMod_Grid_Out, CoordY, STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_) then
            stop 'ReadHydro - ModuleDelft3D_2_Mohid - ERR30' 
        endif     
        
    end subroutine ReadHydro
    
    !---------------------------------------------------------------------------
    
  
    subroutine ModIni(CurrentTime)


        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL


        !----------------------------------------------------------------------
        
        call OpenModOutFile_Ini

        call ReadModGridXYZ
        
        !velocity SSH        
        call ReadMod2D (CurrentTime     = CurrentTime,                              &
                        ObjField4D      = Me%ObjField4D_SSH,                        &
                        PropIDNumber    = Me%SSHProp%IDNumber,                      &
                        ValueIni        = Me%SSHValueIni)
            
        !velocity X
        call ReadMod3D (CurrentTime     = CurrentTime,                              &
                        ObjField4D      = Me%ObjField4D_VelX,                       &
                        PropIDNumber    = Me%VelXProp%IDNumber,                     &
                        ValueIni        = Me%VelxValueIni)
                        
        !velocity Y
        call ReadMod3D (CurrentTime     = CurrentTime,                              &
                        ObjField4D      = Me%ObjField4D_VelY,                       &
                        PropIDNumber    = Me%VelYProp%IDNumber,                     &
                        ValueIni        = Me%VelYValueIni)
        !salinity
        call ReadMod3D (CurrentTime     = CurrentTime,                              &
                        ObjField4D      = Me%ObjField4D_Sal,                        &
                        PropIDNumber    = Me%SalProp%IDNumber,                      &
                        ValueIni        = Me%SalValueIni)

        !temperature
        call ReadMod3D (CurrentTime     = CurrentTime,                              &
                        ObjField4D      = Me%ObjField4D_Temp,                       &
                        PropIDNumber    = Me%TempProp%IDNumber,                     &
                        ValueIni        = Me%TempValueIni)

        call UnitsManager(Me%ModIniOutUnit, CLOSE_FILE, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) then
            stop 'ModIni - ModuleDelft3D_2_Mohid - ERR10'
        endif           

    end subroutine ModIni

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillDelft3D_2_Mohid(ObjDelft3D_2_MohidID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjDelft3D_2_MohidID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDelft3D_2_MohidID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mDelft3D_2_Mohid_,  Me%InstanceID)

            if (nUsers == 0) then

                ObjDelft3D_2_MohidID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillDelft3D_2_Mohid
        

    !------------------------------------------------------------------------
    
    

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjDelft3D_2_Mohid_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjDelft3D_2_Mohid_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjDelft3D_2_Mohid_ID > 0) then
            call LocateObjDelft3D_2_Mohid (ObjDelft3D_2_Mohid_ID)
            ready_ = VerifyReadLock (mDelft3D_2_Mohid_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjDelft3D_2_Mohid (ObjDelft3D_2_MohidID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjDelft3D_2_MohidID

        !Local-----------------------------------------------------------------

        Me => FirstObjDelft3D_2_Mohid
        do while (associated (Me))
            if (Me%InstanceID == ObjDelft3D_2_MohidID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleDelft3D_2_Mohid - LocateObjDelft3D_2_Mohid - ERR01'

    end subroutine LocateObjDelft3D_2_Mohid

    !--------------------------------------------------------------------------

end module ModuleDelft3D_2_Mohid









