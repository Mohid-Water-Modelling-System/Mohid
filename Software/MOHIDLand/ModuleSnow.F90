!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Snow
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jul 2014
! REVISION      : Eduardo Jauch - v4.0
! DESCRIPTION   : Module which calculates the Snow Melting
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
   
module ModuleSnow

   use ModuleGlobalData
   use ModuleTime
   use ModuleTimeSerie         ,only : StartTimeSerie, StartTimeSerieInput,             &
                                       KillTimeSerie, GetNumberOfTimeSeries,            &
                                       GetTimeSerieInitialData, GetTimeSerieValue,      &
                                       GetTimeSerieLocation, GetTimeSerieName,          &
                                       TryIgnoreTimeSerie, CorrectsCellsTimeSerie
   use ModuleEnterData
   use ModuleHDF5
   use ModuleFunctions         ,only : TimeToString, SetMatrixValue, ChangeSuffix,      &
                                       CHUNK_J, LinearInterpolation,                    &
                                       InterpolateValueInTime, ConstructPropertyID
   use ModuleHorizontalGrid    ,only : GetHorizontalGridSize, GetHorizontalGrid,        &
                                       UnGetHorizontalGrid, WriteHorizontalGrid,        &
                                       GetGridCellArea, GetXYCellZ,                     &
                                       GetCellZInterceptByLine,                         &
                                       GetCellZInterceptByPolygon
   use ModuleHorizontalMap     ,only : GetBoundaries, UngetHorizontalMap
   use ModuleGridData          ,only : GetGridData, UngetGridData, WriteGridData
   use ModuleBasinGeometry     ,only : GetBasinPoints, GetRiverPoints, GetCellSlope,    &
                                       GetDrainageDirection, TargetPoint,               &
                                       UnGetBasin
   use ModuleStopWatch         ,only : StartWatch, StopWatch
   use ModuleFillMatrix        ,only : ConstructFillMatrix, ModifyFillMatrix,           &
                                       KillFillMatrix
   
   implicit none

   private 
   
   !Subroutines---------------------------------------------------------------

   !Constructor
   public  ::  ConstructSnow
   private ::      AllocateInstance
   private ::      ReadDataFile
   private ::      AllocateVariables
   private ::      InitializeVariables
   private ::      ConstructHDF5Output
   
   !Selector
   public  ::  GetSnowMelting
   public  ::  UnGetSnowMelting
   
   !Modifier
   public  ::  ModifySnow
   
   !Destructor
   public  ::  KillSnow     
   
   !Management
   private ::  ReadLockExternalVar
   private ::  ReadUnLockExternalVar
   private ::  Ready
   private ::      LocateObjSnow    
   
   !Interfaces----------------------------------------------------------------
   private :: UnGetSnowMelting2D_R4
   private :: UnGetSnowMelting2D_R8
   interface  UnGetSnowMelting
       module procedure UnGetSnowMelting2D_R4
       module procedure UnGetSnowMelting2D_R8
   end interface  UnGetSnowMelting
   
   !Parameters----------------------------------------------------------------
   !integer, parameter                              :: KinematicWave_   = 1   
   
   !Types--------------------------------------------------------------------- 
   type T_OutPut
      type (T_Time), pointer, dimension(:)       :: OutTime                   => null()
      integer                                    :: NextOutPut                = 1
      logical                                    :: Yes                       = .false.
      type (T_Time), dimension(:), pointer       :: RestartOutTime            => null()
      logical                                    :: WriteRestartFile          = .false.
      logical                                    :: RestartOverwrite          = .false.
      integer                                    :: NextRestartOutput         = 1 
      logical                                    :: BoxFluxes                 = .false.
      logical                                    :: TimeSerie_On              = .false.
      logical                                    :: HDF_On                    = .false.
      integer                                    :: Number                    = 0
   end type T_OutPut
    
   type T_Files
      character(PathLength)                       :: DataFile                 = null_str
      character(PathLength)                       :: InitialFile              = null_str
      character(PathLength)                       :: FinalFile                = null_str
      character(PathLength)                       :: TransientHDF             = null_str
      !character(PathLength)                       :: BoxesFile                = null_str
   end type T_Files     
   
   type T_Evolution
      type(T_Time)                                :: LastCompute
      type(T_Time)                                :: NextCompute      
   end type T_Evolution
   
   type T_Property
      type(T_PropertyID)                          :: ID !From ModuleGlobalData
      type(T_Evolution)                           :: Evolution
      
      logical                                     :: Old       = .false., &
                                                     TimeSerie = .false., &
                                                     OutputHDF = .false.
      
      real, dimension(:,:), pointer               :: ValueOld => null(), &
                                                     Value    => null() 
      
      type(T_Property), pointer                   :: Next     => null(), &
                                                     Prev     => null()
   end type T_Property
   
   type T_ExtVar
      integer, dimension(:,:), pointer            :: BasinPoints              => null()
      real   , dimension(:,:), pointer            :: Topography               => null()
      integer, dimension(:,:), pointer            :: RiverPoints              => null()
      !real   , dimension(:,:), pointer            :: CellSlope                => null()
      type (T_Time)                               :: Now
      real                                        :: DT                       = null_real
   end type T_ExtVar   
   
   type  T_Snow
      integer                                     :: InstanceID               = 0
      character(len=StringLength)                 :: ModelName                = null_str
      
      integer                                     :: ObjBasinGeometry         = 0
      integer                                     :: ObjTime                  = 0
      integer                                     :: ObjHorizontalGrid        = 0
      integer                                     :: ObjHorizontalMap         = 0
      integer                                     :: ObjGridData              = 0
      integer                                     :: ObjHDF5                  = 0
      integer                                     :: ObjIniHDF5               = 0
      integer                                     :: ObjEnterData             = 0  
      integer                                     :: ObjTimeSerie             = 0
      
      type (T_OutPut)                             :: OutPut
      type (T_ExtVar)                             :: ExtVar      
      type (T_Files)                              :: Files
      
      type (T_Time)                               :: BeginTime
      type (T_Time)                               :: EndTime      
      
      !Grid size
      type (T_Size2D)                             :: Size
      type (T_Size2D)                             :: WorkSize
      
      logical                                     :: Continuous               = .false.
      logical                                     :: StopOnWrongDate          = .true.      

      type(T_Snow), pointer                       :: Next                     => null()  
      
      real                                        :: SnowMeltingDT            = 86400. 
      real                                        :: MeltingTemperature       = 0.0
      
      integer                                     :: PropertiesNumber         = 0
      
      type(T_Property), pointer                   :: FirstProperty            => null(), &
                                                     LastProperty             => null(), &
                                                     SnowPack                 => null(), &
                                                     DailyAvgTemp             => null(), &
                                                     Albedo                   => null(), &
                                                     ForestCoverFraction      => null(), &
                                                     SlopeFactor              => null() 
                                                     
      real, dimension(:,:), pointer               :: SnowMeltingFlux          => null(), &
                                                     SnowMelted               => null()
      
   end type T_Snow
   
   !Global Module Variables
   type (T_Snow), pointer                         :: FirstObjSnow             => null()
   type (T_Snow), pointer                         :: Me                       => null()
 
   !--------------------------------------------------------------------------
    
   contains

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
   
   subroutine ConstructSnow(ModelName,                                       &            
                            SnowID,                                          &
                            ComputeTimeID,                                   &
                            HorizontalGridID,                                &
                            HorizontalMapID,                                 &
                            GridDataID,                                      &
                            BasinGeometryID,                                 &
                            STAT)

      !Arguments---------------------------------------------------------------
      character(len=*)                                :: ModelName
      integer                                         :: SnowID
      integer                                         :: ComputeTimeID
      integer                                         :: HorizontalGridID
      integer                                         :: HorizontalMapID
      integer                                         :: GridDataID
      integer                                         :: BasinGeometryID      
      integer, optional, intent(OUT)                  :: STAT           

      !External----------------------------------------------------------------
      integer                                         :: ready_         

      !Local-------------------------------------------------------------------
      integer                                         :: STAT_, STAT_CALL

      !------------------------------------------------------------------------
      STAT_ = UNKNOWN_

      !Assures nullification of the global variable
      if (.not. ModuleIsRegistered(mSnow_)) then
         nullify (FirstObjSnow)
         call RegisterModule (mSnow_) 
      endif

      call Ready(SnowID, ready_)    

cd0 : if (ready_ .EQ. OFF_ERR_) then

         call AllocateInstance
            
         Me%ModelName = ModelName
            
         !Associates External Instances
         Me%ObjTime            = AssociateInstance (mTIME_           , ComputeTimeID     )
         Me%ObjHorizontalGrid  = AssociateInstance (mHORIZONTALGRID_ , HorizontalGridID  )
         Me%ObjHorizontalMap   = AssociateInstance (mHORIZONTALMAP_  , HorizontalMapID   )
         Me%ObjGridData        = AssociateInstance (mGRIDDATA_       , GridDataID        )
         Me%ObjBasinGeometry   = AssociateInstance (mBASINGEOMETRY_  , BasinGeometryID   )

         !Time Stuff
         call GetComputeTimeLimits   (Me%ObjTime, BeginTime = Me%BeginTime,           &
                                       EndTime = Me%EndTime, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructSnow - ModuleSnow - ERR010'

         call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructSnow - ModuleSnow - ERR020'

         call ReadLockExternalVar (StaticOnly = .false.)

         !Gets the size of the grid
         call GetHorizontalGridSize (Me%ObjHorizontalGrid,                            &
                                     Size     = Me%Size,                              &
                                     WorkSize = Me%WorkSize,                          &
                                     STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructSnow - ModuleSnow - ERR030'
                      
         call AllocateVariables
         call ReadDataFile         
         call InitializeVariables
         
         if (Me%Continuous) call OpenInitialFile
         call ConstructProperties
         if (Me%Continuous) call CloseInitialFile
         
         if (Me%OutPut%Yes) then
            call ConstructHDF5Output
         endif

         !call CalculateTotalStoredVolume

         !Output Results
         if (Me%OutPut%Yes) then               
            call SnowOutput
         endif

         call ReadUnLockExternalVar (StaticOnly = .false.)

         !Returns ID
         SnowID = Me%InstanceID
      
         STAT_ = SUCCESS_

      else cd0
            
         stop 'ConstructSnow - ModuleSnow - ERR040' 

      end if cd0

      if (present(STAT)) STAT = STAT_

      !----------------------------------------------------------------------

   end subroutine ConstructSnow
 
   !-------------------------------------------------------------------------   
                              
   subroutine AllocateInstance

      !Arguments-------------------------------------------------------------
                                                    
      !Local-----------------------------------------------------------------
      type (T_Snow), pointer :: NewObjSnow
      type (T_Snow), pointer :: PreviousObjSnow

      !Allocates new instance
      allocate (NewObjSnow)
      nullify  (NewObjSnow%Next)

      !Insert New Instance into list and makes Current point to it
      if (.not. associated(FirstObjSnow)) then
         FirstObjSnow => NewObjSnow
         Me           => NewObjSnow
      else
         PreviousObjSnow => FirstObjSnow
         Me              => FirstObjSnow%Next
         do while (associated(Me))
               PreviousObjSnow => Me
               Me              => Me%Next
         enddo
         Me                   => NewObjSnow
         PreviousObjSnow%Next => NewObjSnow
      endif

      Me%InstanceID = RegisterNewInstance (mSnow_)


   end subroutine AllocateInstance

   !-------------------------------------------------------------------------

   subroutine ReadDataFile

      !Arguments-------------------------------------------------------------

      !Local-----------------------------------------------------------------        
      integer                                     :: STAT_CALL
      type(T_PropertyID)                          :: InitialSnowColumnID
      type(T_PropertyID)                          :: DailyAverageTemperatureID
      type(T_PropertyID)                          :: AlbedoID
      type(T_PropertyID)                          :: ForestCoverFractionID
      type(T_PropertyID)                          :: SlopeFactorID
      integer                                     :: iflag
      integer                                     :: ClientNumber
      logical                                     :: BlockFound
      integer                                     :: i, j
      real                                        :: dummy

      !Reads the name of the data file from nomfich
      call ReadFileName ('SNOW_DATA', Me%Files%DataFile, "Snow Data File", STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR010'

      !Reads the name of the transient HDF file from nomfich
      call ReadFileName ('SNOW_HDF', Me%Files%TransientHDF, "Snow HDF File", STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR020'

      call ReadFileName ('SNOW_FIN', Me%Files%FinalFile, Message = "Snow Final File", STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR030'

      !Constructs the DataFile
      call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR040'

      !Continuous Computation
      call GetData (Me%Continuous,               &
                    Me%ObjEnterData, iflag,      &
                    SearchType   = FromFile,     &
                    keyword      = 'CONTINUOUS', &
                    default      = .false.,      &
                    ClientModule = 'ModuleSnow', &
                    STAT         = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR050'
      
      !if (Me%Continuous) then
      !   call ReadFileName ('SNOW_INI', Me%Files%InitialFile, Message = "Snow Initial File", STAT = STAT_CALL)
      !   if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR060'
      !endif      
      
      !Gets DT for computing Snow Melting. Strongly advise to use a day (86400 seconds)
      call GetData (Me%SnowMeltingDT,                 &
                    Me%ObjEnterData, iflag,           &
                    SearchType   = FromFile,          &
                    keyword      = 'SNOW_MELTING_DT', &
                    default      = 86400.0,           &
                    ClientModule = 'ModuleSnow',      &
                    STAT         = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR070'      
      
      !if (.not. Me%Continuous) then
      !   !Gets Initial snow column Block 
      !   call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,       &
      !                              '<BeginInitialSnowColumn>',           &
      !                              '<EndInitialSnowColumn>', BlockFound, &
      !                              STAT = STAT_CALL)
      !   if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR080'
      !
      !   if (BlockFound) then
      !      call ConstructFillMatrix (PropertyID       = InitialSnowColumnID,   &
      !                                EnterDataID      = Me%ObjEnterData,       &
      !                                TimeID           = Me%ObjTime,            &
      !                                HorizontalGridID = Me%ObjHorizontalGrid,  &
      !                                ExtractType      = FromBlock,             &
      !                                PointsToFill2D   = Me%ExtVar%BasinPoints, &
      !                                Matrix2D         = Me%SnowPack%OldValue,  &
      !                                TypeZUV          = TypeZ_,                &
      !                                STAT             = STAT_CALL)
      !      if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR090'
      !
      !      call KillFillMatrix(InitialWaterColumnID%ObjFillMatrix, STAT = STAT_CALL)
      !      if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR100'
      !   else
      !      write(*,*)'Missing Block <BeginInitialSnowColumn> / <EndInitialSnowColumn>' 
      !      stop      'ReadDataFile - ModuleSnow - ERR110'
      !   endif
      !endif
      
   end subroutine ReadDataFile
   
   !--------------------------------------------------------------------------

   subroutine AllocateVariables

      !Arguments-------------------------------------------------------------
      !Local----------------------------------------------------------------- 
      !Begin-----------------------------------------------------------------       

      allocate (Me%SnowMeltingFlux (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
      Me%SnowMeltingFlux = 0.0
           
   end subroutine AllocateVariables
       
   !-------------------------------------------------------------------------
   
   subroutine InitializeVariables

      !Arguments-------------------------------------------------------------        
      !Local-----------------------------------------------------------------         
      !Begin-----------------------------------------------------------------       
   
   end subroutine InitializeVariables
       
   !-------------------------------------------------------------------------
   
   subroutine ConstructHDF5Output

      !Arguments-------------------------------------------------------------

      !Local-----------------------------------------------------------------
      integer                                             :: ILB, IUB, JLB, JUB    
      integer                                             :: STAT_CALL
      integer                                             :: HDF5_CREATE

      !Bounds
      ILB = Me%WorkSize%ILB
      IUB = Me%WorkSize%IUB

      JLB = Me%WorkSize%JLB
      JUB = Me%WorkSize%JUB

      call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

      !Opens HDF File
      call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%TransientHDF)//"5", HDF5_CREATE, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSnow - ERR010'

      !Write the Horizontal Grid
      call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSnow - ERR020'

      !Sets limits for next write operations
      call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSnow - ERR030'

      !Writes the Grid
      call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                           Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSnow - ERR040'

      call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",                   &
                           Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSnow - ERR050'

      !Writes the River Points
      call HDF5WriteData   (Me%ObjHDF5, "/Grid", "RiverPoints", "-",                   &
                           Array2D = Me%ExtVar%RiverPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSnow - ERR060'

      !Flushes All pending HDF5 commands
      call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSnow - ERR070'

   end subroutine ConstructHDF5Output

   !-------------------------------------------------------------------------   
   
   subroutine ConstructProperties
   
      !Arguments-------------------------------------------------------------
   
      !Local-----------------------------------------------------------------
      integer                             :: ClientNumber
      integer                             :: STAT_CALL
      logical                             :: BlockFound
      type (T_Property), pointer          :: NewProperty

      !Begin-----------------------------------------------------------------
      
do1 : do      
         call ExtractBlockFromBuffer (Me%ObjEnterData,                   &
                                      ClientNumber    = ClientNumber,    &
                                      block_begin     = "BeginProperty", &
                                      block_end       = "EndProperty",   &
                                      BlockFound      = BlockFound,      &
                                      STAT            = STAT_CALL)
cd1 :    if (STAT_CALL .EQ. SUCCESS_) then    

cd2 :       if (BlockFound) then                                                  

               !Construct a New Property 
               Call ConstructProperty (NewProperty)

               !Add new Property to the SoilProperties List 
               Call AddProperty (NewProperty)   
   
            else cd2

               call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
               if (STAT_CALL .NE. SUCCESS_) &
                  stop 'ConstructProperties - ModuleSnow - ERR010'
               exit do1    !No more blocks            
                
            end if cd2

         else 
                
            stop 'ConstructProperties - ModuleSnow - ERR020'
            
         end if cd1 
      enddo do1
      
      !Check if all the properties were provided
      if (.not. associated(Me%SnowPack)) &
         stop 'ConstructProperties - ModuleSnow - ERR030'
      if (.not. associated(Me%Albedo)) &
         stop 'ConstructProperties - ModuleSnow - ERR040'
      if (.not. associated(Me%DailyAvgTemp)) &
         stop 'ConstructProperties - ModuleSnow - ERR050'
      if (.not. associated(Me%ForestCoverFraction)) &
         stop 'ConstructProperties - ModuleSnow - ERR060'
      if (.not. associated(Me%SlopeFactor)) &
         stop 'ConstructProperties - ModuleSnow - ERR070'
      
   end subroutine ConstructProperties
   
   !-------------------------------------------------------------------------   
   
   subroutine AddProperty (NewProperty)

      !Arguments-------------------------------------------------------------
      type(T_Property), pointer :: NewProperty

      !----------------------------------------------------------------------

      if (.not.associated(Me%FirstProperty)) then
         Me%PropertiesNumber     = 1
         Me%FirstProperty        => NewProperty
         Me%LastProperty         => NewProperty
      else
         NewProperty%Prev        => Me%LastProperty
         Me%LastProperty%Next    => NewProperty
         Me%LastProperty         => NewProperty
         Me%PropertiesNumber     = Me%PropertiesNumber + 1
      end if 

   end subroutine AddProperty 

   !--------------------------------------------------------------------------    
   
   subroutine ConstructProperty (NewProperty)

      !Arguments-------------------------------------------------------------
      type(T_Property), pointer           :: NewProperty

      !External--------------------------------------------------------------
      integer                             :: STAT_CALL

      !----------------------------------------------------------------------
             
      allocate (NewProperty, STAT = STAT_CALL)            
      if(STAT_CALL .NE. SUCCESS_) stop 'ConstructProperty - ModuleSnow - ERR010'
              
      call ConstructPropertyID        (NewProperty%ID, Me%ObjEnterData, FromBlock)
      call ConstructPropertyState     (NewProperty)
      call ConstructPropertyEvolution (NewProperty)
      call ConstructPropertyValues    (NewProperty)
      call ConstructPropertyOutPut    (NewProperty)
      
      select case (NewProperty%ID%IDNumber)
      case (SnowPack_)
         Me%SnowPack => NewProperty
      case (Albedo_)
         Me%Albedo => NewProperty
      case (DailyAvgTemp_)
         Me%DailyAvgTemp => NewProperty
      case (ForestCoverFraction_)
         Me%ForestCoverFraction => NewProperty
      case (SnowSlopeFactor_)
         Me%SlopeFactor => NewProperty
      end select

   end subroutine ConstructProperty
    
   !-------------------------------------------------------------------------  
   
   subroutine ConstructPropertyState (NewProperty)

      !Arguments-------------------------------------------------------------
      type(T_property), pointer       :: NewProperty

      !External--------------------------------------------------------------
      integer                         :: STAT_CALL, iflag
      !----------------------------------------------------------------------        

      !!<BeginKeyword>
      !   !Keyword          : PARTICULATE
      !   !<BeginDescription>
      !   !<EndDescription>
      !   !Type             : logical   
      !   !Default          : Dissolved
      !   !File keyword     : SEDPROP
      !   !Multiple Options : 1 (.true.), 0 (.false.)
      !   !Search Type      : From Block
      !   !Begin Block      : <beginproperty>
      !   !End Block        : <endproperty>
      !!<EndKeyword>
      !
      !call GetData(NewProperty%Particulate,                                            &
      !            Me%ObjEnterData,  iflag,                                            &
      !            SearchType   = FromBlock,                                           &
      !            keyword      = 'PARTICULATE',                                       &
      !            ClientModule = 'ModuleSnow',                       &
      !            STAT         = STAT_CALL)
      !if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModuleSnow - ERR01'
      !  
      !if (NewProperty%Particulate)then
      !   if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
      !         write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is not'
      !         write(*,*) 'recognised as PARTICULATE'
      !         stop 'Construct_PropertyState - ModuleSnow - ERR03'
      !   end if
      !endif
                
   end subroutine ConstructPropertyState

   !--------------------------------------------------------------------------   
   
   subroutine ConstructPropertyEvolution (NewProperty)

      !Arguments-------------------------------------------------------------
      type(T_property), pointer                   :: NewProperty

      !External--------------------------------------------------------------       
      !Local-----------------------------------------------------------------
      !----------------------------------------------------------------------
      
      NewProperty%Evolution%NextCompute = Me%ExtVar%Now
            
   end subroutine ConstructPropertyEvolution     

   !--------------------------------------------------------------------------   
   
   subroutine ConstructPropertyValues (NewProperty)

      !Arguments-------------------------------------------------------------
      type(T_property),              pointer      :: NewProperty

      !External--------------------------------------------------------------
      integer                                     :: STAT_CALL, i, j

      !Local-----------------------------------------------------------------
      integer                                     :: iflag
      integer                                     :: ILB,IUB
      integer                                     :: JLB,JUB
      integer                                     :: WorkSizeILB, WorkSizeIUB
      integer                                     :: WorkSizeJLB, WorkSizeJUB
      
      !Begin-----------------------------------------------------------------        
      !Boundaries
      ILB = Me%Size%ILB
      IUB = Me%Size%IUB
      JLB = Me%Size%JLB
      JUB = Me%Size%JUB

      WorkSizeILB = Me%WorkSize%ILB
      WorkSizeIUB = Me%WorkSize%IUB
      WorkSizeJLB = Me%WorkSize%JLB
      WorkSizeJUB = Me%WorkSize%JUB

      allocate (NewProperty%Value (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
      if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModuleSnow - ERR010'
      NewProperty%Value(:,:) = FillValueReal

      allocate (NewProperty%ValueOld (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
      if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModuleSnow - ERR020'
      NewProperty%ValueOld(:,:) = FillValueReal

      !This variable is a logic one is true if the property is old
      !and the user wants to continue the run with results of a previous run.
      call GetData (NewProperty%Old,                          &
                    Me%ObjEnterData, iflag,                   &
                    keyword      = 'OLD',                     &
                    Default      = (.not. Me%Continuous),     &                        
                    SearchType   = FromBlock,                 &
                    ClientModule = 'ModuleSnow',              &
                    STAT         = STAT_CALL)              
      if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModuleSnow - ERR030'
          
      !The property can't be "OLD" if it's not a continuation run
      if ((.not. Me%Continuous) .and. NewProperty%Old) then
         write (*,*) 'Property ', trim(NewProperty%ID%Name), &
            ' has OLD set to TRUE, but the CONTINOUS file keyword is missing or set to FALSE'
         stop 'ConstructPropertyValues - ModuleSnow - ERR031'
      endif
      
      ! if the property is not 'OLD' the property values in the domain are initialized
      ! if it's true ('OLD') this same values are read from the final file of the previous run
      if (.not. NewProperty%Old) then

         !Get water points
         call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleSnow - ERR040'

         call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                    EnterDataID          = Me%ObjEnterData,                  &
                                    TimeID               = Me%ObjTime,                       &
                                    HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                    ExtractType          = FromBlock,                        &
                                    PointsToFill2D       = Me%ExtVar%BasinPoints,            &
                                    Matrix2D             = NewProperty%Value,                &
                                    TypeZUV              = TypeZ_,                           &
                                    STAT                 = STAT_CALL)
         if (STAT_CALL /= SUCCESS_)                                                          &
               stop 'ConstructPropertyValues - ModuleSnow - ERR050'

         if(.not. NewProperty%ID%SolutionFromFile)then

               call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
               if (STAT_CALL /= SUCCESS_)&
                  stop 'ConstructPropertyValues - ModuleSnow - ERR060'
         end if

         call UnGetBasin(Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL) 
         if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModuleSnow - ERR070'

      else

         ! If the property is old then the program is going to try to find a property
         ! with the same name in the Water properties initial file written in HDF format  
         call ReadOldValueFromHDF (NewProperty)

      end if   

   end subroutine ConstructPropertyValues

   !-------------------------------------------------------------------------   
   
   subroutine ConstructPropertyOutPut (NewProperty)

      !Arguments-------------------------------------------------------------
      type(T_Property),    pointer        :: NewProperty

      !Local-----------------------------------------------------------------
      integer                             :: STAT_CALL, iflag

      !Begin-----------------------------------------------------------------

      call GetData (NewProperty%TimeSerie,         &
                    Me%ObjEnterData, iflag,        &
                    Keyword      = 'TIME_SERIE',   &
                    ClientModule = 'ModuleSnow',   &
                    Default      = .false.,        &
                    SearchType   = FromBlock,      &
                    STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) &
         stop 'ConstructPropertyOutPut - ModuleSnow - ERR010'
      
      call GetData (NewProperty%OutputHDF,         &
                    Me%ObjEnterData, iflag,        &
                    Keyword      = 'OUTPUT_HDF',   &
                    ClientModule = 'ModuleSnow',   &
                    Default      = .false.,        &
                    SearchType   = FromBlock,      &
                    STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) &
         stop 'ConstructPropertyOutPut - ModuleSnow - ERR020'      
        
      !call GetData(NewProperty%BoxTimeSerie,                                           &
      !            Me%ObjEnterData, iflag,                                             &
      !            Keyword      = 'BOX_TIME_SERIE',                                    &
      !            Default      = .false.,                                             &
      !            SearchType   = FromBlock,                                           &
      !            ClientModule = 'ModuleSnow',                            &
      !            STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) &
      !   stop 'Construct_PropertyOutPut - ModuleSnow - ERR02'
      !
      !if (NewProperty%BoxTimeSerie) then
      !   Me%Output%Boxes_ON = .true.
      !   Me%NumberPropForBoxes = Me%NumberPropForBoxes + 1
      !endif
      !
      !call GetData(NewProperty%BoxTimeSerie2D,                                           &
      !            Me%ObjEnterData, iflag,                                               &
      !            Keyword      = 'BOX_TIME_SERIE2D',                                    &
      !            Default      = .false.,                                               &
      !            SearchType   = FromBlock,                                             &
      !            ClientModule = 'ModuleSnow',                              &
      !            STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) &
      !   stop 'Construct_PropertyOutPut - ModuleSnow - ERR03'
        
   end subroutine ConstructPropertyOutPut
   
   !------------------------------------------------------------------------- 
   
   subroutine ReadOldValueFromHDF (NewProperty)

      !Arguments-------------------------------------------------------------
      type(T_Property), pointer                   :: NewProperty

      !External--------------------------------------------------------------
      integer                                     :: STAT_CALL

      !Local-----------------------------------------------------------------
      character (Len=StringLength)                :: PropertyName                               

      !----------------------------------------------------------------------

      PropertyName = trim(adjustl(NewProperty%ID%Name))

      call HDF5ReadData (Me%ObjIniHDF5, "/Results/"//trim(adjustl(NewProperty%ID%Name)), &
                         trim(adjustl(NewProperty%ID%Name)),                       &
                         Array2D = NewProperty%Value,                              &
                         STAT    = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) &
         stop 'ReadOldValueFromHDF - ModuleSnow - ERR010'

   end subroutine ReadOldValueFromHDF

   !-------------------------------------------------------------------------   
   
   subroutine OpenInitialFile

      !Arguments-------------------------------------------------------------

      !Local-----------------------------------------------------------------
      integer                                     :: STAT_CALL
      integer                                     :: WorkILB, WorkIUB
      integer                                     :: WorkJLB, WorkJUB
      integer                                     :: HDF5_READ

      !----------------------------------------------------------------------

      WorkILB = Me%WorkSize%ILB 
      WorkIUB = Me%WorkSize%IUB 
      WorkJLB = Me%WorkSize%JLB 
      WorkJUB = Me%WorkSize%JUB       
      
      call ReadFileName ('SNOW_INI', Me%Files%InitialFile, "Snow Initial File", STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenInitialFile - ModuleSnow - ERR010'      

      !Gets File Access Code
      call GetHDF5FileAccess (HDF5_READ = HDF5_READ)

      Me%ObjIniHDF5 = 0

      !Opens HDF5 File
      call ConstructHDF5 (Me%ObjIniHDF5,              &
                          trim(Me%Files%InitialFile), &
                          HDF5_READ, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) &
            stop 'OpenInitialFile - ModuleSnow - ERR020'

      ! Reads from HDF file the Property concentration and open boundary values
      call HDF5SetLimits (Me%ObjIniHDF5,     &
                          WorkILB, WorkIUB,  &
                          WorkJLB, WorkJUB,  &
                          STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) &
            stop 'OpenInitialFile - ModuleSnow - ERR030'      
            
   end subroutine OpenInitialFile
   
   !-------------------------------------------------------------------------
   
   subroutine CloseInitialFile
   
      !Arguments-------------------------------------------------------------

      !Local-----------------------------------------------------------------
      integer                                     :: STAT_CALL   
   
      !----------------------------------------------------------------------
      
      call KillHDF5 (Me%ObjIniHDF5, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) &
            stop 'CloseInitialFile - ModuleSnow - ERR010'
   
   end subroutine CloseInitialFile
   
   !-------------------------------------------------------------------------
   
   subroutine ConstructTimeSerie

      !Arguments-------------------------------------------------------------

      !Local-----------------------------------------------------------------
      character(len=StringLength), dimension(:), pointer  :: PropertyList
      integer                                             :: nProperties
      integer                                             :: STAT_CALL
      integer                                             :: iflag
      character(len=StringLength)                         :: TimeSerieLocationFile
      type (T_Property), pointer                          :: PropertyX
      integer                                             :: n
      integer                                             :: TimeSerieNumber, dn, Id, Jd
      real                                                :: CoordX, CoordY
      logical                                             :: CoordON, IgnoreOK
      character(len=StringLength)                         :: TimeSerieName
        
      !Begin------------------------------------------------------------------
        
      !Counts the number of Properties which has timeserie option set to true
      PropertyX => Me%FirstProperty
      nProperties = 0
      do while (associated(PropertyX))
         if (PropertyX%TimeSerie) then
               nProperties = nProperties + 1
         endif
         PropertyX => PropertyX%Next
      enddo
      
      !Allocates PropertyList
      allocate(PropertyList(nProperties))
        
      !Property names
      n=1
      PropertyX  => Me%FirstProperty
      do while (associated(PropertyX))
         if (PropertyX%TimeSerie) then
               PropertyList(n)  = trim(PropertyX%ID%Name)//'['//trim(PropertyX%ID%Units)//']'
               n=n+1
         endif
         PropertyX=>PropertyX%Next
      enddo
        
      call GetData (TimeSerieLocationFile,                &
                    Me%ObjEnterData, iflag,               &
                    SearchType   = FromFile,              &
                    keyword      = 'TIME_SERIE_LOCATION', &
                    ClientModule = 'ModuleSnow',          &
                    Default      = Me%Files%DataFile,     &
                    STAT         = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR010' 

      if (iflag == 1) then
         Me%OutPut%TimeSerie_ON = .true.
      else
         Me%OutPut%TimeSerie_ON = .false.
      endif
        
      !Get water points
      call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR020'

      !Constructs TimeSerie
      call StartTimeSerie (Me%ObjTimeSerie, Me%ObjTime,           &
                           TimeSerieLocationFile,                 &
                           PropertyList, "srsn",                  &
                           WaterPoints2D = Me%ExtVar%BasinPoints, &
                           STAT          = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR030' 

      call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR040'

      !Deallocates PropertyList
      deallocate(PropertyList)

      !Corrects if necessary the cell of the time serie based in the time serie coordinates
      call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR050'

      do dn = 1, TimeSerieNumber

         call GetTimeSerieLocation (Me%ObjTimeSerie, dn, &  
                                    CoordX   = CoordX,   &
                                    CoordY   = CoordY,   & 
                                    CoordON  = CoordON,  &
                                    STAT     = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR060'
            
         call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR070'
            
i1:      if (CoordON) then
            call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR080'

            if (Id < 0 .or. Jd < 0) then
                
               call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
               if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR090'

               if (IgnoreOK) then
                  write(*,*) 'Time Serie outside the domain - ',trim(TimeSerieName)
                  cycle
               else
                  stop 'ConstructTimeSerie - ModuleSnow - ERR100'
               endif

            endif

            call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR110'

         endif i1

         call GetTimeSerieLocation(Me%ObjTimeSerie, dn,     &  
                                    LocalizationI   = Id,   &
                                    LocalizationJ   = Jd,   & 
                                    STAT     = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR120'

         call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR130'

         if (Me%ExtVar%BasinPoints(Id, Jd) /= WaterPoint) then
            write(*,*) 'Time Serie in a cell outside basin - ',trim(TimeSerieName)
         endif
            
         call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSnow - ERR0140'

      enddo
       
   end subroutine ConstructTimeSerie

   !-------------------------------------------------------------------------

   subroutine ConstructHDF
        
      !External--------------------------------------------------------------
      type (T_Property), pointer                  :: CurrentProperty
      logical                                     :: OutputON
      integer :: STAT_CALL

      !Begin-----------------------------------------------------------------

      nullify(Me%OutPut%OutTime)

      OutputON = OFF

      CurrentProperty => Me%FirstProperty
      do while (associated(CurrentProperty))            
         if(CurrentProperty%OutputHDF) OutputON = ON
         CurrentProperty => CurrentProperty%Next
      enddo

      if(OutputON)then
  
         call GetOutPutTime (Me%ObjEnterData,                   &
                             CurrentTime = Me%BeginTime,        &
                             EndTime     = Me%EndTime,          &
                             keyword     = 'OUTPUT_TIME',       &
                             SearchType  = FromFile,            &
                             OutPutsTime = Me%OutPut%OutTime,   &
                             OutPutsOn   = Me%OutPut%HDF_ON,    &
                             OutPutsNumber = Me%OutPut%Number,  &
                             STAT        = STAT_CALL)

         if (STAT_CALL /= SUCCESS_) &
               stop 'ConstructHDF - ModuleSnow - ERR010' 

         if (Me%OutPut%HDF_ON) then
               Me%OutPut%NextOutPut = 1
               call OpenHDF5OutPutFile
         else
               write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
               write(*,*)'one property has HDF format outputs.'
               stop 'ConstructHDF - ModuleSnow - ERR020'
         endif 

         !Output for restart
         call GetOutPutTime (Me%ObjEnterData,                           &
                             CurrentTime  = Me%ExtVar%Now,              &
                             EndTime      = Me%EndTime,                 &
                             keyword      = 'RESTART_FILE_OUTPUT_TIME', &
                             SearchType   = FromFile,                   &
                             OutPutsTime  = Me%OutPut%RestartOutTime,   &
                             OutPutsOn    = Me%OutPut%WriteRestartFile, &
                             STAT         = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR030'

      endif

   end subroutine ConstructHDF  
   
   !-------------------------------------------------------------------------
 
   subroutine OpenHDF5OutPutFile        

      !Local-----------------------------------------------------------------
      integer                                             :: ILB,IUB,JLB,JUB    
      integer                                             :: STAT_CALL
      integer                                             :: HDF5_CREATE
      
      !Begin-----------------------------------------------------------------

      !Bounds
      ILB = Me%WorkSize%ILB
      IUB = Me%WorkSize%IUB

      JLB = Me%WorkSize%JLB
      JUB = Me%WorkSize%JUB

      !Gets a pointer to Topography
      call GetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR010'

      call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR020'

      call GetHDF5FileAccess (HDF5_CREATE = HDF5_CREATE)

      !Opens HDF File
      call ConstructHDF5 (Me%ObjHDF5, trim(Me%Files%TransientHDF)//"5", HDF5_CREATE, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR030'
      
      !Write the Horizontal Grid
      call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR040'


      !Sets limits for next write operations
      call HDF5SetLimits (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR050'
        
      !Writes the Grid
      call HDF5WriteData (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                           Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR060'

      !WriteBasinPoints
      call HDF5WriteData (Me%ObjHDF5, "/Grid", "BasinPoints", "-",          &
                           Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR070'

      !Flushes All pending HDF5 commands
      call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR080'  

      !Unget
      call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR090'  

      !UnGets Topography
      call UnGetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5OutPutFile - ModuleSnow - ERR100'

   end subroutine OpenHDF5OutPutFile         
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   subroutine GetSnowMelting (ObjSnowID, SnowMelted, STAT)
   
      !Arguments-------------------------------------------------------------
      integer                                         :: ObjSnowID
      real, dimension(:, :), pointer                  :: SnowMelted
      integer, intent(OUT), optional                  :: STAT

      !Local-----------------------------------------------------------------
      integer                                         :: STAT_, ready_

      !---------------------------------------------------------------------- 
      
      STAT_ = UNKNOWN_

      call Ready (ObjSnowID, ready_) 
        
      if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
          (ready_ .EQ. READ_LOCK_ERR_)) then
         
         call Read_Lock(mSnow_, Me%InstanceID) 

         SnowMelted => Me%SnowPack%Value

         STAT_ = SUCCESS_
      else
         STAT_ = ready_
      end if

      if (present(STAT))STAT = STAT_      
      
   end subroutine GetSnowMelting
   
   !--------------------------------------------------------------------------
   
   subroutine UnGetSnowMelting2D_R4 (ObjSnowID, Array, STAT)

      !Arguments-------------------------------------------------------------
      integer                                         :: ObjSnowID
      real(4), dimension(:, :), pointer               :: Array
      integer, intent(OUT), optional                  :: STAT

      !Local-----------------------------------------------------------------
      integer                                         :: STAT_, ready_

      !----------------------------------------------------------------------

      STAT_ = UNKNOWN_

      call Ready(ObjSnowID, ready_)

      if (ready_ .EQ. READ_LOCK_ERR_) then

         nullify(Array)
         call Read_Unlock(mSNOW_, Me%InstanceID, "UnGetSnowMelted2D_R4")

         STAT_ = SUCCESS_
      else               
         STAT_ = ready_
      end if

      if (present(STAT)) STAT = STAT_

   end subroutine UnGetSnowMelting2D_R4

   !--------------------------------------------------------------------------

   subroutine UnGetSnowMelting2D_R8(ObjSnowID, Array, STAT)

      !Arguments-------------------------------------------------------------
      integer                                         :: ObjSnowID
      real(8), dimension(:, :), pointer               :: Array
      integer, intent(OUT), optional                  :: STAT

      !Local-----------------------------------------------------------------
      integer                                         :: STAT_, ready_

      !----------------------------------------------------------------------

      STAT_ = UNKNOWN_

      call Ready(ObjSnowID, ready_)

      if (ready_ .EQ. READ_LOCK_ERR_) then

         nullify(Array)
         call Read_Unlock(mSNOW_, Me%InstanceID, "UnGetSnowMelted2D_R8")

         STAT_ = SUCCESS_
      else               
         STAT_ = ready_
      end if

      if (present(STAT)) STAT = STAT_

   end subroutine UnGetSnowMelting2D_R8
        
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
   subroutine ModifySnow(SnowID, STAT)

      !Arguments-------------------------------------------------------------
      integer                                     :: SnowID
      integer, intent(OUT), optional              :: STAT

      !Local-----------------------------------------------------------------
      integer                                     :: STAT_, ready_
      integer                                     :: STAT_CALL
      real                                        :: M
      logical                                     :: Restart
      integer                                     :: Niter, iter
      integer                                     :: n_restart
      integer                                     :: ILB, IUB, JLB, JUB, I, J
      logical                                     :: IsFinalFile      
      !----------------------------------------------------------------------

      STAT_ = UNKNOWN_

      call Ready(SnowID, ready_)

      if (ready_ .EQ. IDLE_ERR_) then

         if (MonitorPerformance) call StartWatch ("ModuleSnow", "ModifySnow")

         ILB = Me%WorkSize%ILB
         IUB = Me%WorkSize%IUB
         JLB = Me%WorkSize%JLB
         JUB = Me%WorkSize%JUB         
                  
         !Time Stuff
         call GetComputeCurrentTime (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ModifySnow - ModuleSnow - ERR010'
         
         call GetComputeTimeStep (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'ModifySnow - ModuleSnow - ERR020'

         !Checks to see if it's time to compute a new SnowMeltingFlux
         if(Me%ExtVar%Now .GE. Me%SnowPack%Evolution%NextCompute) then

            !read the properties values 
            call ModifyProperties
            
            do I = ILB, IUB
            do J = JLB, JUB
            
               if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                      
                  M = 4.0 * (1 - Me%Albedo%Value(i, j)) * exp(-4 * Me%ForestCoverFraction%Value(i, j)) * Me%SlopeFactor%Value(i, j)
                  Me%SnowMeltingFlux(i, j) = max(0.0, M * (Me%DailyAvgTemp%Value(i, j) - Me%MeltingTemperature)) / Me%SnowMeltingDT

               else

                  Me%SnowMeltingFlux(i, j) = FillValueReal

               endif

            enddo
            enddo            
            
            Me%SnowPack%Evolution%NextCompute = Me%ExtVar%Now + Me%SnowMeltingDT
            
         endif
                  
         do I = ILB, IUB
         do J = JLB, JUB
            
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                      
               Me%SnowMelted(i, j) = max(0.0, Me%SnowPack%Value(i, j) - Me%SnowMeltingFlux(i, j) * Me%ExtVar%DT)

            else

               Me%SnowMelted(i, j) = FillValueReal

            endif

         enddo
         enddo           
                        
         !Output Results
         if (Me%OutPut%Yes) then                   
            call SnowOutput
         endif
            
!         if (Me%Output%BoxFluxes) then
!               call ComputeBoxesWaterFluxes
!         endif

         !Restart Output
         if (Me%Output%WriteRestartFile .and. .not. (Me%ExtVar%Now == Me%EndTime)) then
               if(Me%ExtVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                  IsFinalFile = .false.
                  call WriteFinalFile(IsFinalFile)
                  Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
               endif
         endif

         STAT_ = SUCCESS_
         if (MonitorPerformance) call StopWatch ("ModuleSnow", "ModifySnow")

      else               
         STAT_ = ready_
      end if

      if (present(STAT)) STAT = STAT_

   end subroutine ModifySnow
    
   !--------------------------------------------------------------------------   
   
   subroutine ModifyProperties
 
      !Local-----------------------------------------------------------------
      integer :: i, j      
      type (T_Property), pointer                  :: PropertyX
      integer                                     :: STAT_CALL
      
      !Begin-----------------------------------------------------------------
        
      PropertyX => Me%FirstProperty
      do while (associated(PropertyX))

         if (PropertyX%ID%SolutionFromFile) then

               call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,   &
                                      Matrix2D       = PropertyX%Value,              &
                                      PointsToFill2D = Me%ExtVar%BasinPoints,        &
                                      STAT           = STAT_CALL)
               if (STAT_CALL /= SUCCESS_) then
                  write (*,*) "ATTENTION"
                  write (*,*) "Was not possible to read property '", trim(PropertyX%ID%Name), "' from file."  
                  stop 'ModifyProperties - ModuleSnow - ERR010'
               endif
         endif
                      
         PropertyX => PropertyX%Next
        
      enddo

   end subroutine ModifyProperties

   !--------------------------------------------------------------------------
    
   subroutine SnowOutput

      !Arguments-------------------------------------------------------------

      !Local-----------------------------------------------------------------
      integer                                     :: STAT_CALL
      integer                                     :: ILB, IUB, JLB, JUB
      real, dimension(6)  , target                :: AuxTime
      real, dimension(:)  , pointer               :: TimePointer       

      if (MonitorPerformance) call StartWatch ("ModuleSnow", "SnowOutput")

      !Bounds
      ILB = Me%WorkSize%ILB
      IUB = Me%WorkSize%IUB

      JLB = Me%WorkSize%JLB
      JUB = Me%WorkSize%JUB

      if (Me%ExtVar%Now >= Me%OutPut%OutTime(Me%OutPut%NextOutPut)) then

         !Writes current time
         call ExtractDate   (Me%ExtVar%Now , AuxTime(1), AuxTime(2),         &
                                             AuxTime(3), AuxTime(4),         &
                                             AuxTime(5), AuxTime(6))
         TimePointer => AuxTime

         call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'SnowOutput - ModuleSnow - ERR010'

         call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                              "YYYY/MM/DD HH:MM:SS",                         &
                              Array1D      = TimePointer,                    &
                              OutputNumber = Me%OutPut%NextOutPut,           &
                              STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'SnowOutput - ModuleSnow - ERR020'

         !Sets limits for next write operations
         call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'SnowOutput - ModuleSnow - ERR030'
         
         !Writes everything to disk
         call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'SnowOutput - ModuleSnow - ERR99'

         Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

      endif

      if (MonitorPerformance) call StopWatch ("ModuleSnow", "SnowOutput")
        
   end subroutine SnowOutput

   !-------------------------------------------------------------------------
    
   subroutine WriteFinalFile(IsFinalFile)
              
      !Arguments
      logical                                     :: IsFinalFile
      !Local-----------------------------------------------------------------
      type (T_Property), pointer                  :: PropertyX
      integer                                     :: STAT_CALL      
      integer                                     :: HDF5_CREATE
      character(LEN = PathLength)                 :: FileName
      integer                                     :: ObjHDF5
      real, dimension(6), target                  :: AuxTime
      real, dimension(:), pointer                 :: TimePtr
      type (T_Time)                               :: Actual           
      real                                        :: Total_Mass_Created
      character (Len = StringLength)              :: str_mass_created, string_to_be_written         
      !Begin----------------------------------------------------------------

      !Gets a pointer to Topography
      call GetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR010'

      call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR020'

         !Gets File Access Code
      call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

      !Checks if it's at the end of the run 
      !or !if it's supposed to overwrite the final HDF file
      !if ((Me%ExtVar%Now == Me%EndTime) .or. Me%Output%RestartOverwrite) then
      if (IsFinalFile .or. Me%Output%RestartOverwrite) then
         filename = trim(Me%Files%FinalFile)
      else
         FileName = ChangeSuffix(Me%Files%FinalFile, &
                           "_"//trim(TimeToString(Me%ExtVar%Now))//".fin")
      endif

      ObjHDF5 = 0
      !Opens HDF5 File
      call ConstructHDF5 (ObjHDF5,        &
                          trim(filename), &
                          HDF5_CREATE, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) &
         stop 'WriteFinalFile - ModuleSnow - ERR030'

      Actual = Me%ExtVar%Now
         
      call ExtractDate (Actual, AuxTime(1), AuxTime(2), AuxTime(3), &
                                AuxTime(4), AuxTime(5), AuxTime(6))
      !Writes Time
      TimePtr => AuxTime
      call HDF5SetLimits (ObjHDF5, 1, 6, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR040'

      call HDF5WriteData (ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS", &
                          Array1D = TimePtr, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR050'

      !Sets limits for next write operations
      call HDF5SetLimits (ObjHDF5,           &
                          Me%WorkSize%ILB,   &
                          Me%WorkSize%IUB,   &
                          Me%WorkSize%JLB,   &
                          Me%WorkSize%JUB,   &
                          STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR060'

      !Write the Horizontal Grid
      call WriteHorizontalGrid(Me%ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR070'

      !Writes the Grid
      call HDF5WriteData   (ObjHDF5, "//Grid/Topography", "Topography", "m", &
                           Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR080'

      !WriteBasinPoints
      call HDF5WriteData   (ObjHDF5, "//Grid/BasinPoints", "BasinPoints", "-", &
                           Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR090'


      PropertyX => Me%FirstProperty
      do while (associated(PropertyX))

         call HDF5SetLimits (ObjHDF5,         &
                             Me%WorkSize%ILB, &
                             Me%WorkSize%IUB, &
                             Me%WorkSize%JLB, &
                             Me%WorkSize%JUB, &
                             STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR100'

         call HDF5WriteData (ObjHDF5,                              &
                             "/Results/"//trim(PropertyX%ID%Name), &
                             trim(PropertyX%ID%Name),              &
                             trim(PropertyX%ID%Units),             &
                             Array2D = PropertyX%Value,            &
                             STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR110'
                
         PropertyX => PropertyX%Next

      enddo

      !Writes everything to disk
      call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR030'

      !Unget
      call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR90'  

      !UnGets Topography
      call UnGetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleSnow - ERR100'
            

   end subroutine WriteFinalFile

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine KillSnow(SnowID, STAT)

      !Arguments---------------------------------------------------------------
      integer                             :: SnowID              
      integer, optional, intent(OUT)      :: STAT

      !External----------------------------------------------------------------
      integer                             :: ready_              

      !Local-------------------------------------------------------------------
      integer                             :: STAT_, nUsers, STAT_CALL    
      !character(len=StringLength)         :: MassErrorFile
      logical                             :: IsFinalFile

      !------------------------------------------------------------------------

      STAT_ = UNKNOWN_

      call Ready(SnowID, ready_)    

cd1 : if (ready_ .NE. OFF_ERR_) then


         nUsers = DeassociateInstance(mSnow_,  Me%InstanceID)

         if (nUsers == 0) then

            !Writes file with final condition
            IsFinalFile = .true.
            call WriteFinalFile(IsFinalFile)

!                !Writes Mass Error
!                call ReadFileName("ROOT_SRT", MassErrorFile, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleSnow - ERR02a'
!                MassErrorFile = trim(adjustl(MassErrorFile))//"MassError.dat"
!                
!                call WriteGridData  (MassErrorFile,                            &
!                     COMENT1          = "MassErrorFile",                       &
!                     COMENT2          = "MassErrorFile",                       &
!                     HorizontalGridID = Me%ObjHorizontalGrid,                  &
!                     FillValue        = -99.0,                                 &
!                     OverWrite        = .true.,                                &
!                     GridData2D_Real  = Me%MassError,                          &
!                     STAT             = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'KillSnow - Snow - ERR00'
!        
!                if(Me%WriteMaxFlowModulus) then
!                    call WriteGridData  (Me%MaxFlowModulusFile,                &
!                         COMENT1          = "MaxFlowModulusFile",              &
!                         COMENT2          = "MaxFlowModulusFile",              &
!                         HorizontalGridID = Me%ObjHorizontalGrid,              &
!                         FillValue        = -99.0,                             &
!                         OverWrite        = .true.,                            &
!                         GridData2D_Real  = Me%MaxFlowModulus,                 &
!                         STAT             = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'KillSnow - Snow - ERR00'
!                endif
!                
!                if (Me%WriteMaxWaterColumn) then
!                    call WriteGridData  (Me%MaxWaterColumnFile,                &
!                         COMENT1          = "MaxWaterColumnFile",              &
!                         COMENT2          = "MaxWaterColumnFile",              &
!                         HorizontalGridID = Me%ObjHorizontalGrid,              &
!                         FillValue        = -99.0,                             &
!                         OverWrite        = .true.,                            &
!                         GridData2D_Real  = Me%MaxWaterColumn,                 &
!                         STAT             = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'KillSnow - Snow - ERR00'
!                endif
!
!
!                if (Me%ObjDrainageNetwork /= 0) then
! 
!!                    if(Me%WriteMaxWaterColumn) call WriteChannelsLevelData
!
!                    nUsers = DeassociateInstance (mDRAINAGENETWORK_, Me%ObjDrainageNetwork)
!                    if (nUsers == 0) stop 'KillSnow - Snow - ERR01'
!                endif

            if (Me%OutPut%Yes) then
               call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
               if (STAT_CALL /= SUCCESS_) stop 'KillSnow - ModuleSnow - ERR01'
            endif
                
            !if (Me%Discharges) then   
            !   call Kill_Discharges(Me%ObjDischarges, STAT = STAT_CALL)
            !   if (STAT_CALL /= SUCCESS_) stop 'KillSnow - ModuleSnow - ERR02'
            !endif

            !if (Me%ImposeBoundaryValue .and. Me%BoundaryImposedLevelInTime) then
            !        
            !   if (Me%ImposedLevelTS%TimeSerie%ObjTimeSerie /= 0) then
            !      call KillTimeSerie(Me%ImposedLevelTS%TimeSerie%ObjTimeSerie, STAT = STAT_CALL)
            !      if (STAT_CALL /= SUCCESS_) stop 'KillSnow - ModulePorousMedia - ERR03' 
            !   endif
            !
            !endif
            !
            !if (Me%Output%BoxFluxes) then
            !   call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
            !   if (STAT_CALL /= SUCCESS_)                               &
            !      stop 'KillSnow - Snow - ERR04'
            !endif
                
            !Deassociates External Instances
            nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
            if (nUsers == 0) stop 'KillSnow - Snow - ERR05'

            nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
            if (nUsers == 0) stop 'KillSnow - Snow - ERR06'

            nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjGridData)
            if (nUsers == 0) stop 'KillSnow - Snow - ERR07'

            nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
            if (nUsers == 0) stop 'KillSnow - Snow - ERR08'

            nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
            if (nUsers == 0) stop 'KillSnow - Snow - ERR09'
                
            !deallocate(Me%myWaterColumnOld)
            !    
            !deallocate (Me%iFlowX)
            !deallocate (Me%iFlowY)
            !deallocate (Me%lFlowX)
            !deallocate (Me%lFlowY)
            !deallocate (Me%iFlowToChannels)
            !deallocate (Me%lFlowToChannels)
            !deallocate (Me%lFlowBoundary)
            !deallocate (Me%iFlowBoundary)
            !deallocate (Me%iFlowRouteDFour)
            !
            !nullify    (Me%iFlowX)
            !nullify    (Me%iFlowY)
            !nullify    (Me%lFlowX)
            !nullify    (Me%lFlowY)
            !nullify    (Me%iFlowToChannels)
            !nullify    (Me%lFlowToChannels)
            !nullify    (Me%lFlowBoundary)
            !nullify    (Me%iFlowBoundary)
            !nullify    (Me%iFlowRouteDFour)


            !Deallocates Instance
            call DeallocateInstance ()

            SnowID   = 0
            STAT_      = SUCCESS_

         end if

      end if cd1

      if (present(STAT)) STAT = STAT_

      !------------------------------------------------------------------------

   end subroutine KillSnow

   !------------------------------------------------------------------------
    
   subroutine DeallocateInstance ()

      !Arguments-------------------------------------------------------------

      !Local-----------------------------------------------------------------
      type (T_Snow), pointer                    :: AuxObjSnow
      type (T_Snow), pointer                    :: PreviousObjSnow

      !Updates pointers
      if (Me%InstanceID == FirstObjSnow%InstanceID) then
         FirstObjSnow => FirstObjSnow%Next
      else
         PreviousObjSnow => FirstObjSnow
         AuxObjSnow      => FirstObjSnow%Next
         do while (AuxObjSnow%InstanceID /= Me%InstanceID)
               PreviousObjSnow => AuxObjSnow
               AuxObjSnow      => AuxObjSnow%Next
         enddo

         !Now update linked list
         PreviousObjSnow%Next => AuxObjSnow%Next

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

   subroutine Ready (SnowID, ready_) 

      !Arguments-------------------------------------------------------------
      integer                                     :: SnowID
      integer                                     :: ready_

      !----------------------------------------------------------------------

      nullify (Me)

cd1:    if (SnowID > 0) then
         call LocateObjSnow (SnowID)
         ready_ = VerifyReadLock (mSNOW_, Me%InstanceID)
      else
         ready_ = OFF_ERR_
      end if cd1

      !----------------------------------------------------------------------

   end subroutine Ready

   !--------------------------------------------------------------------------

   subroutine LocateObjSnow (ObjSnowID)

      !Arguments-------------------------------------------------------------
      integer                                     :: ObjSnowID

      !Local-----------------------------------------------------------------

      Me => FirstObjSnow
      do while (associated (Me))
         if (Me%InstanceID == ObjSnowID) exit
         Me => Me%Next
      enddo

      if (.not. associated(Me)) stop 'ModuleSnow - LocateObjSnow - ERR010'

   end subroutine LocateObjSnow

   !--------------------------------------------------------------------------

   subroutine ReadLockExternalVar (StaticOnly)
        
      !Arguments-------------------------------------------------------------
      logical                                     :: StaticOnly

      !Local-----------------------------------------------------------------
      integer                                     :: STAT_CALL

      !Time Stuff
      call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR01'

      call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR02'

      !Gets Basin Points
      call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR03'
        
      !!Gets cell slope
      !call GetCellSlope   (Me%ObjBasinGeometry, Me%ExtVar%CellSlope, STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR04'
      !
      !!Gets River Points
      !call GetRiverPoints (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR05'
      !
      !!Gets Horizontal Grid
      !call GetHorizontalGrid(Me%ObjHorizontalGrid,                                     &
      !                       DUX    = Me%ExtVar%DUX,    DVY    = Me%ExtVar%DVY,        &
      !                       DXX    = Me%ExtVar%DXX,    DYY    = Me%ExtVar%DYY,        &
      !                       DZX    = Me%ExtVar%DZX,    DZY    = Me%ExtVar%DZY,        &
      !                       XX2D_Z = Me%ExtVar%XX2D_Z, YY2D_Z = Me%ExtVar%YY2D_Z,     &
      !                       STAT   = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR06'
      !
      !call GetGridCellArea  (Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea,             &
      !                       STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR06a'
      !
      !!Gets a pointer to Topography
      !call GetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR07'
      !
      !if (.not. StaticOnly) then
      !
      !    !Gets Boundary Points
      !    call GetBoundaries    (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
      !    if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSnow - ERR10'
      !
      !endif

   end subroutine ReadLockExternalVar
   
   !--------------------------------------------------------------------------
   
   subroutine ReadUnLockExternalVar(StaticOnly)
        
      !Arguments-------------------------------------------------------------
      logical                                     :: StaticOnly
        
      !Local-----------------------------------------------------------------
      integer                                     :: STAT_CALL

      !Unget Basin Points
      call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
      if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR01'

      !!Unget River Points
      !call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR02'
      !
      !!Unget Cell Slope
      !call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%CellSlope, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR02a'
      !
      !!Unget Horizontal Grid
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUX, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR03'
      !
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVY, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR04'
      !
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DXX, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR05'
      !
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DYY, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR06'
      !
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR05'
      !
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR06'
      !
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%XX2D_Z, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR07'
      !
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%YY2D_Z, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR08'
      !
      !call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR09'
      !
      !!Ungets the Topography
      !call UngetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
      !if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR10'
      !
      !if (.not. StaticOnly) then
      !
      !   call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
      !   if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSnow - ERR11'
      !
      !endif 
        
   end subroutine ReadUnLockExternalVar   
   
   !--------------------------------------------------------------------------
    
end module ModuleSnow   