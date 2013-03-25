!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : Basin Geometry
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to Calculate the delimition of geometric variables of a
!               : hidrographic basin
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


Module ModuleBasinGeometry

    use ModuleGlobalData
    use ModuleEnterData          
    use ModuleFunctions,    only : RODAXY
    use ModuleGridData,     only : ConstructGridData, WriteGridData, GetGridData,        &
                                   GetFillValue, UnGetGridData, KillGridData      
    use ModuleHorizontalGrid
    use ModuleHDF5

    implicit none
    private

    !Subroutine----------------------------------------------------------------


    !Constructor
    public  :: ConstructBasinGeometry
    private ::      AllocateInstance
    private ::      ReadDataFile 
    private ::      AllocateVariables
    private ::      DrainageDirection
    private ::          ConvertDirection
    private ::      FindDepressions
    private ::          AddLake
    private ::          AddLakePoint
    private ::      FlatDepressions
    private ::          LakeHasUndrainedPoints
    private ::      CalculateCellSlope
    private ::      DelineateBasin
    private ::      UpStreamAreas
    private ::      NodeDefinition
    private ::      WriteBasinASCII
    private ::      WriteBasinHDF

    !Selector
    public  :: GetBasinPoints
    public  :: GetMicroChannels
    public  :: GetRiverPoints
    public  :: GetDrainageDirection    
    public  :: GetCellSlope
    public  :: UnGetBasin

    public  :: TargetPoint

    !Destructor
    public  ::  KillBasinGeometry
    private ::      DeallocateInstance

    !Management
    private ::  Ready
    private ::      LocateObjBasinGeometry
    private ::  GetExternalVar
    private ::  UnGetExternalVar

    interface UnGetBasin
        module procedure  UnGetBasin2D_I
        module procedure  UnGetBasin2D_R4
        module procedure  UnGetBasin2D_R8
    end interface UnGetBasin

    !Parameter-----------------------------------------------------------------
    integer, parameter                              :: N_           = 1
    integer, parameter                              :: NE_          = 2
    integer, parameter                              :: E_           = 3
    integer, parameter                              :: SE_          = 4
    integer, parameter                              :: S_           = 5
    integer, parameter                              :: SW_          = 6
    integer, parameter                              :: W_           = 7
    integer, parameter                              :: NW_          = 8
    integer, parameter                              :: NODRAINAGE_  = 9

    !To be Moved to GLOBAL DATA
    integer, parameter                              :: NoRiverPoint = 0
    integer, parameter                              :: RiverPoint   = 1

    !Type----------------------------------------------------------------------
    type T_ExternalVar
        real, dimension(:, :), pointer              :: Topography
        real, dimension(:), pointer                 :: XX_Z, YY_Z
        real, dimension(:,:), pointer               :: XX2D_Z, YY2D_Z
        real, dimension(:, :), pointer              :: CellArea
        real                                        :: GridFillValue
    end type T_ExternalVar

   type T_LakePoint
        integer                                     :: I, J
        type (T_LakePoint), pointer                 :: Next
    end type T_LakePoint

    type T_Lake
        integer                                     :: LakeID
        integer                                     :: ExitI, ExitJ
        type (T_LakePoint), pointer                 :: FirstLakePoint
        type (T_Lake), pointer                      :: Next
    end type T_Lake

    type T_LakeExit
        integer                                     :: PointI, PointJ
    end type T_LakeExit

    type T_UpstreamAreas
        integer                                     :: I, J
        real(8)                                     :: Area
        type (T_UpstreamAreas), pointer             :: Next                     => null()
    end type T_UpstreamAreas


    type T_BasinGeometry
        integer                                     :: InstanceID
        type (T_Size2D)                             :: Size
        type (T_Size2D)                             :: WorkSize
        type (T_ExternalVar)                        :: ExtVar
        integer                                     :: OutletI, OutletJ
        real                                        :: RiverTresholdArea
        integer                                     :: nDepressions            = 0
        logical                                     :: TopographyChanged = .false.
        type (T_Lake), pointer                      :: FirstLake
        type (T_UpstreamAreas), pointer             :: FirstUpstreamArea
        integer, dimension(:, :), pointer           :: BasinPoints
        integer, dimension(:, :), pointer           :: RiverPoints
        integer, dimension(:, :), pointer           :: DrainageDirection
        integer, dimension(:, :), pointer           :: BoundaryPoints
        real,    dimension(:, :), pointer           :: UpStreamArea
        integer, dimension(:, :), pointer           :: NodeIDs
        integer, dimension(:, :), pointer           :: LakePoints
        real,    dimension(:, :), pointer           :: NewTopography
        real,    dimension(:, :), pointer           :: CellSlope
        real,    dimension(:, :), pointer           :: MicroChannels
        logical                                     :: Reservoirs
        character(Pathlength)                       :: ReservoirFile
        character(Pathlength)                       :: DataFile
        character(PathLength)                       :: TransientHDF
        logical                                     :: WriteHDF5
        logical                                     :: WriteReaches
        logical                                     :: WriteDelineation
        logical                                     :: WriteDrainedArea
        logical                                     :: WriteDrainageDirection
        logical                                     :: WriteCellSlope
        logical                                     :: WriteBasinPoints
        logical                                     :: RemoveDepressions
        logical                                     :: FlatSolution
        logical                                     :: DelineateBasin
        real                                        :: SlopeDepressions     = 1.e-4
        character(Pathlength)                       :: DrainedAreaFile
        character(Pathlength)                       :: ReachesFile
        character(Pathlength)                       :: DelineationFile
        character(Pathlength)                       :: DrainageDirectionFile
        character(Pathlength)                       :: CellSlopeFile
        character(Pathlength)                       :: NewTopographyFile
        character(Pathlength)                       :: BasinPointsFile
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjTopography        = 0
        type (T_BasinGeometry), pointer             :: Next
    end type T_BasinGeometry


    !Global Module Variables
    type (T_BasinGeometry), pointer                 :: FirstBasinGeometry
    type (T_BasinGeometry), pointer                 :: Me

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructBasinGeometry (BasinGeometryID, GridDataID, HorizontalGridID, &
                                       DataFile, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: BasinGeometryID
        integer                                     :: GridDataID
        integer                                     :: HorizontalGridID
        character(len=*), optional                  :: DataFile
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mBasinGeometry_)) then
            nullify (FirstBasinGeometry)
            call RegisterModule (mBasinGeometry_) 
        endif
 
        call Ready(BasinGeometryID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance       ()

            nullify (Me%BasinPoints      )
            nullify (Me%RiverPoints      )
            nullify (Me%DrainageDirection)
            nullify (Me%BoundaryPoints   )
            nullify (Me%UpStreamArea     )
            nullify (Me%NodeIDs          )
            nullify (Me%LakePoints       )
            nullify (Me%FirstLake        )
            nullify (Me%NewTopography    )
            nullify (Me%CellSlope        ) 
            nullify (Me%MicroChannels    )

            !Associates other instances
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjTopography     = AssociateInstance (mGRIDDATA_,       GridDataID      )

            if (present(DataFile)) then
                Me%DataFile = DataFile
            else
                !Reads the name of the data file from nomfich
                call ReadFileName ('BASIN_GEOMETRY', Me%DataFile, "Basin Data File", STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructBasinGeometry - ModuleBasinGeometry - ERR01'
            endif

            call GetExternalVar

            call AllocateVariables
            
            call ReadDataFile

            call ConstructVariables
            
            !This have to be changed in the future to a more elegant way.
            if (Me%FlatSolution) then
                call SetInteriorBasinPoints
                call UnGetExternalVar

                !Returns ID
                BasinGeometryID = Me%InstanceID

                STAT_ = SUCCESS_
                return
            endif

            Me%TopographyChanged = .TRUE.

            do while (Me%TopographyChanged) 
                
                Me%TopographyChanged     = .FAlSE.
                Me%DrainageDirection     = null_int
                Me%LakePoints            = NoLakePoint

                call DrainageDirection

                call FindDepressions

                call FlatDepressions
                
                call RemoveAllDepressions

                !Stops if Topography contains Depressions
                if (.not. Me%RemoveDepressions .and. Me%TopographyChanged) then
                    write (*,*) 'Topography contains depressions'
                    write (*,*) 'Please remove them first'
                    write (*,*) 'ModuleBasinGeometry - ConstructBasinGeometry - ERR02'
                    stop 99
                end if

            end do

            !Stops if Topography contains Depressions
            if (Me%RemoveDepressions) then

                call WriteGridData  (Me%NewTopographyFile,                                 &
                                     COMENT1          = "New Topography",                  &
                                     COMENT2          = "Topography without depressions",  &
                                     HorizontalGridID = Me%ObjHorizontalGrid,              &
                                     FillValue        = -99.0,                             &
                                     OverWrite        = .true.,                            &
                                     GridData2D_Real  = Me%NewTopography,                  &
                                     STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleBasinGeometry - ConstructBasinGeometry - ERR03'
                write(*,*)'New Topography file written : ', trim(Me%NewTopographyFile)
                call UnGetExternalVar

                !Returns ID
                BasinGeometryID = Me%InstanceID

                return
            endif


            call CalculateCellSlope

            call UpStreamAreas

            call NodeDefinition

            if (Me%DelineateBasin) then
                call DelineateBasin
            else
                call SetInteriorBasinPoints
            endif

            call WriteBasinASCII

            if (Me%WriteHDF5) then
                call WriteBasinHDF
            endif

            call UnGetExternalVar

            !Returns ID
            BasinGeometryID = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            
            stop 'ModuleBasinGeometry - ConstructBasinGeometry - ERR99' 

        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructBasinGeometry 

    !--------------------------------------------------------------------------

    subroutine AllocateInstance ()

        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        type (T_BasinGeometry), pointer             :: NewObjBasinGeometry
        type (T_BasinGeometry), pointer             :: PreviousObjBasinGeometry


        !Allocates new instance
        allocate (NewObjBasinGeometry)
        nullify  (NewObjBasinGeometry%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstBasinGeometry)) then
            FirstBasinGeometry            => NewObjBasinGeometry
            Me                            => NewObjBasinGeometry
        else
            PreviousObjBasinGeometry      => FirstBasinGeometry
            Me                            => FirstBasinGeometry%Next
            do while (associated(Me))
                PreviousObjBasinGeometry  => Me
                Me                        => Me%Next
            enddo
            Me                            => NewObjBasinGeometry
            PreviousObjBasinGeometry%Next => NewObjBasinGeometry
        endif

        Me%InstanceID = RegisterNewInstance (mBASINGEOMETRY_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadDataFile ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ObjEnterData = 0
        integer                                     :: flag         
        integer                                     :: STAT_CALL
        real                                        :: MicroChannelPer
        logical                                     :: UpStreamFound
        integer                                     :: ClientNumber
        type (T_UpstreamAreas), pointer             :: NewUpstreamArea          => null()                      
        type (T_UpstreamAreas), pointer             :: PreviousUpstreamArea     => null()                      
        type (T_UpstreamAreas), pointer             :: CurrentUpstreamArea      => null()                                      
        

        !Opens Data file
        call ConstructEnterData (ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR02'

        !Flat solution (does do anything then setting boundary points)
        call GetData (Me%FlatSolution, ObjEnterData, flag,                               &
                      SearchType   = FromFile,                                           &
                      keyword      ='FLAT_SOLUTION',                                     &
                      default      = .false.,                                            &
                      ClientModule ='ModuleBasinGeometry',                               &
                      STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR05'


        !Remove Depressions
        call GetData (Me%RemoveDepressions, ObjEnterData, flag,                          &
                      SearchType   = FromFile,                                           &
                      keyword      ='REMOVE_DEPRESSIONS',                                &
                      default      = .false.,                                            &
                      ClientModule ='ModuleBasinGeometry',                               &
                      STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR07'

        if (Me%RemoveDepressions) then

            !Delineation File 
            call GetData (Me%NewTopographyFile, ObjEnterData, flag,                      &
                          SearchType   = FromFile,                                       &
                          keyword      ='NEW_TOPOGRAPHY_FILE',                           &
                          ClientModule ='ModuleBasinGeometry',                           &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_ .or. flag == 0) stop 'ReadDataFile - ModuleBasinGeometry - ERR08'

            !Gets Minimum Slope in Depressions
            call GetData (Me%SlopeDepressions, ObjEnterData, flag,                       &
                          SearchType   = FromFile,                                       &
                          keyword      ='SLOPE_DEPRESSIONS',                             &
                          default      = 1.e-3,                                          &
                          ClientModule ='ModuleBasinGeometry',                           &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR08a'

        else

            !Gets minimum area for a grid cell to be marked as river.
            call GetData (Me%RiverTresholdArea, ObjEnterData, flag,                          &
                          SearchType   = FromFile,                                           &
                          keyword      ='TRESHOLD_AREA',                                     &
                          ClientModule ='ModuleBasinGeometry',                               &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_ .or. flag == 0) stop 'ReadDataFile - ModuleBasinGeometry - ERR04z'

            !Micro Channels Percentage
            call GetData (MicroChannelPer, ObjEnterData, flag,                               &
                          SearchType   = FromFile,                                           &
                          keyword      ='MICRO_CHANNEL',                                     &
                          ClientModule ='ModuleBasinGeometry',                               &
                          default      = 1.0,                                                &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR04y'
            
            Me%MicroChannels = MicroChannelPer

            !Consider Reservoirs
            call GetData (Me%Reservoirs, ObjEnterData, flag,                                 &
                          SearchType   = FromFile,                                           &
                          keyword      ='RESERVOIRS',                                        &
                          default      = .false.,                                            &
                          ClientModule ='ModuleBasinGeometry',                               &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR07a'

            if (Me%Reservoirs) then
                call GetData (Me%ReservoirFile, ObjEnterData, flag,                          &
                              SearchType   = FromFile,                                       &
                              keyword      ='RESERVOIRS_FILE',                               &
                              ClientModule ='ModuleBasinGeometry',                           &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR07b'
            endif

            !Write HDF 5
            call GetData (Me%WriteHDF5, ObjEnterData, flag,                                  &
                          SearchType   = FromFile,                                           &
                          keyword      ='OUTPUT_HDF5',                                       &
                          default      = .false.,                                            &
                          ClientModule ='ModuleBasinGeometry',                               &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_ ) stop 'ReadDataFile - ModuleBasinGeometry - ERR04b'


            !Drained Area
            call GetData (Me%WriteDrainedArea, ObjEnterData, flag,                           &
                          SearchType   = FromFile,                                           &
                          keyword      ='DRAINED_AREA',                                      &
                          default      = .false.,                                            &
                          ClientModule ='ModuleBasinGeometry',                               &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_ ) stop 'ReadDataFile - ModuleBasinGeometry - ERR04a'

            !Gets the name of the Drained Area file
            if (Me%WriteDrainedArea) then
                call GetData (Me%DrainedAreaFile, ObjEnterData, flag,                        &
                              SearchType   = FromFile,                                       &
                              keyword      ='DRAINED_AREA_FILE',                             &
                              ClientModule ='ModuleBasinGeometry',                           &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR04b'
            endif

            !Write Drainage Direction
            call GetData (Me%WriteDrainageDirection, ObjEnterData, flag,                     &
                          SearchType   = FromFile,                                           &
                          keyword      ='DRAINAGE_DIRECTION',                                &
                          default      = .false.,                                            &
                          ClientModule ='ModuleBasinGeometry',                               &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_ ) stop 'ReadDataFile - ModuleBasinGeometry - ERR04a'

            !Gets the name of the Drained Area file
            if (Me%WriteDrainageDirection) then
                call GetData (Me%DrainageDirectionFile, ObjEnterData, flag,                  &
                              SearchType   = FromFile,                                       &
                              keyword      ='DRAINAGE_DIRECTION_FILE',                       &
                              ClientModule ='ModuleBasinGeometry',                           &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR04b'
            endif

            !Write Cell Slope
            call GetData (Me%WriteCellSlope, ObjEnterData, flag,                             &
                          SearchType   = FromFile,                                           &
                          keyword      ='CELL_SLOPE',                                        &
                          default      = .false.,                                            &
                          ClientModule ='ModuleBasinGeometry',                               &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_ ) stop 'ReadDataFile - ModuleBasinGeometry - ERR04a'

            !Gets the name of the Drained Area file
            if (Me%WriteCellSlope) then
                call GetData (Me%CellSlopeFile, ObjEnterData, flag,                          &
                              SearchType   = FromFile,                                       &
                              keyword      ='CELL_SLOPE_FILE',                               &
                              ClientModule ='ModuleBasinGeometry',                           &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR04b'
            endif

            !Make Delineation
            call GetData (Me%DelineateBasin, ObjEnterData, flag,                             &
                          SearchType   = FromFile,                                           &
                          keyword      ='DELINEATE_BASIN',                                   &
                          default      = .false.,                                            &
                          ClientModule ='ModuleBasinGeometry',                               &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR07'
            

            !Write Reaches ? 
            call GetData (Me%WriteReaches, ObjEnterData, flag,                               &
                          SearchType   = FromFile,                                           &
                          keyword      ='WRITE_REACHES',                                     &
                          default      = .false.,                                            &
                          ClientModule ='ModuleBasinGeometry',                               &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR09'

            if (Me%WriteReaches) then

                !Reaches File 
                call GetData (Me%ReachesFile, ObjEnterData, flag,                            &
                              SearchType   = FromFile,                                       &
                              keyword      ='REACHES_FILE',                                  &
                              ClientModule ='ModuleBasinGeometry',                           &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR10'

            endif


            if (Me%DelineateBasin) then

                !Gets Outlet Location - I
                call GetData (Me%OutletI, ObjEnterData, flag,                                    &
                              SearchType   = FromFile,                                           &
                              keyword      ='OUTLET_I',                                          &
                              ClientModule ='ModuleBasinGeometry',                               &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_ .or. flag == 0) stop 'ReadDataFile - ModuleBasinGeometry - ERR03'


                !Gets Outlet Location - J
                call GetData (Me%OutletJ, ObjEnterData, flag,                                    &
                              SearchType   = FromFile,                                           &
                              keyword      ='OUTLET_J',                                          &
                              ClientModule ='ModuleBasinGeometry',                               &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_ .or. flag == 0) stop 'ReadDataFile - ModuleBasinGeometry - ERR04'


                !Write Delineation
                call GetData (Me%WriteDelineation, ObjEnterData, flag,                           &
                              SearchType   = FromFile,                                           &
                              keyword      ='WRITE_DELINEATION',                                 &
                              default      = .false.,                                            &
                              ClientModule ='ModuleBasinGeometry',                               &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR05'

                if (Me%WriteDelineation) then

                    !Delineation File 
                    call GetData (Me%DelineationFile, ObjEnterData, flag,                        &
                                  SearchType   = FromFile,                                       &
                                  keyword      ='DELINEATION_FILE',                              &
                                  ClientModule ='ModuleBasinGeometry',                           &
                                  STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR06'

                    call GetData (Me%BasinPointsFile, ObjEnterData, flag,                        &
                                  SearchType   = FromFile,                                       &
                                  keyword      ='BASIN_POINTS_FILE',                             &
                                  ClientModule ='ModuleBasinGeometry',                           &
                                  STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR06a'
                    
                    if (flag == 1) then
                        Me%WriteBasinPoints = .true.
                    else
                        Me%WriteBasinPoints = .false.
                    endif                                        
                    
                endif

            endif

        endif

        !Reads Predefined Upstream Areas
        call RewindBuffer(ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR07'
        
        
doUA:   do
        
            call ExtractBlockFromBuffer(ObjEnterData,                                           &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<BeginUpstreamArea>',                &
                                        block_end       = '<EndUpstreamArea>',                  &
                                        BlockFound      = UpStreamFound,                        &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. UpStreamFound) then
            
                
                allocate(NewUpstreamArea)
                nullify (NewUpstreamArea%Next)
                
                call GetData(NewUpstreamArea%I, ObjEnterData, flag,                             &
                             SearchType   = FromBlock,                                          &
                             keyword      ='I',                                                 &
                             ClientModule ='ModuleBasinGeometry',                               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR08'

                call GetData(NewUpstreamArea%J, ObjEnterData, flag,                             &
                             SearchType   = FromBlock,                                          &
                             keyword      ='J',                                                 &
                             ClientModule ='ModuleBasinGeometry',                               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR09'
            
                call GetData(NewUpstreamArea%Area, ObjEnterData, flag,                          &
                             SearchType   = FromBlock,                                          &
                             keyword      ='AREA',                                              &
                             ClientModule ='ModuleBasinGeometry',                               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR10'

                if (.not. associated(Me%FirstUpstreamArea)) then
                    Me%FirstUpstreamArea => NewUpstreamArea
                else
                    PreviousUpstreamArea => Me%FirstUpstreamArea
                    CurrentUpstreamArea  => PreviousUpstreamArea%Next
                    do while (associated(CurrentUpstreamArea))
                        PreviousUpstreamArea => CurrentUpstreamArea
                        CurrentUpstreamArea  => PreviousUpstreamArea%Next
                    enddo
                    PreviousUpstreamArea%Next => NewUpstreamArea
                endif

            else
            
                exit doUA
            
            endif
            
        
        enddo doUA



        !Closes Data File
        call KillEnterData      (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasinGeometry - ERR99'


    end subroutine ReadDataFile

    !--------------------------------------------------------------------------

    subroutine AllocateVariables ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Gets the size of the grid
        call GetHorizontalGridSize (Me%ObjHorizontalGrid,                  &
                                    Size     = Me%Size,                    &
                                    WorkSize = Me%WorkSize,                &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateVariables - ModuleBasinGeometry - ERR01'

    
        allocate (Me%DrainageDirection  (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%BasinPoints        (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%RiverPoints        (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%UpStreamArea       (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%NodeIDs            (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%LakePoints         (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%NewTopography      (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%CellSlope          (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%MicroChannels      (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%BoundaryPoints     (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        


        Me%DrainageDirection     = null_int
        Me%BasinPoints           = UndefinedPoint
        Me%RiverPoints           = UndefinedPoint
        Me%UpStreamArea          = null_real
        Me%NodeIDs               = null_int
        Me%LakePoints            = NoLakePoint
        Me%NewTopography         = Me%ExtVar%Topography
        Me%CellSlope             = null_real
        Me%MicroChannels         = null_real
        Me%BoundaryPoints        = UndefinedPoint

    end subroutine AllocateVariables 

    !--------------------------------------------------------------------------

    subroutine ConstructVariables ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB


        !WorkSize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Fills By Default all points are none boundary points
        do j = JLB, JUB
        do i = ILB, IUB
            Me%BoundaryPoints(i, j) = Not_Boundary
        enddo
        enddo

        !Along the edges are boundaries
        do i = ILB, IUB
            Me%BoundaryPoints(i, JLB-1) = Boundary
            Me%BoundaryPoints(i, JUB+1) = Boundary
        enddo    
        
        do j = JLB, JUB
            Me%BoundaryPoints(ILB-1, j) = Boundary
            Me%BoundaryPoints(IUB+1, j) = Boundary
        enddo

        
        !Where Grid data contains values equal or lower the fill value, also boundary points are considered
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%Topography(i, j) <= Me%ExtVar%GridFillValue / 2.0) then
                Me%BoundaryPoints(i, j) = Boundary
            endif
        enddo
        enddo

        if (.not. Me%RemoveDepressions .and. Me%DelineateBasin) then

            !Tests if the Outlet is outside the domain
            if (Me%OutletI < Me%WorkSize%ILB .or. Me%OutletI > Me%WorkSize%IUB .or.          &
                Me%OutletJ < Me%WorkSize%JLB .or. Me%OutletJ > Me%WorkSize%JUB) then
                write (*,*)'Outlet lies outside the domain'
                write (*,*)'OutletI, OutletJ : ', Me%OutletI, Me%OutletJ
                write (*,*)'ILB, IUB         : ', ILB, IUB
                write (*,*)'JLB, JUB         : ', JLB, JUB
                stop 'ConstructVariables - ModuleBasinGeometry - ERR01'
            endif

            !Tests if the outlet is not a boundary points
            if (Me%BoundaryPoints(Me%OutletI, Me%OutletJ) == Boundary) then
                write (*,*)'Outlet lies on the Boundary'
                write (*,*)'OutletI, OutletJ : ', Me%OutletI, Me%OutletJ
                stop 'ConstructVariables - ModuleBasinGeometry - ERR02'
            endif

        endif

    end subroutine ConstructVariables

    !--------------------------------------------------------------------------

    subroutine DrainageDirection ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: iCen, jCen
        real                                        :: MajorSlope, Slope
        integer                                     :: Direction, i, j
        real                                        :: CenX, CenY
        real                                        :: AdjX, AdjY
        integer                                     :: iAdj, jAdj
        real                                        :: dist

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


        !Calculates the drainage Direction for every cell
        do jCen = JLB, JUB
        do iCen = ILB, IUB

            MajorSlope = null_real
            Direction  = null_int

            !Coordinates of Center Point
            CenX = XCenter (iCen, jCen)
            CenY = YCenter (iCen, jCen)

            do j = -1, 1
            do i = -1, 1

                if (.not. (i == 0 .and. j == 0)) then
                    
                    iAdj = iCen + i
                    jAdj = jCen + j

                    if (Me%BoundaryPoints(iAdj, jAdj) == Not_Boundary) then

                        !Coordinates of adjcent Point
                        AdjX    = XCenter(iAdj, jAdj)
                        AdjY    = YCenter(iAdj, jAdj)

                        !Distance between points
                        dist    = sqrt ((CenX - AdjX)**2. + (CenY - AdjY)**2.)

                        !Slope
                        Slope   = (Me%NewTopography(iCen, jCen) - Me%NewTopography(iAdj, jAdj)) / dist

                    else
                        Slope   = null_real
                    endif
                endif

                if (Slope > MajorSlope .and. Slope > 0.0) then
                    MajorSlope = Slope
                    Direction  = ConvertDirection (i, j)
                endif

            enddo
            enddo


            if (Direction /= null_int) then
                Me%DrainageDirection(iCen, jCen) = Direction
            else
                Me%DrainageDirection(iCen, jCen) = NODRAINAGE_
            endif

        enddo
        enddo

    end subroutine DrainageDirection

    !--------------------------------------------------------------------------

    subroutine FindDepressions

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: iCen, jCen, iAdj, jAdj
        integer                                     :: i, j
        integer                                     :: iTarget, jTarget
        integer                                     :: iLowest, jLowest
        logical                                     :: FillingLake
        integer, dimension(3, 3)                    :: SmallLakePoints
        real,    dimension(3, 3)                    :: SmallTopography
        real                                        :: LowestValue
        type (T_Lake), pointer                      :: NewLake, CurrentLake, AuxLake
        type (T_LakePoint), pointer                 :: CurrentLakePoint, AuxLakePoint

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Set Points close to the border to drain to the outside
        do jCen = JLB, JUB
        do iCen = ILB, IUB
            
            do j = -1, 1
            do i = -1, 1

                if (.not. (i == 0 .and. j == 0)) then
                    
                    iAdj = iCen + i
                    jAdj = jCen + j

                    if (Me%BoundaryPoints(iAdj, jAdj) == Boundary) then
                        Me%DrainageDirection(iCen, jCen) = ConvertDirection(i, j)
                    endif

                endif

            enddo
            enddo

        enddo
        enddo


        !Create Depressions in every point which has no Drainage
        do jCen = JLB, JUB
        do iCen = ILB, IUB
            if (Me%DrainageDirection(iCen, jCen) == NODRAINAGE_ .and.                   &
                Me%BoundaryPoints   (iCen, jCen) == Not_Boundary) then

                !Allocates a new Lake    
                allocate(NewLake)
                nullify (NewLake%Next)
                nullify (NewLake%FirstLakePoint)

                Me%nDepressions = Me%nDepressions + 1
                NewLake%LakeID  = Me%nDepressions

                !Inserts Lake to the Lake List
                call AddLake      (NewLake)

                !Insert point to LakePoints
                call AddLakePoint (NewLake, iCen, jCen)
                Me%LakePoints (iCen, jCen) = NewLake%LakeID

            endif
        enddo
        enddo

        !For each Lake Searches the exit point
        CurrentLake => Me%FirstLake
        do while (associated(CurrentLake))
            
            FillingLake = .true.

            !Looks for lowest point which does not belong to any lake until
            !the Lake Exit is found
FillLake:   do while (FillingLake)

                LowestValue   = -null_real
                CurrentLakePoint => CurrentLake%FirstLakePoint
                do while (associated(CurrentLakePoint))

                    SmallTopography = Me%NewTopography(CurrentLakePoint%I-1:CurrentLakePoint%I+1, &
                                                       CurrentLakePoint%J-1:CurrentLakePoint%J+1)

                    SmallLakePoints = Me%LakePoints   (CurrentLakePoint%I-1:CurrentLakePoint%I+1, &
                                                       CurrentLakePoint%J-1:CurrentLakePoint%J+1)

                    do i = 1, 3
                    do j = 1, 3
                        if (SmallLakePoints(i, j) /= CurrentLake%LakeID) then
                            if (SmallTopography(i, j) < LowestValue) then
                                LowestValue = SmallTopography(i, j)
                                iLowest     = CurrentLakePoint%I + i - 2
                                jLowest     = CurrentLakePoint%J + j - 2
                            endif
                        endif
                    enddo
                    enddo

                    CurrentLakePoint => CurrentLakePoint%Next
           
                enddo

                
                !If at iLowest, jLowest exists already a lake, merge the two Depressions
                if (Me%LakePoints(iLowest, jLowest) /= NoLakePoint) then

                    AuxLake => Me%FirstLake
doRemove:           do while (associated(AuxLake))
                        
                        if (AuxLake%LakeID == Me%LakePoints(iLowest, jLowest)) then

                            AuxLakePoint => AuxLake%FirstLakePoint
                            do while (associated(AuxLakePoint))
                                
                                !Insert point to LakePoints
                                call AddLakePoint (CurrentLake, AuxLakePoint%I, AuxLakePoint%J)
                                Me%LakePoints (AuxLakePoint%I, AuxLakePoint%J) = CurrentLake%LakeID

                                AuxLakePoint => AuxLakePoint%Next
                            enddo

                            CurrentLake%ExitI = AuxLake%ExitI
                            CurrentLake%ExitJ = AuxLake%ExitJ

                            call RemoveLake (AuxLake)

                            exit doRemove
                        endif
                        AuxLake => AuxLake%Next
                    enddo doRemove

                endif


                !The following options can occur
                !1. The lowest point is on the border
                !   The Exit point is found. Stop Filling the lake.
                !2. The lowest point is NOT on the border
                !2.1 The lowest point drains to the border
                !    The Exit point is found. Stop Filling the lake.   
                !2.2 The lowest point drains to the Current Lake -> 
                !    Add the point to the current lake and continue filling the lake
                !2.2 The lowest point drains to an Undrained point (different from the
                !    Current Lake (which was intercepted by 2.1). 
                !    The Exit Point is found. Stop Filling the lake

                !1
                if (Me%BoundaryPoints(iLowest, jLowest) == Boundary) then
                    CurrentLake%ExitI = iLowest
                    CurrentLake%ExitJ = jLowest
                    FillingLake       = .false.
                else
                    
                    !Target point of the lowest point
                    call TargetPoint (Me%DrainageDirection(iLowest, jLowest), iLowest, jLowest, iTarget, jTarget)
                    
doDrain:            do
    
                        !2.1
                        if (Me%BoundaryPoints(iTarget, jTarget) == Boundary) then
                            CurrentLake%ExitI = iLowest
                            CurrentLake%ExitJ = jLowest
                            FillingLake       = .false.
                            exit doDrain
                        endif

                        !2.2
                        if (Me%LakePoints(iTarget, jTarget) == CurrentLake%LakeID) then
                            call AddLakePoint (CurrentLake, iLowest, jLowest)
                            Me%LakePoints (iLowest, jLowest) = CurrentLake%LakeID
                            exit doDrain
                        endif

                        !2.2
                        if (Me%DrainageDirection(iTarget, jTarget) == NODRAINAGE_) then

                            CurrentLake%ExitI = iLowest
                            CurrentLake%ExitJ = jLowest
                            FillingLake       = .false.
                            exit doDrain

                        endif
                        
                        !Target point of the target point
                        call TargetPoint (Me%DrainageDirection(iTarget, jTarget), iTarget, jTarget, iTarget, jTarget)

                    enddo doDrain

                endif

            enddo FillLake

            CurrentLake => CurrentLake%Next

        enddo

    end subroutine FindDepressions

    !--------------------------------------------------------------------------

    subroutine FlatDepressions

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: iCen, jCen, i, j, iAdj, jAdj
        real                                        :: LakeExitDepth, dx, dy, dz
        type (T_Lake), pointer                      :: CurrentLake
        type (T_LakePoint), pointer                 :: CurrentLakePoint

        !Flats the depressions
        CurrentLake => Me%FirstLake
        do while (associated(CurrentLake))

            LakeExitDepth = Me%NewTopography(CurrentLake%ExitI, CurrentLake%ExitJ)

            CurrentLakePoint => CurrentLake%FirstLakePoint
            do while (associated(CurrentLakePoint))

                if (Me%NewTopography(CurrentLakePoint%I, CurrentLakePoint%J) /= LakeExitDepth) then
                    Me%NewTopography(CurrentLakePoint%I, CurrentLakePoint%J) = LakeExitDepth
                    Me%TopographyChanged = .true.
                endif

                CurrentLakePoint => CurrentLakePoint%Next
            enddo

            CurrentLake => CurrentLake%Next
        enddo

        !Sets Drainage Direction inside Flat Areas
        CurrentLake => Me%FirstLake
        do while (associated(CurrentLake))

            !Set Drainage Direction to NODRIANAGE
            CurrentLakePoint => CurrentLake%FirstLakePoint
            do while (associated(CurrentLakePoint))
                Me%DrainageDirection(CurrentLakePoint%I, CurrentLakePoint%J) = NODRAINAGE_
                CurrentLakePoint => CurrentLakePoint%Next
            enddo

            !Set Drainage Direction close to lake exit to lake exit
            CurrentLakePoint => CurrentLake%FirstLakePoint
            do while (associated(CurrentLakePoint))

                iCen = CurrentLakePoint%I
                jCen = CurrentLakePoint%J

                do j = -1, 1
                do i = -1, 1

                    if (.not. (i == 0 .and. j == 0)) then
                        iAdj = iCen + i
                        jAdj = jCen + j
                        if (iAdj == CurrentLake%ExitI .and. jAdj == CurrentLake%ExitJ) then
                            Me%DrainageDirection(iCen, jCen) = ConvertDirection(i, j)
                            dx = XCenter (iCen, jCen) - XCenter (iAdj, jAdj)
                            dy = YCenter (iCen, jCen) - YCenter (iAdj, jAdj)
                            dz = sqrt(dx**2.+dy**2.) * Me%SlopeDepressions
                            Me%NewTopography    (iCen, jCen) = Me%NewTopography(iAdj, jAdj) + dz
                            Me%TopographyChanged = .true.
                        endif
                    endif
                enddo
                enddo

                CurrentLakePoint => CurrentLakePoint%Next
            enddo

            do while (LakeHasUndrainedPoints(CurrentLake))

                !Set Drainage Direction inside lake to lake exit
                CurrentLakePoint => CurrentLake%FirstLakePoint
                do while (associated(CurrentLakePoint))

                    iCen = CurrentLakePoint%I
                    jCen = CurrentLakePoint%J

                    do j = -1, 1
                    do i = -1, 1

                        if (.not. (i == 0 .and. j == 0)) then
                    
                            iAdj = iCen + i
                            jAdj = jCen + j
                            
                            if (Me%LakePoints(iAdj, jAdj) == CurrentLake%LakeID) then
                            if (Me%DrainageDirection(iAdj, jAdj) /= NODRAINAGE_ .and.    &
                                Me%DrainageDirection(iCen, jCen) == NODRAINAGE_) then
                                Me%DrainageDirection(iCen, jCen) = ConvertDirection(i, j)

                                dx = XCenter (iCen, jCen) - XCenter (iAdj, jAdj)
                                dy = YCenter (iCen, jCen) - YCenter (iAdj, jAdj)
                                dz = sqrt(dx**2.+dy**2.) * Me%SlopeDepressions

                                Me%NewTopography    (iCen, jCen) = Me%NewTopography(iAdj, jAdj) + dz
                                Me%TopographyChanged = .true.
                            endif
                            endif
                        endif
                    enddo
                    enddo

                    CurrentLakePoint => CurrentLakePoint%Next
                
                enddo
                
            enddo

            CurrentLake => CurrentLake%Next

        enddo

    end subroutine FlatDepressions

    !--------------------------------------------------------------------------

    logical function LakeHasUndrainedPoints (CurrentLake)

        !Arguments-------------------------------------------------------------
        type (T_Lake), pointer                      :: CurrentLake
    
        !Local-----------------------------------------------------------------
        type (T_LakePoint), pointer                 :: CurrentLakePoint

        LakeHasUndrainedPoints = .false.
        CurrentLakePoint => CurrentLake%FirstLakePoint
        do while (associated(CurrentLakePoint))
            if (Me%DrainageDirection(CurrentLakePoint%I, CurrentLakePoint%J) == NODRAINAGE_) then
                LakeHasUndrainedPoints = .true.
                return
            endif
            CurrentLakePoint => CurrentLakePoint%Next
        enddo

    end function LakeHasUndrainedPoints

    !--------------------------------------------------------------------------

    subroutine CalculateCellSlope

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: iCen, jCen
        real                                        :: MajorSlope, Slope
        integer                                     :: Direction, i, j
        real                                        :: CenX, CenY
        real                                        :: AdjX, AdjY
        integer                                     :: iAdj, jAdj
        real                                        :: dist

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


        !Calculates the cell slope (maior slope off all adjecent cells)
        do jCen = JLB, JUB
        do iCen = ILB, IUB

            MajorSlope = null_real
            Direction  = null_int

            !Coordinates of Center Point
            CenX = XCenter (iCen, jCen)
            CenY = YCenter (iCen, jCen)

            do j = -1, 1
            do i = -1, 1

                if (.not. (i == 0 .and. j == 0)) then
                    
                    iAdj = iCen + i
                    jAdj = jCen + j

                    if (Me%BoundaryPoints(iAdj, jAdj) == Not_Boundary) then

                        !Coordinates of adjcent Point
                        AdjX    = XCenter (iAdj, jAdj)
                        AdjY    = YCenter (iAdj, jAdj)

                        !Distance between points
                        dist    = sqrt ((CenX - AdjX)**2. + (CenY - AdjY)**2.)

                        !Slope
                        Slope   = (Me%NewTopography(iCen, jCen) - Me%NewTopography(iAdj, jAdj)) / dist

                    else
                        Slope   = null_real
                    endif
                endif

                if (Slope > MajorSlope .and. Slope > 0) then
                    MajorSlope               = Slope
                endif

            enddo
            enddo

            Me%CellSlope(iCen, jCen) = MajorSlope

        enddo
        enddo

    end subroutine CalculateCellSlope

    !--------------------------------------------------------------------------

    subroutine DelineateBasin ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: iCen, jCen
        integer                                     :: iTarget, jTarget


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Sets all Border Points as NoBasinPoints
        do jCen = JLB, JUB
        do iCen = ILB, IUB
            if (Me%BoundaryPoints(iCen, jCen) == Boundary) then
                Me%BasinPoints(iCen, jCen) = NoBasinPoint
            endif
        enddo
        enddo

        !Sets the outlet as belonging to the basin
        Me%BasinPoints(Me%OutletI, Me%OutletJ) = BasinPoint

        !For each cell looks if the drainage leeds to the outlet
        do jCen = JLB, JUB
        do iCen = ILB, IUB

            call TargetPoint(Me%DrainageDirection(iCen, jCen), iCen, jCen, iTarget, jTarget)

do1:        do 

                !Target Point is the outlet?
                if (iTarget == Me%OutletI .and. jTarget == Me%OutletJ) then
                    Me%BasinPoints(iCen, jCen) = BasinPoint
                    exit do1
                endif

                !If the target knows where it drains to, the current will drain to the same place
                if (Me%BasinPoints(iTarget, jTarget) /= UndefinedPoint) then
                    Me%BasinPoints(iCen, jCen) = Me%BasinPoints(iTarget, jTarget)
                    exit do1
                endif

                !Target is Boundary => Outside Basin
                if (Me%BoundaryPoints(iTarget, jTarget) == Boundary) then
                    Me%BasinPoints(iCen, jCen) = NoBasinPoint
                    exit do1
                endif                    

                !New Target
                call TargetPoint (Me%DrainageDirection(iTarget, jTarget), iTarget, jTarget, iTarget, jTarget)

            enddo do1

        enddo
        enddo

        !Sets the outlet as belonging to the basin - must be repeated here
        Me%BasinPoints(Me%OutletI, Me%OutletJ) = BasinPoint

        !Removes all River Points which lies outside the delineation
        do jCen = JLB, JUB
        do iCen = ILB, IUB
            if (Me%BasinPoints(iCen, jCen) == NoBasinPoint) then
                Me%RiverPoints(iCen, jCen)  = NoRiverPoint
            endif
        enddo
        enddo


    end subroutine DelineateBasin

    !--------------------------------------------------------------------------

    subroutine SetInteriorBasinPoints

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, i, j


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        Me%BasinPoints = NoBasinPoint

        !Sets all Border Points as NoBasinPoints
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%BoundaryPoints(i, j) == Not_Boundary) then
                Me%BasinPoints(i, j) = BasinPoint
            endif
        enddo
        enddo

    end subroutine SetInteriorBasinPoints

    !--------------------------------------------------------------------------

    subroutine AddLake(NewLake)

        !Arguments-------------------------------------------------------------
        type (T_Lake), pointer                      :: NewLake

        !Local-----------------------------------------------------------------
        type (T_Lake), pointer                      :: PreviousLake
        !type (T_Lake), pointer                      :: CurrentLake


        !Insert New Lake into list
        if (.not. associated(Me%FirstLake)) then
            Me%FirstLake      => NewLake
        else

            PreviousLake => Me%FirstLake
            Me%FirstLake => NewLake
            NewLake%Next => PreviousLake

            !PreviousLake      => Me%FirstLake
            !CurrentLake       => Me%FirstLake%Next
            !do while (associated(CurrentLake))
            !    PreviousLake  => CurrentLake
            !    CurrentLake   => CurrentLake%Next
            !enddo
            !CurrentLake       => NewLake
            !PreviousLake%Next => NewLake
        endif

    end subroutine AddLake

    !--------------------------------------------------------------------------

    subroutine AddLakePoint(Lake, i, j)

        !Arguments-------------------------------------------------------------
        type (T_Lake), pointer                      :: Lake
        integer                                     :: i, j

        !Local-----------------------------------------------------------------
        type (T_LakePoint), pointer                 :: NewLakePoint
        type (T_LakePoint), pointer                 :: CurrentLakePoint
        type (T_LakePoint), pointer                 :: PreviousLakePoint


        !Allocates new Lake Point
        allocate (NewLakePoint)
        nullify  (NewLakePoint%Next)

        NewLakePoint%I = i
        NewLakePoint%J = j

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Lake%FirstLakePoint)) then
            Lake%FirstLakePoint => NewLakePoint
        else
            PreviousLakePoint  => Lake%FirstLakePoint
            CurrentLakePoint   => Lake%FirstLakePoint%Next
            do while (associated(CurrentLakePoint))
                PreviousLakePoint => CurrentLakePoint
                CurrentLakePoint  => CurrentLakePoint%Next
            enddo

            PreviousLakePoint%Next => NewLakePoint
        endif

    end subroutine AddLakePoint

    !--------------------------------------------------------------------------

    subroutine RemoveLake (Lake)

        !Arguments-------------------------------------------------------------
        type (T_Lake), pointer                      :: Lake

        !Local-----------------------------------------------------------------
        type (T_Lake), pointer                      :: PreviousLake
        type (T_Lake), pointer                      :: CurrentLake


        !Updates pointers
        if (Lake%LakeID == Me%FirstLake%LakeID) then
            Me%FirstLake => Me%FirstLake%Next
        else
            PreviousLake => Me%FirstLake
            CurrentLake  => Me%FirstLake%Next
            do while (CurrentLake%LakeID /= Lake%LakeID)
                PreviousLake => CurrentLake
                CurrentLake  => CurrentLake%Next
            enddo

            !Now update linked list
            PreviousLake%Next => CurrentLake%Next

        endif

        call RemoveLakePoints (Lake%FirstLakePoint)
            
        !Deallocate Lake
        deallocate (Lake)
        nullify    (Lake) 

    end subroutine RemoveLake
    
    !--------------------------------------------------------------------------

    subroutine RemoveLakePoints (LakePoint)

        !Arguments-------------------------------------------------------------
        type (T_LakePoint), pointer                 :: LakePoint

        !Local-----------------------------------------------------------------
        type (T_LakePoint), pointer                 :: NextLakePoint

        do while (associated(LakePoint))
            NextLakePoint => LakePoint%Next
            deallocate (LakePoint)
            LakePoint => NextLakePoint
        enddo

    end subroutine RemoveLakePoints

    !--------------------------------------------------------------------------

    subroutine RemoveAllDepressions

        !Local-----------------------------------------------------------------
        type (T_Lake), pointer                      :: CurrentLake


        CurrentLake=> Me%FirstLake

        do while (associated(CurrentLake))
            call RemoveLake (CurrentLake)
            CurrentLake  => Me%FirstLake
        enddo

        Me%nDepressions = 0
        nullify (Me%FirstLake)

    end subroutine RemoveAllDepressions

    !--------------------------------------------------------------------------

    subroutine UpStreamAreas 

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: iCen, jCen
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: iTarget, jTarget
        real                                        :: Area
        type (T_UpstreamAreas), pointer             :: CurrentUpstreamArea

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Area of each cell
        do jCen = JLB, JUB
        do iCen = ILB, IUB
            if (Me%BoundaryPoints (iCen, jCen)  == Not_Boundary) then
                !Area of the current cell
                Me%UpStreamArea(iCen, jCen)  = Me%ExtVar%CellArea(iCen, jCen)
            endif
        enddo
        enddo


        do jCen = JLB, JUB
        do iCen = ILB, IUB

            if (Me%BoundaryPoints (iCen, jCen)  == Not_Boundary) then

                Area = Me%ExtVar%CellArea(iCen, jCen)
    
                !Adds user defined upstream areas
                CurrentUpstreamArea => Me%FirstUpstreamArea
                do while (associated(CurrentUpstreamArea))
                    if (iCen == CurrentUpstreamArea%I .and. jCen == CurrentUpstreamArea%J) then
                        Area = Area + CurrentUpstreamArea%Area
                    endif
                    CurrentUpstreamArea => CurrentUpstreamArea%Next
                enddo

                
                call TargetPoint (Me%DrainageDirection(iCen, jCen), iCen, jCen, iTarget, jTarget)

do1:            do 

                    !Sums Area 
                    Me%UpStreamArea(iTarget, jTarget) = Me%UpStreamArea(iTarget, jTarget) + Area

                    !Target Point is boundary
                    if (Me%BoundaryPoints (iTarget, jTarget) == Boundary) then
                        exit do1
                    endif

                    !New Target
                    call TargetPoint (Me%DrainageDirection(iTarget, jTarget), iTarget, jTarget, iTarget, jTarget)

                enddo do1

            endif

        enddo
        enddo




    end subroutine UpStreamAreas

    !--------------------------------------------------------------------------

    subroutine NodeDefinition

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: iCen, jCen
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: NumberOfNodes
        logical                                     :: IsolatedOutletNode, ConfluenceInOutletNode
        integer                                     :: iTarget, jTarget
        logical                                     :: NodesChanged

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Sets Points equal to NoRiverPoints
        Me%RiverPoints = NoRiverPoint

        !Set Nodes where Drained Area is superior to Treshold Area.
        !Update: for several outlet mode isolated nodes may appear
        !without reach. To account a node it has to have a reach associated. 
        !Also outlets with more than one upstream nodes will not be considered
        !for drainage network consistency - network will stop upstream. David Fev 2013
        
        NumberOfNodes = 0 
        NodesChanged = .false.
        do jCen = JLB, JUB
        do iCen = ILB, IUB
            if (Me%BoundaryPoints (iCen, jCen) == Not_Boundary     .and.    &
                Me%UpStreamArea(iCen, jCen) >= Me%RiverTresholdArea) then
                
                IsolatedOutletNode     = .false.
                ConfluenceInOutletNode = .false.
                
                if (.not. Me%DelineateBasin) then    
                    
                    call TargetPoint (Me%DrainageDirection(iCen, jCen), iCen, jCen, iTarget, jTarget)
                    
                    if (Me%BoundaryPoints(iTarget, jTarget) == Boundary) then
                        !this is an outlet
                        !if outlet is isolated, do not consider it
                        if (NumberOfUpstreamNodes(iCen, jCen) .eq. 0) then
                            IsolatedOutletNode = .true.
                        !if outlet has more than one upstream node, do not consider it
                        elseif (NumberOfUpstreamNodes(iCen, jCen) .gt. 1) then
                            ConfluenceInOutletNode = .true.
                            NodesChanged = .true.                            
                        endif             
                    endif
                endif
                
                if ((.not. IsolatedOutletNode) .and. (.not. ConfluenceInOutletNode)) then
                    NumberOfNodes          = NumberOfNodes + 1
                    Me%NodeIDs(iCen, jCen) = NumberOfNodes
                    Me%RiverPoints (iCen, jCen) = RiverPoint
                endif             
            
            endif
        enddo
        enddo
        
        !if a concluence node was removed maybe the isolated node and confluence node are now upstream
        !everytime a confluence is removed the check has to be redone
        !isolated nodes removal does not take problem upstream since there will be no more nodes
        if (.not. Me%DelineateBasin) then
            do while (NodesChanged)
                NodesChanged = .false.
                do jCen = JLB, JUB
                do iCen = ILB, IUB
                    if (Me%RiverPoints (iCen, jCen) == RiverPoint) then
                        call TargetPoint (Me%DrainageDirection(iCen, jCen), iCen, jCen, iTarget, jTarget)
                        !test the new outlets
                        if (Me%RiverPoints(iTarget, jTarget) == NoRiverPoint) then
                            !removing nodes created a isolated node upstream of boundary
                            if (NumberOfUpstreamNodes(iCen, jCen) .eq. 0) then
                                Me%NodeIDs(iCen, jCen)      = 0
                                Me%RiverPoints (iCen, jCen) = NoRiverPoint
                            !removing nodes created a confluence upstream
                            elseif (NumberOfUpstreamNodes(iCen, jCen) .gt. 1) then
                                Me%NodeIDs(iCen, jCen)      = 0
                                Me%RiverPoints (iCen, jCen) = NoRiverPoint 
                                NodesChanged = .true.                       
                            endif
                        endif
                    endif
                enddo
                enddo
            enddo
        endif
        
    end subroutine NodeDefinition

    !--------------------------------------------------------------------------

    integer function NumberOfUpstreamNodes(iCen, jCen)
        !Arguments-------------------------------------------------------------
        integer                                      :: iCen, jCen
        !Local-----------------------------------------------------------------
        integer                                      :: i, j, iAdj, jAdj, iTarget, jTarget
        
        !go trough the neighbours and count how many point to it and are really nodes (UpstreamArea >= TresholdArea)
        NumberOfUpstreamNodes = 0
        
        do j = -1, 1
        do i = -1, 1

            if (.not. (i == 0 .and. j == 0)) then
                
                iAdj = iCen + i
                jAdj = jCen + j

                if (Me%BoundaryPoints(iAdj, jAdj) == Not_Boundary) then
                    
                    call TargetPoint (Me%DrainageDirection(iAdj, jAdj), iAdj, jAdj, iTarget, jTarget)
                    
                    if ((iTarget .eq. iCen) .and. (jTarget .eq. jCen)                 &
                        .and. (Me%UpstreamArea(iAdj, jAdj) .ge. Me%RiverTresholdArea)) then
                        NumberOfUpstreamNodes = NumberOfUpstreamNodes + 1            
                    endif
                
                endif
                
            endif
        
        enddo
        enddo
                
    end function NumberOfUpstreamNodes
    
    !--------------------------------------------------------------------------
    
    subroutine TargetPoint (Direction, iCen, jCen, TargetI, TargetJ)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Direction
        integer, intent(IN)                         :: iCen,    jCen
        integer, intent(OUT)                        :: TargetI, TargetJ

        !Local-----------------------------------------------------------------

        select case (Direction)

        case (N_ )
            TargetI = iCen + 1
            TargetJ = jCen
        case (NE_)
            TargetI = iCen + 1
            TargetJ = jCen + 1
        case (E_ )
            TargetI = iCen
            TargetJ = jCen + 1
        case (SE_)
            TargetI = iCen - 1
            TargetJ = jCen + 1
        case (S_ )
            TargetI = iCen - 1
            TargetJ = jCen
        case (SW_)
            TargetI = iCen - 1
            TargetJ = jCen - 1
        case (W_ )
            TargetI = iCen
            TargetJ = jCen - 1
        case (NW_)
            TargetI = iCen + 1
            TargetJ = jCen - 1
        case (NODRAINAGE_)
            TargetI = iCen
            TargetJ = jCen
        case default
            stop 'TargetPoint - ModuleBasinGeometry - ERR01' 
        end select

    end subroutine TargetPoint

    !--------------------------------------------------------------------------

    integer function ConvertDirection(i, j)

        !Arguments-------------------------------------------------------------
        integer                                     :: i, j

        !Local-----------------------------------------------------------------


        if      (i ==  1 .and. j ==  0) then
            ConvertDirection = N_

        elseif  (i ==  1 .and. j ==  1) then
            ConvertDirection = NE_

        elseif  (i ==  0 .and. j ==  1) then
            ConvertDirection = E_

        elseif  (i == -1 .and. j ==  1) then
            ConvertDirection = SE_

        elseif  (i == -1 .and. j ==  0) then
            ConvertDirection = S_

        elseif  (i == -1 .and. j == -1) then
            ConvertDirection = SW_

        elseif  (i ==  0 .and. j == -1) then
            ConvertDirection = W_

        elseif  (i ==  1 .and. j == -1) then
            ConvertDirection = NW_

        else
            stop 'ConvertDirection - ModuleBasinGeometry - ERR01' 
        endif

    end function

    !--------------------------------------------------------------------------

    real function XCenter (i, j)

        !Arguments-------------------------------------------------------------
        integer                                     :: i, j

        !Local-----------------------------------------------------------------

        !XCenter = (Me%ExtVar%XX_IE(i, j) + Me%ExtVar%XX_IE(i+1, j) + Me%ExtVar%XX_IE(i, j+1) + Me%ExtVar%XX_IE(i+1, j+1)) / 4.0
        XCenter = Me%ExtVar%XX2D_Z(i, j)

    end function XCenter 

    !--------------------------------------------------------------------------

    real function YCenter (i, j)

        !Arguments-------------------------------------------------------------
        integer                                     :: i, j

        !Local-----------------------------------------------------------------

        !YCenter = (Me%ExtVar%YY_IE(i, j) + Me%ExtVar%YY_IE(i+1, j) + Me%ExtVar%YY_IE(i, j+1) + Me%ExtVar%YY_IE(i+1, j+1)) / 4.0
        YCenter = Me%ExtVar%YY2D_Z(i, j)


    end function YCenter 

    !--------------------------------------------------------------------------

    subroutine WriteBasinHDF ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: HDF5_CREATE
        integer                                     :: ObjHDF5
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Reads the name of the Transient File
        call ReadFileName('BASIN_GEOMETRY_HDF', Me%TransientHDF,                                 &
                          Message = "Transient HDF File of Module Basin Geometry",               &
                          Extension = 'hdf5', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR01'


        call GetHDF5FileAccess (HDF5_CREATE = HDF5_CREATE)

        call ConstructHDF5 (ObjHDF5, trim(Me%TransientHDF), HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR02'

        call WriteHorizontalGrid (Me%ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR03'

        call HDF5SetLimits (ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR04'

        call HDF5WriteData (ObjHDF5, "/Grid", "Bathymetry", "m",                         &
                            Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR04a'

        call HDF5WriteData (ObjHDF5, "/Basin", "Upstream Area", "m2",                    &
                            Array2D = Me%UpstreamArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR05'

        call HDF5WriteData (ObjHDF5, "/Basin", "Drainage Direction", "--",               &
                            Array2D = Me%DrainageDirection, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR06'

        call HDF5WriteData (ObjHDF5, "/Basin", "Basin Points", "--",                     &
                            Array2D = Me%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR07'

        call KillHDF5      (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasin - ModuleBasinGeometry - ERR08'

    end subroutine WriteBasinHDF
    
    !--------------------------------------------------------------------------

    subroutine WriteBasinASCII

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, myID
        integer                                     :: ILB, IUB, JLB, JUB, Unit
        integer                                     :: StartI, StartJ, CurrI, CurrJ, LastDirection
        integer                                     :: iCen, jCen, iTarget, jTarget, iNode
        integer                                     :: iOutletTarget, jOutletTarget
        real                                        :: Xorig, Yorig, Angle, x, y
        logical                                     :: FoundDir
        real, dimension(:), pointer                 :: XX,   YY
        real, dimension(:, :), pointer              :: ReservoirIDS
        real                                        :: ReservoirFillValue
        integer                                     :: ObjGridDataReservoirs

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Gets Origin
        call GetGridOrigin(Me%ObjHorizontalGrid, Xorig, Yorig, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR01'

        !Gets Angle
        call GetGridAngle (Me%ObjHorizontalGrid, Angle, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR01a'

        !Gets XX, YY        
        call GetHorizontalGrid     (Me%ObjHorizontalGrid, XX = XX, YY = YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR02'
        

        !Upstream Area
        if (Me%WriteDrainedArea) then
            call WriteGridData  (Me%DrainedAreaFile,                                   &
                                 COMENT1          = "Temporary File",                  &
                                 COMENT2          = "Upstream Area",                   &
                                 HorizontalGridID = Me%ObjHorizontalGrid,              &
                                 FillValue        = -99.0,                             &
                                 OverWrite        = .true.,                            &
                                 GridData2D_Real  = Me%UpstreamArea,                   &
                                 STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR03'
        endif

        !Drainage Direction
        if (Me%WriteDrainageDirection) then
            call WriteGridData  (Me%DrainageDirectionFile,                             &
                                 COMENT1          = "Temporary File",                  &
                                 COMENT2          = "Drainage Direction",              &
                                 HorizontalGridID = Me%ObjHorizontalGrid,              &
                                 FillValue        = 0.0,                               &
                                 OverWrite        = .true.,                            &
                                 GridData2D_Int   = Me%DrainageDirection,              &
                                 STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR04'
        endif

        !Slope
        if (Me%WriteCellSlope) then
            call WriteGridData  (Me%CellSlopeFile,                                     &
                                 COMENT1          = "Temporary File",                  &
                                 COMENT2          = "CellSlope",                       &
                                 HorizontalGridID = Me%ObjHorizontalGrid,              &
                                 FillValue        = 0.0,                               &
                                 OverWrite        = .true.,                            &
                                 GridData2D_Real  = Me%CellSlope,                      &
                                 STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR04'
        endif

        !Starts GridData for Reservoirs, if any
        if (Me%Reservoirs) then
            call ConstructGridData(ObjGridDataReservoirs, Me%ObjHorizontalGrid,          &
                                   FileName = Me%ReservoirFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'NodeDefinition - ModuleBasinGeometry - ERR01'

            call GetGridData (ObjGridDataReservoirs, ReservoirIDS, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'NodeDefinition - ModuleBasinGeometry - ERR02'

            call GetFillValue(ObjGridDataReservoirs, ReservoirFillValue, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'NodeDefinition - ModuleBasinGeometry - ERR03'
        endif

        


        !Writes Reaches File
        if (Me%WriteReaches) then

            !Reclassifies Nodes so NodeIDs are continues from 1 to the number of nodes
            iNode = 1
            Me%NodeIDs = null_int
            do jCen = JLB, JUB
            do iCen = ILB, IUB

                if (Me%RiverPoints      (iCen, jCen)  == RiverPoint) then

                    !Don't write Reaches which upstream node is inside an reservoir
                    if (Me%Reservoirs) then
                        if (ReservoirIDS(iCen , jCen) > ReservoirFillValue) cycle 
                    endif
                  
                    Me%NodeIDs(iCen, jCen) = iNode
                    iNode = iNode + 1

                endif

            enddo
            enddo


            !Target Point of outlet node is a ghost cell that will be written in Drainage Network file
            if (Me%DelineateBasin) then
                call TargetPoint (Me%DrainageDirection(Me%OutletI, Me%OutletJ), Me%OutletI, &
                                  Me%OutletJ, iOutletTarget, jOutletTarget)
                
                Me%NodeIDs (iOutletTarget, jOutletTarget) = iNode
            
            endif


            call UnitsManager (Unit, OPEN_FILE)
            open (unit = unit, file = trim(Me%ReachesFile), status = 'unknown')

            !Writes Node Definitions
            do jCen = JLB, JUB
            do iCen = ILB, IUB

                if (Me%RiverPoints (iCen, jCen)  == RiverPoint) then
                    
                    write (unit, *)"<BeginNode>"
                    write (unit, *)"ID                 : ", Me%NodeIDs(iCen, jCen)

                    x = Me%ExtVar%XX2D_Z(iCen, jCen)
                    y = Me%ExtVar%YY2D_Z(iCen, jCen)
                    !call RODAXY(Xorig, Yorig, Angle, x, y)  

                    write (unit, *)"COORDINATES        : ", x, y
                    write (unit, *)"CROSS_SECTION_TYPE : ", 0   !Not Defined yet 
                    write (unit, *)"DRAINED_AREA       : ", Me%UpStreamArea     (iCen, jCen)
                    write (unit, *)"BOTTOM_LEVEL       : ", Me%ExtVar%Topography(iCen, jCen)
                    write (unit, *)"TERRAIN_LEVEL      : ", Me%ExtVar%Topography(iCen, jCen)
                    write (unit, *)"GRID_I             : ", iCen
                    write (unit, *)"GRID_J             : ", jCen
                    write (unit, *)"<EndNode>"
                
                end if

            enddo
            enddo

            !Write Outlet Ghost Node Definition (it is not a river point)
            if (Me%DelineateBasin) then

                write (unit, *)"<BeginNode>"
                write (unit, *)"ID                 : ", Me%NodeIDs(iOutletTarget, jOutletTarget)


                x = Me%ExtVar%XX2D_Z(iOutletTarget, jOutletTarget)
                y = Me%ExtVar%YY2D_Z(iOutletTarget, jOutletTarget)

!                x = Me%ExtVar%XX_Z(jOutletTarget)
!                y = Me%ExtVar%YY_Z(iOutletTarget)
!                call RODAXY(Xorig, Yorig, Angle, x, y)  

                write (unit, *)"COORDINATES        : ", x, y
                write (unit, *)"CROSS_SECTION_TYPE : ", 0   !Not Defined yet 
                write (unit, *)"DRAINED_AREA       : ", Me%UpStreamArea     (iOutletTarget, jOutletTarget)
                write (unit, *)"BOTTOM_LEVEL       : ", Me%ExtVar%Topography(iOutletTarget, jOutletTarget)
                write (unit, *)"TERRAIN_LEVEL      : ", Me%ExtVar%Topography(iOutletTarget, jOutletTarget)
                write (unit, *)"GRID_I             : ", iOutletTarget
                write (unit, *)"GRID_J             : ", jOutletTarget
                write (unit, *)"<EndNode>"           
            
            endif



            !Writes Reaches Definition
            myID = 0
            do jCen = JLB, JUB
            do iCen = ILB, IUB

                if (Me%RiverPoints      (iCen, jCen)  == RiverPoint) then

                    !Don't write reach after outlet
                    !if (Me%DelineateBasin .and. iCen == Me%OutletI .and. jCen == Me%OutletJ) cycle

                    !Don't write Reaches which upstream node is inside an reservoir
                    if (Me%Reservoirs) then
                        if (ReservoirIDS(iCen , jCen) > ReservoirFillValue) cycle 
                    endif

                    call TargetPoint (Me%DrainageDirection(iCen, jCen), iCen, jCen, iTarget, jTarget)
                    
                    if (Me%NodeIDs(iCen, jCen) > 0 .and. Me%NodeIDs(iTarget, jTarget) > 0) then

                        myID = myID + 1
                        write (unit, *)"<BeginReach>"
                        write (unit, *)"ID              : ", myID
                        write (unit, *)"UPSTREAM_NODE   : ", Me%NodeIDs(iCen,    jCen   )
                        write (unit, *)"DOWNSTREAM_NODE : ", Me%NodeIDs(iTarget, jTarget)
                        write (unit, *)"ACTIVE          : ", 1
                        write (unit, *)"<EndReach>"
                    endif

                endif

            enddo
            enddo

            call UnitsManager (Unit, CLOSE_FILE)

        endif


        if (Me%Reservoirs) then
            call UnGetGridData(ObjGridDataReservoirs, ReservoirIDS, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'NodeDefinition - ModuleBasinGeometry - ERR04'

            call KillGridData (ObjGridDataReservoirs, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'NodeDefinition - ModuleBasinGeometry - ERR05'
        endif


        !Writes Basin Delimition
        if (Me%WriteDelineation) then

            call UnitsManager (Unit, OPEN_FILE)
            open (unit = unit, file = trim(Me%DelineationFile), status = 'unknown')

            write (unit, *)"<beginpolygon>"

            !Selects Starting Corner of the polygon
            if      (Me%BasinPoints(Me%OutletI, Me%OutletJ-1) == NoBasinPoint) then
                StartI = Me%OutletI
                StartJ = Me%OutletJ
            elseif  (Me%BasinPoints(Me%OutletI+1, Me%OutletJ) == NoBasinPoint) then
                StartI = Me%OutletI + 1
                StartJ = Me%OutletJ
            elseif  (Me%BasinPoints(Me%OutletI, Me%OutletJ+1) == NoBasinPoint) then
                StartI = Me%OutletI 
                StartJ = Me%OutletJ + 1
            else
                StartI = Me%OutletI 
                StartJ = Me%OutletJ
            endif                


            x = XX(StartJ)
            y = YY(StartI)
            call RODAXY(Xorig, Yorig, Angle, x, y)  
            write (unit, *)x, y

            LastDirection = null_int
            CurrI         = StartI
            CurrJ         = StartJ

            do while (.not. (CurrI == StartI .and. CurrJ == StartJ) .or. LastDirection == null_int)

                FoundDir = .false.

                !Test Moving North
!                if (LastDirection /= S_ .and. .not. FoundDir) then
                if (LastDirection /= S_) then
                    if (Me%BasinPoints(CurrI  , CurrJ-1) + Me%BasinPoints(CurrI  , CurrJ  ) == 1) then
                        LastDirection = N_
                        CurrI         = CurrI + 1
                        FoundDir      = .True.
                    endif
                endif

                !Test Moving West
!                if (LastDirection /= E_ .and. .not. FoundDir) then
                if (LastDirection /= E_ ) then
                    if (Me%BasinPoints(CurrI  , CurrJ-1) + Me%BasinPoints(CurrI-1, CurrJ-1) == 1) then
                        LastDirection = W_
                        CurrJ         = CurrJ - 1
                        FoundDir      = .True.
                    endif
                endif

                !Test Moving South
!                if (LastDirection /= N_ .and. .not. FoundDir) then
                if (LastDirection /= N_ ) then
                    if (Me%BasinPoints(CurrI-1, CurrJ-1) + Me%BasinPoints(CurrI-1, CurrJ  ) == 1) then
                        LastDirection = S_
                        CurrI         = CurrI - 1
                        FoundDir      = .True.
                    endif
                endif

                !Test Moving East
!                if (LastDirection /= W_ .and. .not. FoundDir) then
                if (LastDirection /= W_ ) then
                    if (Me%BasinPoints(CurrI  , CurrJ  ) + Me%BasinPoints(CurrI-1, CurrJ  ) == 1) then
                        LastDirection = E_
                        CurrJ         = CurrJ + 1
                        FoundDir      = .True.
                    endif
                endif

                if (.not. FoundDir) then
                    stop 'WriteBasinASCII - ModuleBasinGeometry - ERR06'
                endif
    
                x = XX(CurrJ)
                y = YY(CurrI)
                call RODAXY(Xorig, Yorig, Angle, x, y)  

                write (unit, *)x, y

            enddo
           

            write (unit, *)"<endpolygon>"
            
            call UnitsManager (Unit, CLOSE_FILE)
            
            if (Me%WriteBasinPoints) then
            
            call WriteGridData  (Me%BasinPointsFile,                                        &
                                     COMENT1          = "Created using BasinDelineator.",       &
                                     COMENT2          = "BasinPoints 2D GridData: ",            &
                                     HorizontalGridID = Me%ObjHorizontalGrid,                   &
                                     FillValue        = -99.0,                                  &
                                     OverWrite        = .true.,                                 &
                                     GridData2D_Int   = Me%BasinPoints,                         &
                                     STAT             = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR06a'            
            
            endif
            
        endif


        call UngetHorizontalGrid    (Me%ObjHorizontalGrid, XX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR07'

        call UngetHorizontalGrid    (Me%ObjHorizontalGrid, YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBasinASCII - ModuleBasinGeometry - ERR08'

    end subroutine WriteBasinASCII

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetBasinPoints (BasinGeometryID, BasinPoints, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BasinGeometryID
        integer, dimension(:, :),  pointer              :: BasinPoints
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBASINGEOMETRY_, Me%InstanceID)

            BasinPoints => Me%BasinPoints

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBasinPoints

    !--------------------------------------------------------------------------

    subroutine GetRiverPoints (BasinGeometryID, RiverPoints, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BasinGeometryID
        integer, dimension(:, :),  pointer              :: RiverPoints
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBASINGEOMETRY_, Me%InstanceID)

            RiverPoints => Me%RiverPoints

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetRiverPoints

   
    !--------------------------------------------------------------------------

    subroutine GetDrainageDirection (BasinGeometryID, DrainageDirection, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BasinGeometryID
        integer, dimension(:, :),  pointer              :: DrainageDirection
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBASINGEOMETRY_, Me%InstanceID)

            DrainageDirection => Me%DrainageDirection

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetDrainageDirection

    !--------------------------------------------------------------------------    

    subroutine GetCellSlope (BasinGeometryID, CellSlope, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BasinGeometryID
        real, dimension(:, :),  pointer                 :: CellSlope
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBASINGEOMETRY_, Me%InstanceID)

            CellSlope => Me%CellSlope

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetCellSlope

    !--------------------------------------------------------------------------

    subroutine GetMicroChannels (BasinGeometryID, MicroChannels, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BasinGeometryID
        real, dimension(:, :),  pointer                 :: MicroChannels
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBASINGEOMETRY_, Me%InstanceID)

            MicroChannels => Me%MicroChannels

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetMicroChannels
        
    !--------------------------------------------------------------------------

    subroutine UnGetBasin2D_I(BasinGeometryID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BasinGeometryID
        integer, dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBASINGEOMETRY_, Me%InstanceID, "UnGetBasin2D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBasin2D_I


    !--------------------------------------------------------------------------

    subroutine UnGetBasin2D_R4(BasinGeometryID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BasinGeometryID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBASINGEOMETRY_, Me%InstanceID, "UnGetBasin2D_R4")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBasin2D_R4
    

    !--------------------------------------------------------------------------

    subroutine UnGetBasin2D_R8(BasinGeometryID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BasinGeometryID
        real(8), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBASINGEOMETRY_, Me%InstanceID, "UnGetBasin2D_R4")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBasin2D_R8
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillBasinGeometry(BasinGeometryID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: BasinGeometryID
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_
        integer                                     :: nUsers
        type (T_Lake), pointer                      :: CurrentLake

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BasinGeometryID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mBASINGEOMETRY_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deassociates External Instances
                nUsers = DeassociateInstance (mGRIDDATA_,  Me%ObjTopography)
                if (nUsers == 0) stop 'KillBasinGeometry - BasinGeometry - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillBasinGeometry - BasinGeometry - ERR03'

                deallocate (Me%BasinPoints      )
                deallocate (Me%DrainageDirection)
                deallocate (Me%UpStreamArea     )
                deallocate (Me%NodeIDs          )
                deallocate (Me%LakePoints       )

                CurrentLake => Me%FirstLake
                do while (associated(CurrentLake))
                    call RemoveLake(CurrentLake)
                    CurrentLake => Me%FirstLake
                enddo

                call DeallocateInstance

                BasinGeometryID = 0
                STAT_           = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
              
        !----------------------------------------------------------------------

    end subroutine KillBasinGeometry

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_BasinGeometry), pointer          :: AuxObjBasinGeometry
        type (T_BasinGeometry), pointer          :: PreviousObjBasinGeometry

        !Updates pointers
        if (Me%InstanceID == FirstBasinGeometry%InstanceID) then
            FirstBasinGeometry => FirstBasinGeometry%Next
        else
            PreviousObjBasinGeometry => FirstBasinGeometry
            AuxObjBasinGeometry      => FirstBasinGeometry%Next
            do while (AuxObjBasinGeometry%InstanceID /= Me%InstanceID)
                PreviousObjBasinGeometry => AuxObjBasinGeometry
                AuxObjBasinGeometry      => AuxObjBasinGeometry%Next
            enddo

            !Now update linked list
            PreviousObjBasinGeometry%Next => AuxObjBasinGeometry%Next

        endif
            
        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

    end subroutine DeallocateInstance


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetExternalVar () 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Gets Topography
        call GetGridData           (Me%ObjTopography,  Me%ExtVar%Topography,        &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetExternalVar - ModuleBasinGeometry - ERR01'

        !Gets Topography Fill Value
        call GetFillValue          (Me%ObjTopography,  Me%ExtVar%GridFillValue,     &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetExternalVar - ModuleBasinGeometry - ERR01a'


        !Gets Grid
        call GetHorizontalGrid     (Me%ObjHorizontalGrid, XX2D_Z = Me%ExtVar%XX2D_Z,  &
                                    YY2D_Z = Me%ExtVar%YY2D_Z,                        &
                                    XX_Z = Me%ExtVar%XX_Z,                          &
                                    YY_Z = Me%ExtVar%YY_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetExternalVar - ModuleBasinGeometry - ERR02'

        call GetGridCellArea       (Me%ObjHorizontalGrid, Me%ExtVar%CellArea,       &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetExternalVar - ModuleBasinGeometry - ERR03'


    end subroutine GetExternalVar 

    !--------------------------------------------------------------------------

    subroutine UnGetExternalVar () 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        call UngetGridData          (Me%ObjTopography, Me%ExtVar%Topography,        &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UnGetExternalVar - ModuleBasinGeometry - ERR01'

        call UngetHorizontalGrid    (Me%ObjHorizontalGrid, Me%ExtVar%XX2D_Z,         &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UnGetExternalVar - ModuleBasinGeometry - ERR02'

        call UngetHorizontalGrid    (Me%ObjHorizontalGrid, Me%ExtVar%YY2D_Z,         &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UnGetExternalVar - ModuleBasinGeometry - ERR03'

        call UngetHorizontalGrid    (Me%ObjHorizontalGrid, Me%ExtVar%XX_Z,          &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UnGetExternalVar - ModuleBasinGeometry - ERR04'

        call UngetHorizontalGrid    (Me%ObjHorizontalGrid, Me%ExtVar%YY_Z,          &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UnGetExternalVar - ModuleBasinGeometry - ERR05'

        call UngetHorizontalGrid    (Me%ObjHorizontalGrid, Me%ExtVar%CellArea,      &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UnGetExternalVar - ModuleBasinGeometry - ERR06'


    end subroutine UnGetExternalVar 

    !--------------------------------------------------------------------------

    subroutine Ready (ObjBasinGeometry_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBasinGeometry_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjBasinGeometry_ID > 0) then
            call LocateObjBasinGeometry(ObjBasinGeometry_ID)
            ready_ = VerifyReadLock (mBASINGEOMETRY_,  Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjBasinGeometry (ObjBasinGeometryID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBasinGeometryID

        !Local-----------------------------------------------------------------

        Me => FirstBasinGeometry
        do while (associated (Me))
            if (Me%InstanceID == ObjBasinGeometryID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleBasinGeometry - LocateObjBasinGeometry - ERR01'

    end subroutine LocateObjBasinGeometry

    !--------------------------------------------------------------------------


end module ModuleBasinGeometry

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
