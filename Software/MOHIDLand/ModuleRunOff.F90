!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : RunOff
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jan 2004
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module which calculates the Surface RunOff
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

Module ModuleRunOff

    use ModuleGlobalData
    use ModuleTime
    use ModuleTimeSerie         ,only : StartTimeSerieInput, KillTimeSerie,              &
                                        GetTimeSerieInitialData, GetTimeSerieValue,      &
                                        StartTimeSerie, WriteTimeSerieLine,              &
                                        GetNumberOfTimeSeries, GetTimeSerieLocation,     &
                                        TryIgnoreTimeSerie, CorrectsCellsTimeSerie,      &
                                        GetTimeSerieName, WriteTimeSerie
    use ModuleEnterData
    use ModuleHDF5
    use ModuleFunctions         ,only : TimeToString, SetMatrixValue, ChangeSuffix,      &
                                        CHUNK_J, LinearInterpolation, InterpolateValueInTime
    use ModuleHorizontalGrid    ,only : GetHorizontalGridSize, GetHorizontalGrid,        &
                                        UnGetHorizontalGrid, WriteHorizontalGrid,        &
                                        GetGridCellArea, GetXYCellZ,                     &
                                        GetXYInsideDomain,                               &
                                        GetCellZInterceptByLine,                         &
                                        GetCellZInterceptByPolygon, GetGridRotation,     &
                                        GetGridAngle, GetCheckDistortion,                & 
                                        GetCellIDfromIJ, GetCellIJfromID, GetCellZ_XY
    use ModuleHorizontalMap     ,only : GetBoundaries, UngetHorizontalMap
    use ModuleGridData          ,only : GetGridData, UngetGridData, WriteGridData
    use ModuleBasinGeometry     ,only : GetBasinPoints, GetRiverPoints, GetCellSlope,    &
                                        GetDrainageDirection, TargetPoint,               &
                                        UnGetBasin
    use ModuleStopWatch         ,only : StartWatch, StopWatch
    use ModuleFillMatrix        ,only : ConstructFillMatrix, ModifyFillMatrix,           &
                                        KillFillMatrix
    use ModuleDrainageNetwork   ,only : GetChannelsWaterLevel, GetChannelsSurfaceWidth,  &
                                        GetChannelsBankSlope, GetChannelsNodeLength,     &
                                        GetChannelsBottomLevel, UnGetDrainageNetwork,    &
                                        GetChannelsID,GetChannelsVolume,                 &
                                        GetChannelsMaxVolume, GetChannelsActiveState,    &
                                        GetChannelsTopArea, GetChannelsVelocity,         &
                                        GetHasTwoGridPoints
    use ModuleDischarges        ,only : Construct_Discharges, GetDischargesNumber,       &
                                        GetDischargesGridLocalization,                   &
                                        GetDischargeWaterFlow, GetDischargesIDName,      &
                                        TryIgnoreDischarge, GetDischargeSpatialEmission, &
                                        CorrectsCellsDischarges, Kill_Discharges,        &
                                        GetByPassON, GetDischargeFlowDistribuiton,       &
                                        UnGetDischarges, SetLocationCellsZ,              &
                                        CorrectsBypassCellsDischarges
    use ModuleBoxDif,           only :  StartBoxDif, GetBoxes, GetNumberOfBoxes, UngetBoxDif, &
                                        BoxDif, KillBoxDif                                                      
    use ModuleDrawing
    use iso_c_binding
    
    implicit none

    private 

#ifdef _SEWERGEMSENGINECOUPLER_

    interface

        integer(c_int) function SewerGEMSEngine_open(inFile, rptFile, outFile) bind(C, name='swmm_open')
            use iso_c_binding
            character(kind = c_char) :: inFile(*)
            character(kind = c_char) :: rptFile(*)
            character(kind = c_char) :: outFile(*)
        end function SewerGEMSEngine_open
        
       integer(c_int) function SewerGEMSEngine_start(saveResults) bind(C, name='swmm_start')
            use iso_c_binding
            integer(c_int) :: saveResults
        end function SewerGEMSEngine_start
        
        integer(c_int) function SewerGEMSEngine_end() bind(C, name='swmm_end')
            use iso_c_binding
        end function SewerGEMSEngine_end

        integer(c_int) function SewerGEMSEngine_close() bind(C, name='swmm_close')
            use iso_c_binding
        end function SewerGEMSEngine_close

        integer(c_int) function SewerGEMSEngine_getNumberOfNodes(nNodes) bind(C, name='swmm_getNumberOfNodes')
            use iso_c_binding
            integer(c_int) :: nNodes
        end function SewerGEMSEngine_getNumberOfNodes
        
        integer(c_int) function SewerGEMSEngine_getNodeXY(id, xx, yy) bind(C, name='swmm_getNodeXY')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: xx
            real(c_double) :: yy
        end function SewerGEMSEngine_getNodeXY
        
        integer(c_int) function SewerGEMSEngine_getNodeName(id, nName) bind(C, name='swmm_getNodeName')
            use iso_c_binding
            integer(c_int) :: id
            character(kind = c_char) :: nName(*)
        end function SewerGEMSEngine_getNodeName

        integer(c_int) function SewerGEMSEngine_getdt(Dt) bind(C, name='swmm_getdt')
            use iso_c_binding
            real(c_double) :: Dt
        end function SewerGEMSEngine_getdt
        
        integer(c_int) function SewerGEMSEngine_step_imposed_dt(elapsedTime, imposedDt) bind(C, name='swmm_step_imposed_dt')
            use iso_c_binding
            real(c_double) :: elapsedTime
            real(c_double) :: imposedDt
        end function SewerGEMSEngine_step_imposed_dt

        integer(c_int) function SewerGEMSEngine_getOutfallFlow(id, outflow) bind(C, name='swmm_getOutfallFlow')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: outflow
        end function SewerGEMSEngine_getOutfallFlow

        integer(c_int) function SewerGEMSEngine_setOutfallWaterLevel(id, level) bind(C, name='swmm_setOutfallWaterLevel')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: level
        end function SewerGEMSEngine_setOutfallWaterLevel

        integer(c_int) function SewerGEMSEngine_getNodeDownstreamLinkID(id, nodeDownstreamLinkID) bind(C, name='swmm_getNodeDownstreamLinkID')
            use iso_c_binding
            integer(c_int) :: id
            integer(c_int) :: nodeDownstreamLinkID
        end function SewerGEMSEngine_getNodeDownstreamLinkID

        integer(c_int) function SewerGEMSEngine_setHeadwallWaterDepth(id, depth) bind(C, name='swmm_setHeadwallWaterDepth')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: depth
        end function SewerGEMSEngine_setHeadwallWaterDepth

        integer(c_int) function SewerGEMSEngine_getHeadwallDownstreamFlow(id, flow) bind(C, name='swmm_getHeadwallDownstreamFlow')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: flow
        end function SewerGEMSEngine_getHeadwallDownstreamFlow


        integer(c_int) function SewerGEMSEngine_getNodeWaterLevel(id, level) bind(C, name='swmm_getNodeWaterLevel')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: level
        end function SewerGEMSEngine_getNodeWaterLevel

        integer(c_int) function SewerGEMSEngine_getNodeWaterDepth(id, depth) bind(C, name='swmm_getNodeWaterDepth')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: depth
        end function SewerGEMSEngine_getNodeWaterDepth

        integer(c_int) function SewerGEMSEngine_setNodeSurfaceDepth(id, level) bind(C, name='swmm_setNodeSurfaceDepth')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: level
        end function SewerGEMSEngine_setNodeSurfaceDepth

        integer(c_int) function SewerGEMSEngine_setSurfaceLinkFlow(id, inflow) bind(C, name='swmm_setSurfaceLinkFlow')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: inflow
        end function SewerGEMSEngine_setSurfaceLinkFlow

        integer(c_int) function SewerGEMSEngine_getDates(startDate, startTime, endDate, endTime) bind(C, name='swmm_getTimeLimits')
            use iso_c_binding
            character(kind = c_char) :: startDate(*)
            character(kind = c_char) :: startTime(*)
            character(kind = c_char) :: endDate(*)
            character(kind = c_char) :: endTime(*)
        end function SewerGEMSEngine_getDates

        integer(c_int) function SewerGEMSEngine_getNodeType(id, nType) bind(C, name='swmm_getNodeType')
            use iso_c_binding
            integer(c_int) :: id
            integer(c_int) :: nType
        end function SewerGEMSEngine_getNodeType

        integer(c_int) function  SewerGEMSEngine_getIsNodeOpenChannel(id, isOpen) bind(C, name='swmm_getIsNodeOpenChannel')
            use iso_c_binding
            integer(c_int) :: id
            integer(c_int) :: isOpen
        end function SewerGEMSEngine_getIsNodeOpenChannel

        subroutine ConvertSewerGEMSEngineToDrainageNetwork(f1, f2, f3, f4, f5) bind(C, name='ConvertSwmmToDrainageNetworkCaller')
            use iso_c_binding
            character(kind = c_char) :: f1(*)
            character(kind = c_char) :: f2(*)
            character(kind = c_char) :: f3(*)
            character(kind = c_char) :: f4(*)
            character(kind = c_char) :: f5(*)
        end subroutine ConvertSewerGEMSEngineToDrainageNetwork

        integer(c_int) function SewerGEMSEngine_setInletPotentialFlow(id, potentialflow) bind(C, name='swmm_setInletPotentialFlow')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: potentialflow
        end function SewerGEMSEngine_setInletPotentialFlow

        integer(c_int) function  SewerGEMSEngine_getIsNodeSurcharged(id, isOpen) bind(C, name='swmm_isNodeSurcharged')
            use iso_c_binding
            integer(c_int) :: id
            integer(c_int) :: isOpen
        end function SewerGEMSEngine_getIsNodeSurcharged

        integer(c_int) function SewerGEMSEngine_getNodeOverflow(id, flow) bind(C, name='swmm_getNodeOverflow')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: flow
        end function SewerGEMSEngine_getNodeOverflow

        integer(c_int) function SewerGEMSEngine_setNodeSurchargeDepth(id, level) bind(C, name='swmm_setNodeSurchargeDepth')
            use iso_c_binding
            integer(c_int) :: id
            real(c_double) :: level
        end function SewerGEMSEngine_setNodeSurchargeDepth

        integer(c_int) function SewerGEMSEngine_getTotalVolume(totalVolume) bind(C, name='swmm_getTotalVolume')
            use iso_c_binding
            real(c_double) :: totalVolume
        end function SewerGEMSEngine_getTotalVolume
        
    end interface

#endif _SEWERGEMSENGINECOUPLER_

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  ::  ConstructRunOff
    private ::      AllocateInstance
    private ::      AllocateVariables
    private ::      ReadDataFile
    private ::      InitializeVariables
    private ::      ConstructOverLandCoefficient
    private ::      ConstructSewerStormWater
    private ::      ConstructHDF5Output
    private ::      ConstructTimeSeries
    private ::      AllocateRiverGridPointArrays

    !Selector
    public  ::  GetOverLandFlow
    public  ::  GetManning
    public  ::  GetManningDelta
    public  ::  GetFlowToChannels
    public  ::  GetBoundaryImposed
    public  ::  GetRouteDFour
    public  ::  GetRouteDFourCells
    public  ::  GetRouteDFourNeighbours
    public  ::  GetRouteDFourFlux
    public  ::  GetBoundaryFlux
    public  ::  GetBoundaryCells
    public  ::  GetFlowDischarge
    public  ::  GetRunOffTotalDischargeFlowVolume
    public  ::  GetRunoffWaterLevel 
    public  ::  GetRunoffWaterColumn        !Final WaterColumn 
    public  ::  GetRunoffWaterColumnOld     !Initial WaterColumn
    public  ::  GetRunoffWaterColumnAT      !WaterColumn After Transport (For RP) 
    public  ::  GetRunoffCenterVelocity
    public  ::  GetRunoffTotalStoredVolume
    public  ::  GetRunOffStoredVolumes
    public  ::  GetRunOffBoundaryFlowVolume
    public  ::  GetMassError
    public  ::  GetNextRunOffDT
    public  ::  SetBasinColumnToRunoff
    public  ::  UnGetRunOff
    

    !Modifier
    public  ::  ModifyRunOff
    private ::      LocalWaterColumn
    private ::      IntegrateFlow
    private ::      ComputeNextDT
    private ::      RunOffOutput
    private ::      OutputTimeSeries
    private ::  AdjustSlope

    !Destructor
    public  ::  KillRunOff                                                     
    public  ::      SetBasinStatsToRunOff

    !Management
    private ::  ReadLockExternalVar
    private ::  ReadUnLockExternalVar
    private ::  Ready
    private ::      LocateObjRunOff 

    !Interfaces----------------------------------------------------------------
    private :: UnGetRunOff2D_R4
    private :: UnGetRunOff2D_R8
    interface  UnGetRunOff
        module procedure UnGetRunOff2D_R4
        module procedure UnGetRunOff2D_R8
    end interface  UnGetRunOff
    
    !Parameters----------------------------------------------------------------
    integer, parameter                              :: KinematicWave_   = 1
    integer, parameter                              :: DiffusionWave_   = 2
    integer, parameter                              :: DynamicWave_     = 3      

    integer, parameter                              :: UnitMax          = 80
    
    !water column computation in faces
    integer, parameter                              :: WCMaxBottom_     = 1
    integer, parameter                              :: WCAverageBottom_ = 2
    
    !Boundary flux
    integer, parameter                              :: ComputeFlow_       = 1
    integer, parameter                              :: InstantaneousFlow_ = 2

    !Route D4 flux
    integer, parameter                              :: Celerity_          = 1
    integer, parameter                              :: Manning_           = 2    
    
    !Restart fiels format
    integer, parameter                              :: BIN_                 = 1
    integer, parameter                              :: HDF_                 = 2        
    
    !Inlet types
    integer, parameter                              :: Weir_                = 1
    integer, parameter                              :: FlowCapture_         = 2  
    integer, parameter                              :: DepthFlowRatingCurve_= 3    
    integer, parameter                              :: FlowFlowRatingCurve_ = 4    

    !Open channel link types
    integer, parameter                              :: Direct_              = 1
    integer, parameter                              :: Weighted_            = 2  
    integer, parameter                              :: OutfallLink_         = 3  
    integer, parameter                              :: PondLink_            = 4  

    
    !SewerGEMS node types
    !Outside domai or, not connected to 2D RunOff bolted e.g. manhole
    integer, parameter                              :: NotCoupled_          = 0
    !Manhole (can flow from SewerGEMS to 2D RunOff)
    integer, parameter                              :: Manhole_             = 1
    !Inlet/Catch Basin  (can flow from SewerGEMS to 2D RunOff and from 2D to SewerGEMS)
    integer, parameter                              :: Inlet_               = 2
    !Cross Section (can flow from SewerGEMS to 2D RunOff and from 2D to SewerGEMS)
    integer, parameter                              :: CrossSection_        = 3
    !Outfall (can flow from SewerGEMS to 2D RunOff and from 2D to SewerGEMS)
    integer, parameter                              :: Outfall_             = 4
    !Pond (can flow from SewerGEMS to 2D RunOff and from 2D to SewerGEMS)
    integer, parameter                              :: Pond_                = 5
    !Headwall (can flow from SewerGEMS to 2D RunOff and from 2D to SewerGEMS)
    integer, parameter                              :: Headwall_            = 6

    !SWMM node types
    integer, parameter                              :: SWMMJunction_        = 0
    integer, parameter                              :: SWMMOutfall_         = 1
    integer, parameter                              :: SWMMStorage_         = 2
    
    !Types---------------------------------------------------------------------
    
    !TODO: Use inheritance to make these types less code management and repititions
    !GridPoint on top of river 1D node (will recieve water level from 1D model)
    type T_NodeGridPoint
        integer                                     :: ID                   = null_int
        integer                                     :: GridI                = null_int
        integer                                     :: GridJ                = null_int
        real                                        :: RiverLevel           = null_real
        type(T_NodeGridPoint), pointer              :: Next                 => null()
        type(T_NodeGridPoint), pointer              :: Prev                 => null()        
    end type
    type :: NodeGridPointPtr_class                   ! because foooooortraaaaaaan doesn't allow dynamic pointer arrays
        class(T_NodeGridPoint), pointer :: ptr => null() ! the actual pointer
    end type NodeGridPointPtr_class
    !GridPoint on top of river right and left banks (will recieve water level from associated NodeGridPoint 
    !and will be the cells where is computed river interaction flow)
    type T_BankGridPoint
        integer                                     :: ID                   = null_int
        integer                                     :: GridI                = null_int
        integer                                     :: GridJ                = null_int
        integer                                     :: NGPId                = null_int  !associated NodeGridPoint
        integer                                     :: NGPIDidx             = null_int
        real                                        :: RiverLevel           = null_real
        type(T_BankGridPoint), pointer              :: Next                 => null()
        type(T_BankGridPoint), pointer              :: Prev                 => null()           
    end type
    type :: BankGridPointPtr_class                   ! because foooooortraaaaaaan doesn't allow dynamic pointer arrays
        class(T_BankGridPoint), pointer :: ptr => null() ! the actual pointer
    end type BankGridPointPtr_class
    !GridPoint on margins (will recieve water level interpolated from associated BankGridPoint's
    !and will integrate their river interaction flow into closest BGP - based on interpolation X)
    type T_MarginGridPoint
        integer                                     :: ID                   = null_int
        integer                                     :: GridI                = null_int
        integer                                     :: GridJ                = null_int
        integer                                     :: BGPUpId              = null_int  !associated BankGridPoint upstream
        integer                                     :: BGPDownId            = null_int  !associated BankGridPoint downstream        
        integer                                     :: BGPUpIdidx           = null_int
        integer                                     :: BGPDownIdidx         = null_int
        real                                        :: InterpolationFraction  = null_real !x fraction from BGP upstream to downstream segment (at cell center)
        integer                                     :: BGPIntegrateFluxId   = null_int  !associated BGP where to route computed flow (BGP upstream or downstream)
        integer                                     :: NGPIntegrateFluxId   = null_int  !associated NGP where to route computed flow (NGP associated to BGP)
        integer                                     :: GridIIntegrateFlux   = null_int  !I where to integrate flux (BGP in case of DN or NGP in case OpenMI)
        integer                                     :: GridJIntegrateFlux   = null_int  !J where to integrate flux (BGP in case of DN or NGP in case OpenMI)
        real                                        :: RiverLevel           = null_real
        type(T_MarginGridPoint), pointer            :: Next                 => null()
        type(T_MarginGridPoint), pointer            :: Prev                 => null()           
    end type
    type :: MarginGridPointPtr_class                   ! because foooooortraaaaaaan doesn't allow dynamic pointer arrays
        class(T_MarginGridPoint), pointer :: ptr => null() ! the actual pointer
    end type MarginGridPointPtr_class

    type T_SewerGEMSOpenChannelLink
        integer                                     :: ID                   = null_int !1 to number of open channel links
        integer                                     :: I
        integer                                     :: J
        integer                                     :: TypeOf               = null_int  !1 - Direct_ 2 - Weighted_ 3-OutfallLink_ 4-PondLink_
        character(Stringlength)                     :: LinkNodeName         = null_str
        character(Stringlength)                     :: SecondLinkNodeName   = null_str
        integer                                     :: LinkID               = null_int
        integer                                     :: SecondLinkID         = null_int  !if water level is interpolated between 2 nodes
        real                                        :: Weight               = 1.0       !water level interpolation fraction (0 < W < 1)
        integer                                     :: CrossSectionID       = null_int  !internal ID of active cross sections
        integer                                     :: OutfallID            = null_int  !internal ID of active outfalls
        integer                                     :: PondID               = null_int  
        integer                                     :: HeadwallID           = null_int  
        real                                        :: WaterLevel           = null_real
        real                                        :: Flow                 = null_real
        real                                        :: FluxWidth            = null_real
        real                                        :: CellWidth            = null_real
    end type T_SewerGEMSOpenChannelLink

    type T_SewerGEMSPond
        integer                                     :: ID                   = null_int !1 to ponds
        integer                                     :: SWMM_ID              = null_int !SWMM node id
        character(Stringlength)                     :: Name                 = null_str
        real                                        :: WaterLevel           = null_real
        real                                        :: Flow                 = null_real
        integer, dimension(:), allocatable          :: I 
        integer, dimension(:), allocatable          :: J 
        character(Pathlength)                       :: LineFileName         = null_str
        character(Pathlength)                       :: CellsFileName        = null_str
        type (T_Lines),   pointer                   :: Line                 => null()
        integer                                     :: nCells               = null_int 
    end type T_SewerGEMSPond

    type T_SewerGEMSCrossSection
        integer                                     :: ID                   = null_int !1 to number of cross-sections
        integer                                     :: SWMM_ID              = null_int !SWMM node id
        character(Stringlength)                     :: Name                 = null_str
        integer                                     :: I                    = null_int
        integer                                     :: J                    = null_int
        real                                        :: WaterLevel           = null_real
        real                                        :: Flow                 = null_real
    end type T_SewerGEMSCrossSection

    type T_SewerGEMSOutfall
        integer                                     :: ID                   = null_int !1 to number of outfalls
        integer                                     :: SWMM_ID              = null_int !SWMM node id
        character(Stringlength)                     :: Name                 = null_str
        integer                                     :: I                    = null_int
        integer                                     :: J                    = null_int
        real                                        :: WaterLevel           = null_real
        real                                        :: Flow                 = null_real
    end type T_SewerGEMSOutfall
    
    type T_SewerGEMSManhole
        integer                                     :: ID                   = null_int !1 to number of manholes
        integer                                     :: SWMM_ID              = null_int !SWMM node id
        character(Stringlength)                     :: Name                 = null_str
        integer                                     :: I                    = null_int
        integer                                     :: J                    = null_int
        real                                        :: Outflow              = null_real
    end type T_SewerGEMSManhole

    type T_SewerGEMSInlet
        integer                                     :: ID                   = null_int !1 to number of inlets
        integer                                     :: SWMM_ID              = null_int !SWMM node id
        character(Stringlength)                     :: Name                 = null_str
        integer                                     :: I                    = null_int
        integer                                     :: J                    = null_int
        integer                                     :: TypeOf               = null_int
        real                                        :: Width                = null_real
        real                                        :: CaptureFraction      = null_real
        real                                        :: FlowEnteringCell     = 0.0
        real                                        :: PotentialFlow        = 0.0
        real                                        :: EffectiveFlow        = 0.0
        real                                        :: WaterLevel           = 0.0
        integer                                     :: CellID               = null_int
        integer                                     :: nInletsInGridCell    = null_int
        character(Pathlength)                       :: RatingCurveFileName  = null_str
        real, dimension(:), allocatable             :: RatingCurveStage
        real, dimension(:), allocatable             :: RatingCurveFlow
        integer                                     :: RatingCurve_nValues  = null_int 
        real                                        :: RatingCurveBelowMin  = null_real 
        real                                        :: RatingCurveAboveMax  = null_real
        logical                                     :: OutputResults        = .false.
        type (T_Time)                               :: NextOutPutTime 
        real                                        :: OutputTimeStep       = null_real 
        real                                        :: OutputTime           = null_real 
        integer                                     :: OutputUnit           = null_int 
    end type T_SewerGEMSInlet

    type T_SewerGEMSHeadwall
        integer                                     :: ID                   = null_int !1 to number of headwalls
        integer                                     :: SWMM_ID              = null_int !SWMM node id
        integer                                     :: SWMM_DownstreamLinkID= null_int !SWMM downstream link id
        character(Stringlength)                     :: Name                 = null_str
        integer                                     :: I                    = null_int
        integer                                     :: J                    = null_int
        real                                        :: Flow                 = 0.0
        real                                        :: FlowEnteringCell     = 0.0
        logical                                     :: OutputResults        = .false.
        type (T_Time)                               :: NextOutPutTime 
        real                                        :: OutputTimeStep       = null_real 
        real                                        :: OutputTime           = null_real 
        integer                                     :: OutputUnit           = null_int 
    end type T_SewerGEMSHeadwall

    type T_OutPutRunOff
        type (T_Time), pointer, dimension(:)        :: OutTime              => null()
        integer                                     :: NextOutPut           = 1
        logical                                     :: Yes                  = .false.
        type (T_Time), dimension(:), pointer        :: RestartOutTime       => null()
        logical                                     :: WriteRestartFile     = .false.
        logical                                     :: RestartOverwrite     = .false.
        integer                                     :: NextRestartOutput    = 1 
        integer                                     :: RestartFormat         = BIN_
        
        logical                                     :: BoxFluxes            = .false.
        logical                                     :: OutputFloodRisk      = .false.
        real                                        :: FloodRiskVelCoef     = null_real
        
        logical                                     :: WriteMaxFlowModulus  = .false.
        character(Pathlength)                       :: MaxFlowModulusFile   = null_str
        real, dimension(:,:), pointer               :: MaxFlowModulus       => null()

        logical                                     :: WriteMaxWaterColumn  = .false.        
        character(Pathlength)                       :: MaxWaterColumnFile   = null_str
        character(Pathlength)                       :: MaxWaterLevelFile    = null_str
        character(Pathlength)                       :: TimeOfMaxWaterColumnFile = null_str
        real, dimension(:,:), pointer               :: MaxWaterColumn       => null()
        real, dimension(:,:), pointer               :: TimeOfMaxWaterColumn => null()        

        logical                                     :: WriteVelocityAtMaxWaterColumn  = .false.        
        character(Pathlength)                       :: VelocityAtMaxWaterColumnFile   = null_str
        real, dimension(:,:), pointer               :: VelocityAtMaxWaterColumn       => null()        

        logical                                     :: WriteMaxFloodRisk              = .false.        
        character(Pathlength)                       :: MaxFloodRiskFile               = null_str
        real, dimension(:,:), pointer               :: MaxFloodRisk                   => null()            
        
        logical                                     :: WriteFloodPeriod               = .false.        
        character(Pathlength)                       :: FloodPeriodFile                = null_str
        real, dimension(:,:), pointer               :: FloodPeriod                    => null()        
        real                                        :: FloodPeriodWaterColumnLimit    = null_real
        character(Pathlength), dimension(:), pointer:: FloodPeriodFiles               => null()
        real, dimension(:), pointer                 :: FloodPeriodWaterColumnLimits   => null()
        real, dimension(:,:,:), pointer             :: FloodPeriods                   => null()
        integer                                     :: nFloodPeriodLimits             = null_int        

        logical                                     :: WriteFloodArrivalTime          = .false.
        character(Pathlength)                       :: FloodArrivalTimeFile           = null_str
        real, dimension(:,:), pointer               :: FloodArrivalTime               => null()
        real                                        :: FloodArrivalWaterColumnLimit   = null_real        
        real                                        :: TotalFloodedArea               = null_real
        real                                        :: MaxTotalFloodedArea            = null_real
        real                                        :: TimeOfMaxTotalFloodedArea      = null_real

        logical                                     :: TimeSeries                     = .false.
        logical                                     :: TimeSerieDischON               = .false. 
        integer                                     :: DischargesNumber               = null_int 
        integer, dimension(:),   pointer            :: TimeSerieDischID               => null()     
        real,    dimension(:,:), pointer            :: TimeSerieDischProp             => null()
        integer                                     :: TS_Numb_DischProp              = null_int 
        type (T_Time)                               :: NextOutPutDisch    
        real                                        :: OutPutDischDT
        character(len=PathLength)                   :: TimeSerieLocationFile, DiscTimeSerieLocationFile

        
    end type T_OutPutRunOff


    type T_FilesRunOff
        character(PathLength)                       :: DataFile             = null_str
        character(PathLength)                       :: InitialFile          = null_str
        character(PathLength)                       :: FinalFile            = null_str
        character(PathLength)                       :: TransientHDF         = null_str
        character(PathLength)                       :: BoxesFile            = null_str
        character(PathLength)                       :: SWMMdat              = null_str    
        character(PathLength)                       :: SWMMinp              = null_str    
        character(PathLength)                       :: SWMMrpt              = null_str
        character(PathLength)                       :: SWMMout              = null_str
        character(PathLength)                       :: SWMMUncoupledElements= null_str
        character(PathLength)                       :: SWMMHDF              = null_str
        character(PathLength)                       :: SWMMTimeSeries       = null_str
        character(PathLength)                       :: SWMMTimeSeriesDir    = null_str
        character(PathLength)                       :: MassErrorFile        = null_str
        character(PathLength)                       :: RunOffLogFile        = null_str
    end type T_FilesRunOff    

    type T_ExtVarRunOff
        integer, dimension(:,:), pointer            :: BasinPoints              => null()
        real(8), dimension(:,:), pointer            :: WaterColumn              => null()
        real   , dimension(:,:), pointer            :: GridCellArea             => null()
        real   , dimension(:,:), pointer            :: DUX, DUY                 => null()
        real   , dimension(:,:), pointer            :: DVX, DVY                 => null()
        real   , dimension(:,:), pointer            :: DXX, DYY                 => null()
        real   , dimension(:,:), pointer            :: DZX, DZY                 => null()
        real   , dimension(:,:), pointer            :: XX2D_Z, YY2D_Z           => null()
        real   , dimension(:,:), pointer            :: Topography               => null()
        integer, dimension(:,:), pointer            :: BoundaryPoints2D         => null()
        integer, dimension(:,:), pointer            :: RiverPoints              => null()
        real   , dimension(:,:), pointer            :: CellSlope                => null()
        logical                                     :: Distortion               = .false.
        real   , dimension(:,:), pointer            :: RotationX, RotationY     => null()
        real                                        :: GridRotation             = null_real
        type (T_Time)                               :: Now
        real                                        :: DT                       = null_real
    end type T_ExtVarRunOff

    type T_Converge
        integer                                     :: MinIterations                = 1               
        integer                                     :: MaxIterations                = 1024
        logical                                     :: Stabilize                    = .false.
        real                                        :: StabilizeFactor              = 0.01        
        real                                        :: DTFactorUp                   = 1.25
        real                                        :: DTFactorDown                 = 1.25
        real                                        :: StabilizeHardCutLimit        = 128
        real                                        :: DTSplitFactor                = 2.0               
        real                                        :: CurrentDT                    = null_real  
        real                                        :: NextDT                       = null_real
        integer                                     :: LastGoodNiteration           = 1
        integer                                     :: NextNiteration               = 1               
        logical                                     :: LimitDTCourant               = .false.        
        real                                        :: MaxCourant                   = 1.0  
        integer                                     :: MinToRestart                 = 0  
        real                                        :: MinimumValueToStabilize      = 0.001
        logical                                     :: CheckDecreaseOnly            = .false.        
        logical                                     :: CorrectDischarge             = .false.
        logical                                     :: CorrectDischargeByPass       = .true.   !Default true Paulo suggestion
    end type T_Converge

    type     T_FromTimeSerieRunOff
        integer                                     :: ObjTimeSerie         = 0
        character(len=PathLength)                   :: FileName             = null_str
        integer                                     :: DataColumn           = null_int
    end type T_FromTimeSerieRunOff
    
    !level imposed as time serie
    type     T_ImposedLevelTS
        type(T_FromTimeSerieRunOff)                 :: TimeSerie
    end type T_ImposedLevelTS

    type T_BoundaryLine
        real                                        :: WaterLevelValue      = null_real
        character(PathLength)                       :: FileName             = null_str
        type (T_Lines),   pointer                   :: Line                 => null()
        logical                                     :: Variable             = .false.
        type(T_FromTimeSerieRunOff)                 :: TimeSerie
        integer                                     :: nCells               = null_int
        integer, dimension(:), allocatable          :: I
        integer, dimension(:), allocatable          :: J
    end type T_BoundaryLine

  
    type  T_RunOff
        integer                                     :: InstanceID               = 0
        character(len=StringLength)                 :: ModelName                = null_str
        integer                                     :: ObjBasinGeometry         = 0
        integer                                     :: ObjTime                  = 0
        integer                                     :: ObjHorizontalGrid        = 0
        integer                                     :: ObjHorizontalMap         = 0
        integer                                     :: ObjGridData              = 0
        integer                                     :: ObjHDF5                  = 0
        integer                                     :: ObjDrainageNetwork       = 0
        integer                                     :: ObjDischarges            = 0
        integer                                     :: ObjEnterData             = 0
        integer                                     :: ObjBoxDif                = 0
        integer                                     :: ObjTimeSerie             = 0
        type (T_OutPutRunOff)                       :: OutPut
        type (T_ExtVarRunOff)                       :: ExtVar
        type (T_FilesRunOff)                        :: Files
        type (T_Time)                               :: BeginTime
        type (T_Time)                               :: EndTime
        real(8), dimension(:,:), pointer            :: myWaterLevel             => null()
        real(8), dimension(:,:), pointer            :: myWaterColumn            => null()
        real,    dimension(:,:), pointer            :: InitialWaterColumn       => null()
        real,    dimension(:,:), pointer            :: InitialWaterLevel        => null()
        logical                                     :: PresentInitialWaterColumn = .false.
        logical                                     :: PresentInitialWaterLevel  = .false.
        real(8), dimension(:,:), pointer            :: myWaterVolume            => null() 
        real(8), dimension(:,:), pointer            :: myWaterColumnOld         => null() !OldColumn from Basin
        real(8), dimension(:,:), pointer            :: myWaterColumnAfterTransport => null() !for property transport
        real(8), dimension(:,:), pointer            :: myWaterVolumePred        => null() !to avoid negative collumns
        real(8), dimension(:,:), pointer            :: myWaterVolumeOld         => null()
        real,    dimension(:,:), pointer            :: lFlowToChannels          => null() !Instantaneous Flow To Channels
        real,    dimension(:,:), pointer            :: iFlowToChannels          => null() !Integrated    Flow
        real,    dimension(:,:), pointer            :: iFlowRouteDFour          => null() !Integrated Route D4 flux
        real,    dimension(:,:), pointer            :: lFlowBoundary            => null() !Instantaneous Flow to impose BC
        real,    dimension(:,:), pointer            :: iFlowBoundary            => null() !Integrated    Flow to impose BC
        real,    dimension(:,:), pointer            :: lFlowDischarge           => null() !Instantaneous Flow of discharges
        real,    dimension(:,:), pointer            :: iFlowDischarge           => null() !Integrated    Flow of discharges
        real(8), dimension(:,:), pointer            :: lFlowX, lFlowY           => null() !Instantaneous OverLandFlow (LocalDT   )
        real(8), dimension(:,:), pointer            :: iFlowX, iFlowY           => null() !Integrated    OverLandFlow (AfterSumDT)
        real(8), dimension(:,:), pointer            :: FlowXOld, FlowYOld       => null() !Flow From previous iteration
        real(8), dimension(:,:), pointer            :: InitialFlowX, InitialFlowY => null() !Initial Flow of convergence
        real(8), dimension(:,:), pointer            :: VelModFaceU, VelModFaceV => null() !Flow From previous iteration
        real,    dimension(:,:), pointer            :: AreaU, AreaV             => null()
        integer, dimension(:,:), pointer            :: ComputeFaceU             => null()
        integer, dimension(:,:), pointer            :: ComputeFaceV             => null()
        integer, dimension(:,:), pointer            :: ComputeAdvectionU        => null()
        integer, dimension(:,:), pointer            :: ComputeAdvectionV        => null()
        real,    dimension(:,:), pointer            :: NoAdvectionPoints        => null()

        real                                        :: GridCosAngleX, GridCosAngleY, GridSinAngleX, GridSinAngleY

        integer, dimension(:,:), pointer            :: OpenPoints               => null() !Mask for gridcells above min watercolumn
        real,    dimension(:,:), pointer            :: OverLandCoefficient      => null() !Manning or Chezy
        real,    dimension(:,:), pointer            :: OverLandCoefficientDelta => null() !For erosion/deposition
        real,    dimension(:,:), pointer            :: OverLandCoefficientX     => null() !Manning or Chezy
        real,    dimension(:,:), pointer            :: OverLandCoefficientY     => null() !Manning or Chezy
        real,    dimension(:,:), pointer            :: MassError                => null() !Contains mass error
        real,    dimension(:,:), pointer            :: CenterFlowX              => null()
        real,    dimension(:,:), pointer            :: CenterFlowY              => null()
        real,    dimension(:,:), pointer            :: CenterVelocityX          => null()
        real,    dimension(:,:), pointer            :: CenterVelocityY          => null()
        real,    dimension(:,:), pointer            :: FlowModulus              => null()
        real,    dimension(:,:), pointer            :: VelocityModulus          => null()
        integer, dimension(:,:), pointer            :: LowestNeighborI          => null() !Lowest Neighbor in the surroundings
        integer, dimension(:,:), pointer            :: LowestNeighborJ          => null() !Lowest Neighbor in the surroundings       
        integer, dimension(:,:), pointer            :: DFourSinkPoint           => null() !Point which can't drain with in X/Y only
        integer, dimension(:,:), pointer            :: StabilityPoints          => null() !Points where models check stability
        type(T_PropertyID)                          :: OverLandCoefficientID, NoAdvectionZonesID
        
        logical                                     :: StormWaterModel          = .false. !If connected to SWMM
        real                                        :: StormWaterModelDT        = -null_real

        type(T_SewerGEMSInlet),           dimension(:), allocatable :: Inlets
        type(T_SewerGEMSManhole),         dimension(:), allocatable :: Manholes
        type(T_SewerGEMSOutfall),         dimension(:), allocatable :: Outfalls
        type(T_SewerGEMSCrossSection),    dimension(:), allocatable :: CrossSections
        type(T_SewerGEMSPond),            dimension(:), allocatable :: Ponds
        type(T_SewerGEMSOpenChannelLink), dimension(:), allocatable :: OpenChannelLinks
        type(T_SewerGEMSHeadwall),        dimension(:), allocatable :: Headwalls


        integer                                     :: NumberOfInlets           = null_int
        integer                                     :: NumberOfHeadwalls        = null_int
        integer                                     :: NumberOfNodes            = null_int
        integer                                     :: NumberOfIgnoredNodes     = null_int
        character(len=99), dimension(:), allocatable:: IgnoredNodes

        real, dimension(:), allocatable             :: NodesX, NodesY
        integer, dimension(:), allocatable          :: NodesI, NodesJ, NodesID, NodesCellID
        integer, dimension(:), allocatable          :: NodesType
        character(len=99), dimension(:), allocatable:: NodesNames

        integer                                     :: ActiveNodes              = null_int
        integer                                     :: InactiveNodes            = null_int
        integer                                     :: NumberOfManholes         = null_int
        integer                                     :: NumberOfOutfalls         = null_int
        integer                                     :: NumberOfCrossSections    = null_int
        integer                                     :: NumberOfPonds            = null_int
        integer                                     :: NumberOfOpenChannelLinks = null_int

        real, dimension(:,:), pointer               :: NodeRiverLevel           => null() !river level at river points (from DN or external model)
        integer, dimension(:,:), pointer            :: NodeRiverMapping         => null() !mapping of river points where interaction occurs (for external model)
        integer, dimension(:,:), pointer            :: RiverNodeMap             => null() !i,j indexes of grid cell of river points where interaction occurs (for external model)
        real, dimension(:,:), pointer               :: MarginRiverLevel         => null() !river level at margin points
        real, dimension(:,:), pointer               :: MarginFlowToChannels     => null() !flow to channels at margin points
        
        real                                        :: MinSlope                 = null_real
        logical                                     :: AdjustSlope              = .false.
        logical                                     :: Stabilize                = .false.
        logical                                     :: Discharges               = .false.
        logical                                     :: MomentumDischarges       = .false.
        logical                                     :: RouteDFourPoints         = .false.
        logical                                     :: RouteDFourPointsOnDN     = .false.
        integer                                     :: RouteDFourMethod         = null_int

        logical                                     :: Buildings                = .false.
!        real                                        :: StabilizeFactor       = null_real
!        real                                        :: StabilizeHardCutLimit = 0.1
        integer                                     :: HydrodynamicApproximation= DiffusionWave_
        logical                                     :: CalculateAdvection       = .true.
        logical                                     :: NoAdvectionZones         = .false.
        logical                                     :: CalculateDiffusion       = .false.
        real                                        :: Viscosity_Coef           = null_real
        logical                                     :: CalculateCellMargins     = .true.
        logical                                     :: ImposeMaxVelocity        = .false.
        real                                        :: ImposedMaxVelocity       = 0.1
!        integer                                     :: LastGoodNiter           = 1
!        integer                                     :: NextNiter               = 1
!        real                                        :: InternalTimeStepSplit   = 1.5
!        integer                                     :: MinIterations           = 1
!        integer                                     :: MinToRestart            = 0
        real                                        :: MinimumWaterColumn       = null_real
        real                                        :: MinimumWaterColumnAdvection = null_real
!        real                                        :: MinimumWaterColumnStabilize = null_real
!        real                                        :: NextDT               = null_real
!        real                                        :: DTFactor             = null_real
!        real                                        :: DTFactorUp           = null_real
!        real                                        :: DTFactorDown         = null_real
!        real                                        :: CurrentDT            = null_real
!        logical                                     :: LimitDTCourant       = .false.
!        logical                                     :: LimitDTVariation     = .true.
!        real                                        :: MaxCourant           = 1.0        
        logical                                     :: ImposeBoundaryValue      = .false.
        logical                                     :: AllowBoundaryInflow      = .false.
        logical                                     :: BoundaryImposedLevelInTime = .false.
        real                                        :: BoundaryValue            = null_real
        type(T_ImposedLevelTS)                      :: ImposedLevelTS      
        real                                        :: MaxDtmForBoundary        = null_real
        integer                                     :: BoundaryMethod           = null_int
        logical                                     :: HasBoundaryLines         = .false.
        integer, dimension(:,:), pointer            :: BoundaryCells            => null()
        integer                                     :: NumberOfBoundaryLines    = null_int
        type(T_BoundaryLine), dimension(:), allocatable :: BoundaryLines
        real, dimension(:,:), allocatable           :: WaterLevelBoundaryValue

!        integer                                     :: MaxIterations        = 5
        logical                                     :: SimpleChannelInteraction = .false.
        logical                                     :: ChannelHasTwoGridPoints  = .false.
        logical                                     :: LimitToCriticalFlow      = .true.
        integer                                     :: FaceWaterColumn          = WCMaxBottom_
!        real                                        :: MaxVariation            = null_real
        integer                                     :: OverlandChannelInteractionMethod = null_int
!        logical                                     :: CheckDecreaseOnly       = .false.

        type(T_Converge)                            :: CV !Convergence options
        
        real(8)                                     :: BoundaryFlowVolume        = 0.0 !m3 => positive if flow is towards boundary.          
        real(8)                                     :: VolumeStoredInSurface     = 0.0
        real(8)                                     :: VolumeStoredInStormSystem = 0.0
        real(8)                                     :: TotalDischargeFlowVolume  = 0.0  
        real(8)                                     :: TotalBoundaryFlowVolume   = 0.0     
        real(8)                                     :: TotalBoundaryInflowVolume = 0.0     
        real(8)                                     :: TotalBoundaryOutFlowVolume= 0.0     
        real(8)                                     :: TotalInfiltrationVolume   = 0.0     
        real(8)                                     :: TotalRainfallVolume       = 0.0  
        real(8)                                     :: TotalStoredVolume         = 0.0
        real                                        :: InitialTotalVolume        = 0.0
        real                                        :: TotalStormWaterVolume     = 0.0
        real                                        :: TotalInletsVolume         = 0.0
        real                                        :: TotalManholesVolume       = 0.0
        real                                        :: TotalHeadwallsVolume      = 0.0
        real                                        :: TotalPondsVolume          = 0.0
        real                                        :: TotalOpenChannelVolume    = 0.0
        real                                        :: TotalOutfallsVolume       = 0.0
        real                                        :: Total1DVolume             = 0.0
        real                                        :: Total1D2DVolume           = 0.0
        real                                        :: MaxTotal1DVolume          = 0.0
        real                                        :: MaxTotal2DVolume          = 0.0
        real                                        :: MaxTotal1D2DVolume        = 0.0
        real                                        :: TimeOfMaxTotal1DVolume    = 0.0
        real                                        :: TimeOfMaxTotal2DVolume    = 0.0
        real                                        :: TimeOfMaxTotal1D2DVolume  = 0.0

        real(8)                                     :: AvrgAccInfiltrationDepth  = 0.0   
        real(8)                                     :: AvrgAccInfiltrationVolume = 0.0   
   
        
        logical                                     :: Continuous          = .false.
        logical                                     :: StopOnWrongDate     = .true.
        

        integer                                     :: BasinCellsCount    = 0

        !Grid size
        type (T_Size2D)                             :: Size
        type (T_Size2D)                             :: WorkSize
        
        
        type(T_NodeGridPoint    ), pointer          :: FirstNodeGridPoint        => null()
        type(T_NodeGridPoint    ), pointer          :: LastNodeGridPoint         => null()
        integer                                     :: NodeGridPointNumber     = 0        
        type(NodeGridPointPtr_class), dimension(:), allocatable :: NodeGridPointArray
        
        type(T_MarginGridPoint    ), pointer        :: FirstMarginGridPoint        => null()
        type(T_MarginGridPoint    ), pointer        :: LastMarginGridPoint         => null()
        integer                                     :: MarginGridPointNumber     = 0
        type(MarginGridPointPtr_class), dimension(:), allocatable :: MarginGridPointArray
        
        type(T_BankGridPoint    ), pointer          :: FirstBankGridPoint        => null()
        type(T_BankGridPoint    ), pointer          :: LastBankGridPoint         => null()
        integer                                     :: BankGridPointNumber     = 0
        type(BankGridPointPtr_class), dimension(:), allocatable :: BankGridPointArray
        
        logical                                     :: Use1D2DInteractionMapping = .false.
        integer                                     :: Counter = 0
        
        type(T_RunOff), pointer                     :: Next                 => null()
    end type  T_RunOff


    !Global Module Variables
    type (T_RunOff), pointer                        :: FirstObjRunOff       => null()
    type (T_RunOff), pointer                        :: Me                   => null()

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructRunOff(ModelName,                                       &            
                               RunOffID,                                        &
                               ComputeTimeID,                                   &
                               HorizontalGridID,                                &
                               HorizontalMapID,                                 &
                               GridDataID,                                      &
                               BasinGeometryID,                                 &
                               DrainageNetworkID,                               &
                               DischargesID,                                    &
                               STAT)

        !Arguments---------------------------------------------------------------
        character(len=*)                                :: ModelName
        integer                                         :: RunOffID
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: GridDataID
        integer                                         :: BasinGeometryID
        integer                                         :: DrainageNetworkID
        integer, optional, intent(OUT)                  :: STAT     
        integer, intent (OUT)                           :: DischargesID

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL

        !------------------------------------------------------------------------
        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mRunOff_)) then
            nullify (FirstObjRunOff)
            call RegisterModule (mRunOff_) 
        endif

        call Ready(RunOffID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%ModelName = ModelName
            
            !Associates External Instances
            Me%ObjTime            = AssociateInstance (mTIME_           , ComputeTimeID     )
            Me%ObjHorizontalGrid  = AssociateInstance (mHORIZONTALGRID_ , HorizontalGridID  )
            Me%ObjHorizontalMap   = AssociateInstance (mHORIZONTALMAP_  , HorizontalMapID   )
            Me%ObjGridData        = AssociateInstance (mGRIDDATA_       , GridDataID        )
            Me%ObjBasinGeometry   = AssociateInstance (mBASINGEOMETRY_  , BasinGeometryID   )

            if (DrainageNetworkID /= 0) then
                Me%ObjDrainageNetwork   = AssociateInstance (mDRAINAGENETWORK_, DrainageNetworkID)
            endif
            
            !Time Stuff
            call GetComputeTimeLimits   (Me%ObjTime, BeginTime = Me%BeginTime,           &
                                         EndTime = Me%EndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunOff - ModuleRunOff - ERR010'

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunOff - ModuleRunOff - ERR011'
        
            Me%CV%NextNiteration = 1
            Me%CV%CurrentDT = Me%ExtVar%DT

            call ReadLockExternalVar (StaticOnly = .false.)
            
            call CheckHorizontalGridRotation

            !Gets the size of the grid
            call GetHorizontalGridSize (Me%ObjHorizontalGrid,                            &
                                        Size     = Me%Size,                              &
                                        WorkSize = Me%WorkSize,                          &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunOff - ModuleRunOff - ERR020'
            
            call AllocateVariables
            
            call ReadDataFile

            call InitializeVariables

            call ConstructOverLandCoefficient

            call ModifyGeometryAndMapping

            !Checks if River Network is consistent with the one previously constructed
            if (DrainageNetworkID /= 0) then
                call CheckRiverNetWorkConsistency
            endif

            !Constructs Discharges
            if (Me%Discharges) then
                call ConstructDischarges
            endif
           
            !Constructs SewerStormWaterNodesMap
            if (Me%StormWaterModel) then
                call ConstructSewerStormWater     
            endif
            
            !Constructs Boundary Cells
            if (Me%ImposeBoundaryValue) then

                call ConstructWaterLevelBoundaryConditions

                if (Me%BoundaryImposedLevelInTime)then
                    call ModifyBoundaryLevel
                endif
                
            endif
            
            !Reads conditions from previous run
            if (Me%Continuous) then
                if (Me%OutPut%RestartFormat == BIN_) then
                    call ReadInitialFile_Bin
                else if (Me%OutPut%RestartFormat == HDF_) then
                    call ReadInitialFile_Hdf
                endif
            endif
            
            if (Me%OutPut%Yes) then
                call ConstructHDF5Output
            endif

            call CalculateTotalStoredVolume

            Me%InitialTotalVolume = Me%TotalStoredVolume

            !Output Results
            if (Me%OutPut%Yes .or. Me%OutPut%TimeSeries) then
                call ComputeCenterValues      
            endif
            
            if(Me%OutPut%Yes)then
                call RunOffOutput
            endif
            
            if(Me%OutPut%TimeSeries) then
                call OutputTimeSeries
            endif

            if (Me%Output%WriteMaxWaterColumn .or. Me%Output%WriteMaxFloodRisk) then
                call OutputFlooding
            endif

            if (Me%Output%WriteFloodArrivalTime) then            
                call OutputFloodArrivalTime            
            endif

            call ComputeNextDT (Me%CV%NextNiteration)                                    

            call ReadUnLockExternalVar (StaticOnly = .false.)

            !Returns ID
            RunOffID          = Me%InstanceID
            DischargesID      = Me%ObjDischarges

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleRunOff - ConstructRunOff - ERR030' 

        end if cd0

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructRunOff
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_RunOff), pointer                         :: NewObjRunOff
        type (T_RunOff), pointer                         :: PreviousObjRunOff


        !Allocates new instance
        allocate (NewObjRunOff)
        nullify  (NewObjRunOff%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjRunOff)) then
            FirstObjRunOff          => NewObjRunOff
            Me                      => NewObjRunOff
        else
            PreviousObjRunOff       => FirstObjRunOff
            Me                      => FirstObjRunOff%Next
            do while (associated(Me))
                PreviousObjRunOff   => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjRunOff
            PreviousObjRunOff%Next  => NewObjRunOff
        endif

        Me%InstanceID = RegisterNewInstance (mRUNOFF_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadDataFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        type(T_PropertyID)                          :: InitialWaterColumnID, InitialWaterLevelID
        type(T_PropertyID)                          :: OverLandCoefficientDeltaID
        character(len=PathLength)                   :: MappingFileName, RootSRT
        character(len=PathLength)                   :: IgnoredNodesFileName, OpenChannelLinksFileName
        character(len=PathLength)                   :: InletsFileName, PondsFileName, HeadwallsFileName
        integer                                     :: iflag, ClientNumber
        logical                                     :: BlockFound
        integer                                     :: i, j, n
        logical                                     :: DynamicAdjustManning
        real                                        :: dummy
        character(len = 1)                          :: Char_n
        !Reads the name of the data file from nomfich
        call ReadFileName ('RUNOFF_DATA', Me%Files%DataFile, "RunOff Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR010'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('RUNOFF_HDF', Me%Files%TransientHDF, "RunOff HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR020'

        call ReadFileName('RUNOFF_FIN', Me%Files%FinalFile,                             &
                           Message = "RunOff Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR030'
        
        call ReadFileName('ROOT_SRT', RootSRT,                                          &
                           Message = "Mass Error File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR031'
        Me%Files%MassErrorFile = trim(adjustl(RootSRT))//"MassError.dat"
        Me%Files%RunOffLogFile = trim(adjustl(RootSRT))//"RunOff.log"

        !Constructs the DataFile
        call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR040'

        !Initial Water Column
        call GetData(dummy,                                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'INITIAL_WATER_COLUMN',                              &
                     default      = 0.0,                                                 & 
                     ClientModule = 'ModuleRunOff',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR050'
        
        if (iflag /= 0) then
            write(*,*)'The keyword INITIAL_WATER_COLUMN is obsolete.'
            write(*,*)'Please use the block <BeginInitialWaterColumn> / <EndInitialWaterColumn>'
            stop 'ReadDataFile - ModuleRunOff - ERR060'
        endif        

        !Gets Block 
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                    &
                                    '<BeginInitialWaterColumn>',                      &
                                    '<EndInitialWaterColumn>', BlockFound,            &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR070'

        if (BlockFound) then
            call ConstructFillMatrix  ( PropertyID       = InitialWaterColumnID,         &
                                        EnterDataID      = Me%ObjEnterData,              &
                                        TimeID           = Me%ObjTime,                   &
                                        HorizontalGridID = Me%ObjHorizontalGrid,         &
                                        ExtractType      = FromBlock,                    &
                                        PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                        Matrix2D         = Me%InitialWaterColumn,        &
                                        TypeZUV          = TypeZ_,                       &
                                        STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR080'

            call KillFillMatrix(InitialWaterColumnID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR081'
            
            Me%PresentInitialWaterColumn = .true.
            
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR090'

            

        else

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                    &
                                        '<BeginInitialWaterLevel>',                       &
                                        '<EndInitialWaterLevel>', BlockFound,             &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR091'

            if (BlockFound) then
                call ConstructFillMatrix  ( PropertyID       = InitialWaterLevelID,          &
                                            EnterDataID      = Me%ObjEnterData,              &
                                            TimeID           = Me%ObjTime,                   &
                                            HorizontalGridID = Me%ObjHorizontalGrid,         &
                                            ExtractType      = FromBlock,                    &
                                            PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                            Matrix2D         = Me%InitialWaterLevel,         &
                                            TypeZUV          = TypeZ_,                       &
                                            STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR092'

                call KillFillMatrix(InitialWaterLevelID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR093'  
                
                Me%PresentInitialWaterLevel = .true.
                
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR094'
                
            else
                write(*,*)
                write(*,*)'Missing Block <BeginInitialWaterColumn> / <EndInitialWaterColumn>' 
                write(*,*)'or <BeginInitialWaterLevel> / <EndInitialWaterLevel>' 
                stop      'ReadDataFile - ModuleRunOff - ERR100'
            endif
        endif

         !Gets Minimum Slope 
        call GetData(Me%MinSlope,                                               &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'MIN_SLOPE',                                &
                     default      = 0.0,                                        &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0110'

        if (Me%MinSlope < 0.0 .or. Me%MinSlope >= 1.) then
            write (*,*) 'Invalid Minimum Slope [MIN_SLOPE]'
            stop 'ReadDataFile - ModuleRunOff - ERR0120'
        end if

        !Adjusts Slope according to
        !http://www.hkh-friend.net.np/rhdc/training/lectures/HEGGEN/Tc_3.pdf
        call GetData(Me%AdjustSlope,                                            &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'ADJUST_SLOPE',                             &
                     default      = .true.,                                     &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0130'


        !Gets Routing method
        call GetData(Me%HydrodynamicApproximation,                              &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'HYDRODYNAMIC_APROX',                       &
                     default      = DiffusionWave_,                             &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0140'

        if (Me%HydrodynamicApproximation /= KinematicWave_ .and.                &
            Me%HydrodynamicApproximation /= DiffusionWave_ .and.                &
            Me%HydrodynamicApproximation /= DynamicWave_) then
            write (*,*) 'Invalid Hydrodynamic Approximation [HYDRODYNAMIC_APROX]'
            stop 'ReadDataFile - ModuleRunOff - ERR0150'
        end if     
   
        if (Me%HydrodynamicApproximation == DynamicWave_) then

            !Gets if advection is to be calculated
            call GetData(Me%CalculateAdvection,                                 &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromFile,                               &
                         keyword      = 'ADVECTION',                            &
                         default      = .true.,                                 &
                         ClientModule = 'ModuleRunOff',                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0160'

            
            if (Me%CalculateAdvection) then    
                
                !Minimum Water Column for advection computation
                call GetData(Me%MinimumWaterColumnAdvection,                                     &
                             Me%ObjEnterData, iflag,                                             &
                             SearchType   = FromFile,                                            &
                             keyword      = 'MIN_WATER_COLUMN_ADVECTION',                        &
                             default      = 0.0,                                                 &
                             ClientModule = 'ModuleRunOff',                                      &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0170'
                

                call GetData(Me%NoAdvectionZones,                                                &
                             Me%ObjEnterData, iflag,                                             &
                             SearchType   = FromFile,                                            &
                             keyword      = 'NO_ADVECTION_ZONES',                                &
                             default      = .false.,                                             &
                             ClientModule = 'ModuleRunOff',                                      &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0170a'
                
                if(Me%NoAdvectionZones)then
                    
                    call ConstructAdvectionZones
                    
                else
                    
                    call SetMatrixValue(Me%ComputeAdvectionU, Me%Size, 1)
                    call SetMatrixValue(Me%ComputeAdvectionU, Me%Size, 1)
                    
                endif
                
            else

                call SetMatrixValue(Me%ComputeAdvectionU, Me%Size, 0)
                call SetMatrixValue(Me%ComputeAdvectionU, Me%Size, 0)
             
            endif
            
        endif
        
        

        !Method for computing water column in the face (1 - Using max level and max bottom; 
        !2- using max level and average of bottom)
        call GetData(Me%FaceWaterColumn,                                    &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'WATER_COLUMN_FACE',                    &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = WCMaxBottom_,                           &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR175'
        
        if (Me%FaceWaterColumn /= WCMaxBottom_ .and. Me%FaceWaterColumn /= WCAverageBottom_) then
            write(*,*) 'Unknown option for WATER_COLUMN_FACE'
            stop 'ReadDataFile - ModuleRunOff - ERR176'
        endif        
        
        if (Me%FaceWaterColumn == WCMaxBottom_) then
            !Gets if compute "margins" aside of adjacent cells that produce friction
            call GetData(Me%CalculateCellMargins,                               &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromFile,                               &
                         keyword      = 'HYDRAULIC_RADIUS_MARGINS',             &
                         default      = .true.,                                 &
                         ClientModule = 'ModuleRunOff',                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR180'        
        endif
        
        !Gets if solution is limited by an maximum velocity
        call GetData(Me%ImposeMaxVelocity,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'IMPOSE_MAX_VELOCITY',                      &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR190'

        if (Me%ImposeMaxVelocity) then
        
            !Gets if solution is limited by an maximum velocity
            call GetData(Me%ImposedMaxVelocity,                                     &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'MAX_VELOCITY',                             &
                         default      = 0.1,                                        &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR200'
        
        endif


        !Gets if Manning Coeficient is increased with water depth
        call GetData(DynamicAdjustManning,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'DYNAMIC_ADJUST_MANNING',                   &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR210'

        if (iflag > 0 .and. .not. DynamicAdjustManning) then
            write(*,*)'The option DynamicAdjustManning (DYNAMIC_ADJUST_MANNING) has been removed.'
            write(*,*)'Please review your runoff data file!'
        endif

        if (DynamicAdjustManning) then
            write(*,*)'The option DynamicAdjustManning (DYNAMIC_ADJUST_MANNING) has been removed.'
            write(*,*)'Please review your runoff data file!'
            stop
        endif

        !Minimum Water Column for overland flow
        call GetData(Me%MinimumWaterColumn,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MIN_WATER_COLUMN',                                  &
                     ClientModule = 'ModuleRunOff',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR220'
        if (iflag == 0) then
            write(*,*)'MIN_WATER_COLUMN must be defined in module Runoff'
            stop 'ReadDataFile - ModuleRunOff - ERR230'
        endif

        !Continuous Computation
        call GetData(Me%Continuous,                                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CONTINUOUS',                                        &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleRunoff',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR240'

        if (Me%Continuous) then
            call ReadFileName('RUNOFF_INI', Me%Files%InitialFile,                         &
                               Message = "Runoff Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0250'
        endif

        call GetData(Me%StopOnWrongDate,                                    &
                     Me%ObjEnterData, iflag,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'STOP_ON_WRONG_DATE',                   &
                     default      = .true.,                                 &
                     ClientModule = 'ModuleRunOff',                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR260'

        !Impose Boundary Value
        call GetData(Me%ImposeBoundaryValue,                                    &
                     Me%ObjEnterData, iflag,                                    &  
                     keyword      = 'IMPOSE_BOUNDARY_VALUE',                    &
                     ClientModule = 'ModuleRunOff',                             &
                     SearchType   = FromFile,                                   &
                     Default      = .false.,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR350'        
        
        if (Me%ImposeBoundaryValue) then
            
            !verify if the user wants to allow water to go in the domain (true) if
            !boundary level higher than water level or not (false) and the level imposed
            !behaves like a wall, only exits if higher and does not allow to get inside
            call GetData(Me%AllowBoundaryInflow,                                    &
                         Me%ObjEnterData, iflag,                                    &  
                         keyword      = 'ALLOW_BOUNDARY_INFLOW',                    &
                         ClientModule = 'ModuleRunOff',                             &
                         SearchType   = FromFile,                                   &
                         Default      = .false.,                                    &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR363'  

            call ReadBoundaryConditions
            
        endif
        
        
        !Discharges
        call GetData(Me%Discharges,                                         &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'DISCHARGES',                           &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR370'     
        
        if (Me%Discharges) then
            
            !Momentum discharges
            call GetData(Me%MomentumDischarges,                             &
                         Me%ObjEnterData, iflag,                            &
                         keyword    = 'MOMENTUM_DISCHARGE',                 &
                         Default    = .false.,                              &
                         SearchType = FromFile,                             &
                         ClientModule ='ModuleRunOff',                      &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR371'
            
            if (Me%MomentumDischarges) then
                write (*,*) "MOMENTUM_DISCHARGE is only allowed when HYDRODYNAMIC_APROX is 3 (Dynamic Wave)"
                stop 'ReadDataFile - ModuleRunOff - ERR372'
            end if
        
            !Discharges output time series
            call GetData(Me%Output%TimeSerieDischON,                        &
                         Me%ObjEnterData, iflag,                            &
                         keyword    = 'TIME_SERIE_DISCHARGES',              &
                         Default    = .false.,                              &
                         SearchType = FromFile,                             &
                         ClientModule ='ModuleRunOff',                      &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR375'

            if (Me%Output%TimeSerieDischON) then
            
                call GetData(Me%Output%DiscTimeSerieLocationFile,               &
                             Me%ObjEnterData,iflag,                             &
                             SearchType   = FromFile,                           &
                             keyword      = 'DISCHARGE_TIME_SERIE_LOCATION',    &
                             ClientModule = 'ModuleRunOff',                     &
                             Default      = Me%Files%DataFile,                  &
                             STAT         = STAT_CALL)
                if (STAT_CALL/=SUCCESS_)stop 'ReadDataFile - ModuleRunOff - ERR377'

                Me%Output%NextOutPutDisch = Me%BeginTime

                call GetData(Me%Output%OutPutDischDT,                           &
                             Me%ObjEnterData, iflag,                            &
                             keyword    = 'DISCHARGE_DT_OUTPUT_TIME',           &
                             Default    = FillValueReal,                        &
                             SearchType = FromFile,                             &
                             ClientModule ='ModuleRunOff',                      &
                             STAT       = STAT_CALL)            
                if (STAT_CALL/=SUCCESS_)stop 'ReadDataFile - ModuleRunOff - ERR378'

                if (iflag == 0) then
                    call GetComputeTimeStep(Me%ObjTime,                         &
                                            Me%Output%OutPutDischDT,            &
                                            STAT = STAT_CALL)
                    if (STAT_CALL/=SUCCESS_)stop 'ReadDataFile - ModuleRunOff - ERR379'

                endif                

            endif        
        endif    

        !Discharges
        call GetData(Me%SimpleChannelInteraction,                           &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'SIMPLE_CHANNEL_FLOW',                  &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR380'        
        
        !Check if DN link is made with more than one margin
        if (Me%SimpleChannelInteraction  .and. Me%ObjDrainageNetwork /= 0) then
            call GetHasTwoGridPoints(Me%ObjDrainageNetwork, Me%ChannelHasTwoGridPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0385'
        endif
        
        !Routes D4 Points
        call GetData(Me%RouteDFourPoints,                                   &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'ROUTE_D4',                             &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR390'
        
        if (Me%RouteDFourPoints) then
            !Routes D4 Points
            call GetData(Me%RouteDFourPointsOnDN,                               &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'ROUTE_D4_ON_DN',                       &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = .false.,                                &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR391'

            call GetData(Me%RouteDFourMethod,                                   &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'ROUTE_D4_METHOD',                      &
 !                        Default      = Celerity_,                              &
                         Default      = Manning_,                               &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR392'        

            if (Me%RouteDFourMethod /= Celerity_ .and. Me%RouteDFourMethod /= Manning_) then
                write(*,*)'ROUTE_D4_METHOD must be or 1 - Celerity based or 2 - Manning Equation'
                stop 'ReadDataFile - ModuleRunOff - ERR0393'
            endif
            
        endif

        !Limits Flow to critical
        call GetData(Me%LimitToCriticalFlow,                                &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'LIMIT_TO_CRITICAL_FLOW',               &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .true.,                                 &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR394'

        !If Buildings are to be simulated (flow ocuation in urban areas)
        call GetData(Me%Buildings,                                          &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'BUILDINGS',                            &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR430'

        if(Me%Buildings)then
            write(*,*) 
            write(*,*)"BUILDINGS in no longer an available option" 
            write(*,*)"Please review your model setup"
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR431'
        endif
        
        !If Connected to a StormWater model
        call GetData(Me%StormWaterModel,                                    &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'STORM_WATER_MODEL_LINK',               &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR440'

        !Gets Output Time 
        call GetOutPutTime(Me%ObjEnterData,                                              &
                           CurrentTime = Me%ExtVar%Now,                                  &
                           EndTime     = Me%EndTime,                                     &
                           keyword     = 'OUTPUT_TIME',                                  &
                           SearchType  = FromFile,                                       &
                           OutPutsTime = Me%OutPut%OutTime,                              &
                           OutPutsOn   = Me%OutPut%Yes,                                  &
                           STAT        = STAT_CALL)
        Me%OutPut%NextOutPut = 1    


        !Output for restart
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%ExtVar%Now,                                &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoff - ERR450'

        call GetData(Me%OutPut%RestartFormat,                                           &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_FORMAT',                              &
                     Default      = BIN_,                                               &
                     ClientModule = 'ModuleRunoff',                                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoff - ERR452'        
        if (Me%OutPut%RestartFormat /= BIN_ .and. Me%OutPut%RestartFormat /= HDF_) then
            write (*,*)
            write (*,*) 'RESTART_FILE_FORMAT options are: 1 - Binary or 2 - HDF'
            stop 'ReadDataFile - ModuleRunoff - ERR455'            
        endif        
        
        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleRunOff',                                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleRunoff - ERR460'



        call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0470'

        !Gets Block for OverLand Coef
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                      &
                                    '<BeginOverLandCoefficient>',                       &
                                    '<EndOverLandCoefficient>', BlockFound,             &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0480'
        if (BlockFound) then
            call ConstructFillMatrix  ( PropertyID       = Me%OverLandCoefficientID,     &
                                        EnterDataID      = Me%ObjEnterData,              &
                                        TimeID           = Me%ObjTime,                   &
                                        HorizontalGridID = Me%ObjHorizontalGrid,         &
                                        ExtractType      = FromBlock,                    &
                                        PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                        Matrix2D         = Me%OverLandCoefficient,       &
                                        TypeZUV          = TypeZ_,                       &
                                        STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0490'

            call KillFillMatrix(Me%OverLandCoefficientID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0500'
            
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0501'



            !Check that manning values entered are not zero or negative
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    if (.not. Me%OverLandCoefficient(i,j) .gt. 0.0) then
                        write(*,*) 'Found Manning Overland coefficient zero or negative in input'
                        write(*,*) 'in cell', i, j
                        stop 'ReadDataFile - ModuleRunoff - ERR0510'
                    endif
                
                
                endif
                
            enddo
            enddo

        else
            write(*,*)'Missing Block <BeginOverLandCoefficient> / <EndOverLandCoefficient>' 
            stop      'ReadDataFile - ModuleRunOff - ERR0520'
        endif
        

        
        !Gets Block for OverLand Coef Difference 
        !To compute overland resistance in bottom for shear computation (erosion/deposition).
        !This process was created to remove from manning the resistance given by 
        !aerial vegetation parts that affect flow but do not affect bottom shear. Without that,
        !a manning increase (e.g. forestation) in one cell increases water depth (and reduces velocity)
        !but may increase shear stress (because water height increase is transformed in bottom resistance 
        !using manning - chezy see module runoff properties)
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                    '<BeginOverLandCoefficientDelta>',           &
                                    '<EndOverLandCoefficientDelta>', BlockFound, &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0530'
        if (BlockFound) then
            call ConstructFillMatrix  ( PropertyID       = OverLandCoefficientDeltaID,   &
                                        EnterDataID      = Me%ObjEnterData,                 &
                                        TimeID           = Me%ObjTime,                   &
                                        HorizontalGridID = Me%ObjHorizontalGrid,         &
                                        ExtractType      = FromBlock,                    &
                                        PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                        Matrix2D         = Me%OverLandCoefficientDelta,  &
                                        TypeZUV          = TypeZ_,                       &
                                        STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0540'

            call KillFillMatrix(OverLandCoefficientDeltaID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0550'
            
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0551'


            !Check that final manning values are not zero or negative
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    if (.not. (Me%OverLandCoefficient(i,j) - Me%OverLandCoefficientDelta(i,j)) .gt. 0.0) then
                        write(*,*) 'Manning Overland coefficient delta found zero or negative in input'
                        write(*,*) 'in cell', i, j
                        stop 'ReadDataFile - ModuleRunoff - ERR0560'
                    endif
                
                
                endif
                
            enddo
            enddo

        else
            !Do not remove aerial vegetation effect from manning 
            Me%OverLandCoefficientDelta(:,:) = 0.0
        endif
        
        if (Me%StormWaterModel) then
        
            !Catch-basins/inlets
            call GetData(InletsFileName,                                                    &
                         Me%ObjEnterData,                                                   &
                         iflag,                                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'INLETS_FILENAME',                                  &
                         ClientModule = 'ModuleRunOff',                                     &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleRunoff - ERR696'
                
            if (iflag == 0) then
                write(*,*)
                write(*,*)"INLETS_FILENAME keyword setting the inlets file path was not found"
                write(*,*)"Simulation will not consider inlets/catch-basins"

                Me%NumberOfInlets = 0
                
            else

                call ReadInletsFromFile(InletsFileName, RootSRT)

            endif

            !Headwalls
            call GetData(HeadwallsFileName,                                                 &
                         Me%ObjEnterData,                                                   &
                         iflag,                                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'HEADWALLS_FILENAME',                               &
                         ClientModule = 'ModuleRunOff',                                     &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleRunoff - ERR697'
                
            if (iflag == 0) then
                write(*,*)
                write(*,*)"HEADWALLS_FILENAME keyword setting the inlets file path was not found"
                write(*,*)"Simulation will not consider headwalls"

                Me%NumberOfHeadwalls = 0
                
            else

                call ReadHeadwallsFromFile(HeadwallsFileName, RootSRT)

            endif

            !Ponds
            call GetData(PondsFileName,                                                     &
                         Me%ObjEnterData,                                                   &
                         iflag,                                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'PONDS_FILENAME',                                   &
                         ClientModule = 'ModuleRunOff',                                     &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleRunoff - ERR698'
                
            if (iflag == 0) then
                write(*,*)
                write(*,*)"PONDS_FILENAME keyword setting the ponds file path was not found"
                write(*,*)"Simulation will not consider ponds"

                Me%NumberOfPonds = 0
                
            else    
                call ReadPondsFromFile(PondsFileName)
            endif

            !Open channels links
            call GetData(OpenChannelLinksFileName,                                          &
                         Me%ObjEnterData,                                                   &
                         iflag,                                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'OPEN_CHANNEL_LINKS_FILENAME',                      &
                         ClientModule = 'ModuleRunOff',                                     &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleRunoff - ERR699'
                
            if (iflag == 0) then
                write(*,*)
                write(*,*)"OPEN_CHANNEL_LINKS_FILENAME keyword setting the inlets file path was not found"
                write(*,*)"Simulation will not consider inlets/catch-basins"

                Me%NumberOfOpenChannelLinks = 0
                
            else

                call ReadOpenChannelLinksFromFile(OpenChannelLinksFileName)

            endif


            !Ignored SWMM nodes (nodes inside grid that are not linked to 2D surface)
            call GetData(IgnoredNodesFileName,                                              &
                         Me%ObjEnterData,                                                   &
                         iflag,                                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'IGNORED_NODES_FILENAME',                           &
                         ClientModule = 'ModuleRunOff',                                     &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleRunoff - ERR700'
                
            if (iflag == 0) then
                write(*,*)
                write(*,*)"IGNORED_NODES_FILENAME keyword setting the ponds file path was not found"
                write(*,*)"Simulation will link all SWMM nodes inside the 2D grid"

                Me%NumberOfIgnoredNodes = 0
                
            else    
                call ReadIgnoredNodesFromFile(IgnoredNodesFileName)
            endif
                
        endif

                
        !Get mapping to river in case DN 1D river/2D floodplain model 
        if (Me%ObjDrainageNetwork /= 0) then
            !Get file with 1D interactions (for External 1D Model interpolation of level and integration of computed flow)
            call GetData(MappingFileName,                                          &
                            Me%ObjEnterData,iflag,                                 &
                            SearchType   = FromFile,                               &
                            keyword      = '1D_INTERACTION_MAPPING_FILE',          &
                            ClientModule = 'ModuleRunoff',                         &
                            STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR800' 
                               
            
            if (iflag .EQ. 1) then
                Me%Use1D2DInteractionMapping = .true.
                call Read1DInteractionMapping(MappingFileName)
            endif
            
        endif

        
        
        

        !Write Max Flow Modulus File 
        call GetData(Me%Output%WriteMaxFlowModulus,                             &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'WRITE_MAX_FLOW_FILE',                      &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR720'

        if(Me%Output%WriteMaxFlowModulus) then
            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%Output%MaxFlowModulusFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR730'
            Me%Output%MaxFlowModulusFile = trim(adjustl(Me%Output%MaxFlowModulusFile))//"MaxRunOff.dat"
        end if
              

        !Write all 3 flood layers
        call GetData(Me%Output%OutputFloodRisk,                                 &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'OUTPUT_FLOOD_RISK',                        &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR731'          
        
        if (Me%Output%OutputFloodRisk) then
            Me%Output%WriteMaxWaterColumn           = .true.
            Me%Output%WriteVelocityAtMaxWaterColumn = .true.
            Me%Output%WriteMaxFloodRisk             = .true.  
            Me%Output%WriteFloodPeriod              = .true.  
            Me%Output%WriteFloodArrivalTime         = .true.            

            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%Output%MaxWaterColumnFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0732'
            Me%Output%MaxWaterLevelFile         = trim(adjustl(Me%Output%MaxWaterColumnFile))//"MaxWaterLevel.dat"
            Me%Output%TimeOfMaxWaterColumnFile  = trim(adjustl(Me%Output%MaxWaterColumnFile))//"TimeMaxWaterColumn.dat"
            Me%Output%MaxWaterColumnFile = trim(adjustl(Me%Output%MaxWaterColumnFile))//"MaxWaterColumn.dat"
            
            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%Output%VelocityAtMaxWaterColumnFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0733'
            Me%Output%VelocityAtMaxWaterColumnFile = trim(adjustl(Me%Output%VelocityAtMaxWaterColumnFile))&
                                                        &//"VelocityAtMaxWaterColumn.dat"            

            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%Output%MaxFloodRiskFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0734'
            Me%Output%MaxFloodRiskFile = trim(adjustl(Me%Output%MaxFloodRiskFile))//"MaxFloodRisk.dat"                             

            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%Output%FloodPeriodFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0740'

            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%Output%FloodArrivalTimeFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0740'
            Me%Output%FloodArrivalTimeFile = trim(adjustl(Me%Output%FloodArrivalTimeFile))//"FloodArrivalTime.dat"     

            
        else        
        
            !Write Max water column 
            call GetData(Me%Output%WriteMaxWaterColumn,                             &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'WRITE_MAX_WATER_COLUMN',                   &
                         default      = .true.,                                     &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR740'

            if(Me%Output%WriteMaxWaterColumn) then
                !Gets the root path from the file nomfich.dat
                call ReadFileName("ROOT_SRT", Me%Output%MaxWaterColumnFile, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0750'
                Me%Output%MaxWaterLevelFile         = trim(adjustl(Me%Output%MaxWaterColumnFile))//"MaxWaterLevel.dat"
                Me%Output%TimeOfMaxWaterColumnFile  = trim(adjustl(Me%Output%MaxWaterColumnFile))//"TimeMaxWaterColumn.dat"
                Me%Output%MaxWaterColumnFile        = trim(adjustl(Me%Output%MaxWaterColumnFile))//"MaxWaterColumn.dat"
            
                !Write velocity at maximum water column 
                call GetData(Me%Output%WriteVelocityAtMaxWaterColumn,                   &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromFile,                                   &
                             keyword      = 'WRITE_VELOCITY_AT_MAX_WATER_COLUMN',       &
                             default      = .false.,                                    &
                             ClientModule = 'ModuleRunOff',                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR760'            

                if(Me%Output%WriteVelocityAtMaxWaterColumn) then
                    !Gets the root path from the file nomfich.dat
                    call ReadFileName("ROOT_SRT", Me%Output%VelocityAtMaxWaterColumnFile, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0770'
                    Me%Output%VelocityAtMaxWaterColumnFile = trim(adjustl(Me%Output%VelocityAtMaxWaterColumnFile))&
                                                             &//"VelocityAtMaxWaterColumn.dat"
            
                end if
            endif

            !Write max water column * velocity
            call GetData(Me%Output%WriteMaxFloodRisk,                               &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'WRITE_MAX_FLOOD_RISK',                     &
                         default      = .false.,                                    &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR780'          

        
            if(Me%Output%WriteMaxFloodRisk) then
                !Gets the root path from the file nomfich.dat
                call ReadFileName("ROOT_SRT", Me%Output%MaxFloodRiskFile, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0790'
                Me%Output%MaxFloodRiskFile = trim(adjustl(Me%Output%MaxFloodRiskFile))//"MaxFloodRisk.dat"   
            endif
            
            !Write flood period
            call GetData(Me%Output%WriteFloodPeriod,                                &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'WRITE_FLOOD_PERIOD',                       &
                         default      = .false.,                                    &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR792'          

        
            if(Me%Output%WriteFloodPeriod) then
                !Gets the root path from the file nomfich.dat
                call ReadFileName("ROOT_SRT", Me%Output%FloodPeriodFile, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0793'
            endif
      
        endif
        
        !factor for velocity in flood risk
        if (Me%Output%WriteMaxFloodRisk) then
            call GetData(Me%Output%FloodRiskVelCoef,                                &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'FLOOD_RISK_VEL_COEF',                      &
                         default      = 0.5,                                        &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR794'            
        endif

        !water column limit above which the cell is considered flooded
        if (Me%Output%WriteFloodPeriod) then

            call GetData(Me%Output%nFloodPeriodLimits,                              &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'FLOOD_PERIOD_N_CLASSES',                   &
                         default      = 0,                                          &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR795'

            if(Me%Output%nFloodPeriodLimits > 0)then

                if(Me%Output%nFloodPeriodLimits > 9)then
                    write(*,*)"Maximum allowed number of FLOOD_PERIOD_N_CLASSES is 9." 
                    stop 'ReadDataFile - ModuleRunOff - ERR795a'
                endif

                allocate(Me%Output%FloodPeriodWaterColumnLimits(1:Me%Output%nFloodPeriodLimits))

                call GetData(Me%Output%FloodPeriodWaterColumnLimits,                &
                             Me%ObjEnterData, iflag,                                &
                             SearchType   = FromFile,                               &
                             keyword      = 'FLOOD_PERIOD_CLASSES',                 &
                             ClientModule = 'ModuleRunOff',                         &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    if (STAT_CALL /= SIZE_ERR_) then
                        write(*,*) "Wrong number of depth limits FLOOD_PERIOD_CLASSES"
                    endif
                    stop 'ReadDataFile - ModuleRunOff - ERR796'
                endif

                if(iflag == 0)then
                    write(*,*) "Could not find keyword FLOOD_PERIOD_CLASSES"
                    stop 'ReadDataFile - ModuleRunOff - ERR797'
                endif

                allocate (Me%Output%FloodPeriods(Me%Size%ILB:Me%Size%IUB, &
                                                 Me%Size%JLB:Me%Size%JUB, &
                                                 1:Me%Output%nFloodPeriodLimits)) 
                Me%Output%FloodPeriods = 0.
                
                allocate (Me%Output%FloodPeriodFiles(1:Me%Output%nFloodPeriodLimits)) 

                do n = 1, Me%Output%nFloodPeriodLimits

                    write(Char_n, '(i1)') n

                    Me%Output%FloodPeriodFiles(n) = trim(adjustl(Me%Output%FloodPeriodFile))//&
                                                    "FloodPeriod_"//Char_n//".dat"   
                enddo

            else

                !if multiple flood period water column limits are not found 
                !use only one (the default it 0.05m)
                call GetData(Me%Output%FloodPeriodWaterColumnLimit,                     &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromFile,                                   &
                             keyword      = 'FLOOD_WATER_COLUMN_LIMIT',                 &
                             default      = 0.05,                                       &
                             ClientModule = 'ModuleRunOff',                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR798'

                Me%Output%FloodPeriodFile = trim(adjustl(Me%Output%FloodPeriodFile))//"FloodPeriod.dat"   

            endif

        endif

        if (Me%Output%WriteFloodArrivalTime) then

            call GetData(Me%Output%FloodArrivalWaterColumnLimit,                    &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'FLOOD_ARRIVAL_WATER_COLUMN_LIMIT',         &
                         default      = 0.05,                                       &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR799'    
         
        endif


        call ReadConvergenceParameters
        
        call ConstructTimeSeries
        
        call StartOutputBoxFluxes
                
        !Closes Data File
        call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR800'


    end subroutine ReadDataFile
    
    !--------------------------------------------------------------------------
    
    subroutine Read1DInteractionMapping(filename)
    
        !Arguments-------------------------------------------------------------

        character(len=PathLength)                   :: filename
        !Local----------------------------------------------------------------
        integer                                     :: mapping1DObjEnterData, ClientNumber, STAT_CALL
        integer                                     :: iflag
        type(T_NodeGridPoint), pointer              :: NewNodeGridPoint
        type(T_BankGridPoint), pointer              :: NewBankGridPoint, BankGridPointFlux
        type(T_MarginGridPoint), pointer            :: NewMarginGridPoint
        logical                                     :: BlockFound, FoundFlux
        !Begin----------------------------------------------------------------    
    
        
        mapping1DObjEnterData = 0

        call ConstructEnterData(mapping1DObjEnterData, filename, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR10'  

do1:    do         
            !Constructs Node Grid Point that will have 1D river level
            call ExtractBlockFromBuffer(mapping1DObjEnterData,                                  &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<BeginNodeGridPoint>',               &
                                        block_end       = '<EndNodeGridPoint>',                 &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
        
                
                allocate (NewNodeGridPoint)
                nullify(NewNodeGridPoint%Prev,NewNodeGridPoint%Next)                
                
                call GetData(NewNodeGridPoint%ID,                               &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'ID',                               &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR020'


                call GetData(NewNodeGridPoint%GridI,                            &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'GRID_I',                           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR030'       
            
                call GetData(NewNodeGridPoint%GridJ,                            &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'GRID_J',                           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR040'   
                
                call AddNodeGridPoint(NewNodeGridPoint)

                
            else
                
                call Block_Unlock(mapping1DObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR050'
                           
                exit do1   
                
            endif
            
        enddo do1




do2:    do         
            !Constructs Bank Grid Point that will recieve 1D river level from Node Grid Point
            call ExtractBlockFromBuffer(mapping1DObjEnterData,                                  &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<BeginBankGridPoint>',               &
                                        block_end       = '<EndBankGridPoint>',                 &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
        
                
                allocate (NewBankGridPoint)
                nullify(NewBankGridPoint%Prev,NewBankGridPoint%Next)                
                
                call GetData(NewBankGridPoint%ID,                               &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'ID',                               &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR060'


                call GetData(NewBankGridPoint%GridI,                            &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'GRID_I',                           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR070'       
            
                call GetData(NewBankGridPoint%GridJ,                            &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'GRID_J',                           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR080'   
                
                !link to NodeGridPoint
                call GetData(NewBankGridPoint%NGPId,                            &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'NGP_ID',                           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR090'                
                          
                call AddBankGridPoint(NewBankGridPoint)
                
            else
                
                call Block_Unlock(mapping1DObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR095'
                           
                exit do2  
                
            endif
            
        enddo do2




do3:    do         
            !Constructs Margin Grid Point that will have 1D river level interpolated from Bank River Points
            call ExtractBlockFromBuffer(mapping1DObjEnterData,                                  &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<BeginMarginGridPoint>',               &
                                        block_end       = '<EndMarginGridPoint>',                 &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
        
                
                allocate (NewMarginGridPoint)
                nullify(NewMarginGridPoint%Prev,NewMarginGridPoint%Next)                
                
                call GetData(NewMarginGridPoint%ID,                             &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'ID',                               &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR0100'


                call GetData(NewMarginGridPoint%GridI,                          &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'GRID_I',                           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR0110'       
            
                call GetData(NewMarginGridPoint%GridJ,                          &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'GRID_J',                           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR0120'   
                
                !link to BankGridPoint Upstream
                call GetData(NewMarginGridPoint%BGPUpId,                        &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'BGP_INTERPOLATION_1_ID',           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR0130'     
                
                !link to BankGridPoint Downstream
                call GetData(NewMarginGridPoint%BGPDownId,                      &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'BGP_INTERPOLATION_2_ID',           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR0140'        
                
                !interpolation fraction from Upstream for each MarginGridPoint
                call GetData(NewMarginGridPoint%InterpolationFraction,          &
                             mapping1DObjEnterData, iflag,                      &
                             SearchType   = FromBlock,                          &
                             keyword      = 'INTERPOLATION_FRACTION',           &
                             ClientModule = 'ModuleRunoff',                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR0150'                    
                
                if (NewMarginGridPoint%InterpolationFraction .lt. 0.0 .or. NewMarginGridPoint%InterpolationFraction .gt. 1.0) then
                    write (*,*)
                    write (*,*) 'Margin GridPoint INTERPOLATION_FRACTION can not be negative or higher than 1.0'
                    call SetError(FATAL_, KEYWORD_, "Read1DInteractionMapping - ModuleRunOff - ERR0155")
                endif
            
                !BankGripoint to send the computed flow depends on interpolation fraction
                !if less then 0.5 goes to upstream BGPm otherwise goes to downstream BGP
                if (NewMarginGridPoint%InterpolationFraction .le. 0.5) then
                    NewMarginGridPoint%BGPIntegrateFluxId = NewMarginGridPoint%BGPUpId                   
                else
                    NewMarginGridPoint%BGPIntegrateFluxId = NewMarginGridPoint%BGPDownId
                endif
                
                
                !dont allow duplicates margin points or fluxes would be duplicated
                call FindMarginGridPoint(NewMarginGridPoint%GridI, NewMarginGridPoint%GridJ, FoundFlux)
                if (FoundFlux) then
                    write(*,*)
                    write(*,*)'Found duplicate MarginGridPoint in cell ', NewMarginGridPoint%GridI, NewMarginGridPoint%GridJ
                    call SetError(FATAL_, KEYWORD_, "Read1DInteractionMapping - ModuleRunOff - ERR0156")
                endif                
                
                !associate i and j from BGP (DN) or NGP (OpenMI) to mgp to avoid searching in run-time
                !BGP or NGP where to associate flux
                call FindBankGridPoint(NewMarginGridPoint%BGPIntegrateFluxId, BankGridPointFlux, FoundFlux)
                if (.not. FoundFlux) then
                    write(*,*)
                    write(*,*)'Not found BankGridPoint to receive flux ', NewMarginGridPoint%BGPIntegrateFluxId
                    call SetError(FATAL_, KEYWORD_, "Read1DInteractionMapping - ModuleRunOff - ERR0157")
                else
                    
                    !in case DN flux will be interacted at Bank points. stop here
                    if (Me%ObjDrainageNetwork /= 0) then
                        
                        NewMarginGridPoint%GridIIntegrateFlux = BankGridPointFlux%GridI
                        NewMarginGridPoint%GridJIntegrateFlux = BankGridPointFlux%GridJ
                    
                    endif                    
                endif 
                
                 
                
                call AddMarginGridPoint(NewMarginGridPoint)
                
            else
                
                call Block_Unlock(mapping1DObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Read1DInteractionMapping - ModuleRunoff - ERR0160'
                           
                exit do3   
                
            endif
            
        enddo do3

        call AllocateRiverGridPointArrays()
    
    end subroutine Read1DInteractionMapping

    !--------------------------------------------------------------------------

    subroutine ReadIgnoredNodesFromFile(Filename)

        !Arguments-------------------------------------------------------------
        character(len=PathLength)                   :: Filename

        !Local----------------------------------------------------------------
        integer                                     :: IgnoredNodesObjEnterData, ClientNumber, STAT_CALL
        integer                                     :: iflag, nLine, nNode
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine
        !Begin----------------------------------------------------------------    
        
        IgnoredNodesObjEnterData  = 0
        Me%NumberOfIgnoredNodes   = 0

        call ConstructEnterData(IgnoredNodesObjEnterData, Filename, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadIgnoredNodesFromFile - ModuleRunoff - ERR01'  
        

        !Count number of ignored nodes
        call ExtractBlockFromBuffer(IgnoredNodesObjEnterData,                               &
                                    ClientNumber    = ClientNumber,                         &
                                    block_begin     = '<begin_ignored_nodes>',              &
                                    block_end       = '<end_ignored_nodes>',                &
                                    BlockFound      = BlockFound,                           &
                                    FirstLine       = FirstLine,                            &
                                    LastLine        = LastLine,                             &  
                                    STAT            = STAT_CALL)
        if(STAT_CALL .EQ. SUCCESS_)then
            
            if (BlockFound) then

                Me%NumberOfIgnoredNodes = LastLine - FirstLine - 1

                if(Me%NumberOfIgnoredNodes > 0)then

                    allocate(Me%IgnoredNodes(1:Me%NumberOfIgnoredNodes))
            
                    nNode = 1

                    do nLine = FirstLine + 1, LastLine - 1
                    
                        call GetData(Me%IgnoredNodes(nNode), IgnoredNodesObjEnterData, &
                                     iflag, Buffer_Line = nLine, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ReadIgnoredNodesFromFile - ModuleRunoff - ERR02'

                        nNode = nNode + 1 
                    enddo 

                else

                    write(*,*)
                    write(*,*)"No SWMM ignored nodes were defined."
                    write(*,*)"All SWMM nodes inside the grid will be linked to 2D model."
                    write(*,*)"File: ", trim(adjustl(Filename))
                    write(*,*)"ReadIgnoredNodesFromFile - ModuleRunoff - WRN01"
                    write(*,*)

                endif
                
                call Block_Unlock(IgnoredNodesObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIgnoredNodesFromFile - ModuleRunoff - ERR03'  

            else
                
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                write(*,*) 'Could not find <begin_ignored_nodes>...<end_ignored_nodes> block'
                stop    'ReadIgnoredNodesFromFile - ModuleRunoff - ERR04'

            endif

        elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
            write(*,*)  
            write(*,*) 'Error calling ExtractBlockFromBuffer. '
            stop    'ReadIgnoredNodesFromFile - ModuleRunoff - ERR05'
        
        else
            stop    'ReadIgnoredNodesFromFile - ModuleRunoff - ERR06'
        end if
            
        call KillEnterData(IgnoredNodesObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadIgnoredNodesFromFile - ModuleRunoff - ERR10'  



    end subroutine ReadIgnoredNodesFromFile

    !--------------------------------------------------------------------------

    subroutine ReadPondsFromFile(Filename)
    
        !Arguments-------------------------------------------------------------
        character(len=PathLength)                   :: Filename

        !Local----------------------------------------------------------------
        integer                                     :: PondsObjEnterData, ClientNumber, STAT_CALL
        integer                                     :: iflag, n
        logical                                     :: BlockFound

        !Begin----------------------------------------------------------------    
        
        PondsObjEnterData  = 0
        Me%NumberOfPonds   = 0

        call ConstructEnterData(PondsObjEnterData, Filename, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadPondsFromFile - ModuleRunoff - ERR01'  
        

do1:    do         
            !Count number of ponds
            call ExtractBlockFromBuffer(PondsObjEnterData,                                      &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_pond>',                       &
                                        block_end       = '<end_pond>',                         &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
                
                Me%NumberOfPonds = Me%NumberOfPonds + 1
                
            else
                
                call Block_Unlock(PondsObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadPondsFromFile - ModuleRunoff - ERR02'                             
                
                exit do1   
                
            endif
            
        enddo do1

        if(Me%NumberOfPonds < 1)then
            write(*,*)
            write(*,*)"No ponds were defined."
            write(*,*)"File: ", trim(adjustl(Filename))
            write(*,*)"ReadPondsFromFile - ModuleRunoff - WRN01"
            write(*,*)
            return
        endif

        allocate(Me%Ponds(1:Me%NumberOfPonds))
        
        call RewindBuffer (PondsObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadPondsFromFile - ModuleRunOff - ERR03'
        
        n = 0
        
do2:    do         
            !Count number of ponds
            call ExtractBlockFromBuffer(PondsObjEnterData,                                      &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_pond>',                       &
                                        block_end       = '<end_pond>',                         &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
                
                n = n + 1

                Me%Ponds(n)%ID = n

                call GetData(Me%Ponds(n)%Name,                                                  &
                             PondsObjEnterData, iflag,                                          &
                             Keyword        = 'NAME',                                           &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             STAT           = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadPondsFromFile - ModuleRunOff - ERR04'

                if(iflag == 0)then
                    write(*,*)"Please define NAME for pond"
                    write(*,*)"ID = ", n
                    stop 'ReadPondsFromFile - ModuleRunOff - ERR05'
                endif

            else
                
                call Block_Unlock(PondsObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadPondsFromFile - ModuleRunoff - ERR100'                             
                
                exit do2   
                
            endif
            
        enddo do2

        call KillEnterData(PondsObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadPondsFromFile - ModuleRunoff - ERR900'  
    
    end subroutine ReadPondsFromFile

    !--------------------------------------------------------------------------

    subroutine ReadOpenChannelLinksFromFile(Filename)
    
        !Arguments-------------------------------------------------------------
        character(len=PathLength)                   :: Filename

        !Local----------------------------------------------------------------
        integer                                     :: OCLinksObjEnterData, ClientNumber, STAT_CALL
        integer                                     :: iflag, n, i, j
        logical                                     :: BlockFound

        !Begin----------------------------------------------------------------    
        
        OCLinksObjEnterData         = 0
        Me%NumberOfOpenChannelLinks = 0

        call ConstructEnterData(OCLinksObjEnterData, Filename, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunoff - ERR01'  
        
do1:    do         
            !Count number of ponds
            call ExtractBlockFromBuffer(OCLinksObjEnterData,                                    &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_openchannel_link>',           &
                                        block_end       = '<end_openchannel_link>',             &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
                
                Me%NumberOfOpenChannelLinks = Me%NumberOfOpenChannelLinks + 1
                
            else
                
                call Block_Unlock(OCLinksObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunoff - ERR10'                             
                
                exit do1   
                
            endif
            
        enddo do1

        if(Me%NumberOfOpenChannelLinks < 1)then
            write(*,*)
            write(*,*)"No open channel links were defined."
            write(*,*)"File: ", trim(adjustl(Filename))
            write(*,*)"ReadOpenChannelLinksFromFile - ModuleRunoff - WRN01"
            write(*,*)
            return
        endif

        allocate(Me%OpenChannelLinks(1:Me%NumberOfOpenChannelLinks))
        
        call RewindBuffer (OCLinksObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR20'

        n = 0
        
do2:    do         
            !Count number of ponds
            call ExtractBlockFromBuffer(OCLinksObjEnterData,                                    &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_openchannel_link>',           &
                                        block_end       = '<end_openchannel_link>',             &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
                
                n = n + 1

                Me%OpenChannelLinks(n)%ID = n

                call GetData(Me%OpenChannelLinks(n)%TypeOf,                                     &
                             OCLinksObjEnterData, iflag,                                        &
                             Keyword        = 'TYPE',                                           &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR30'

                if(iflag == 0)then
                    write(*,*)"Please define TYPE for open channel link"
                    write(*,*)"ID = ", n
                    stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR40'
                endif

                if(Me%OpenChannelLinks(n)%TypeOf .ne. Direct_       .and. &
                   Me%OpenChannelLinks(n)%TypeOf .ne. Weighted_     .and. &
                   Me%OpenChannelLinks(n)%TypeOf .ne. OutfallLink_  .and. &
                   Me%OpenChannelLinks(n)%TypeOf .ne. PondLink_)then
                    write(*,*)"Unknown TYPE for open channel link"
                    write(*,*)"ID = ", n
                    stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR50'
                endif

                call GetData(Me%OpenChannelLinks(n)%I,                                          &
                             OCLinksObjEnterData, iflag,                                        &
                             Keyword        = 'GRID_I',                                         &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR60'
                
                if(iflag == 0)then
                    write(*,*)"Please define GRID_I for open channel link"
                    write(*,*)"ID = ", n
                    stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR70'
                endif

                call GetData(Me%OpenChannelLinks(n)%J,                                          &
                             OCLinksObjEnterData, iflag,                                        &
                             Keyword        = 'GRID_J',                                         &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR80'

                if(iflag == 0)then
                    write(*,*)"Please define GRID_J for open channel link"
                    write(*,*)"ID = ", n
                    stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR90'
                endif

                if(Me%ExtVar%BasinPoints(Me%OpenChannelLinks(n)%I, Me%OpenChannelLinks(n)%J) .ne. 1)then
                    write(*,*)"Open channel link is not located in active grid point"
                    write(*,*)"ID     = ", n
                    write(*,*)"GRID_I = ", Me%OpenChannelLinks(n)%I
                    write(*,*)"GRID_J = ", Me%OpenChannelLinks(n)%J
                    stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR91'
                end if  


                i = Me%OpenChannelLinks(n)%I
                j = Me%OpenChannelLinks(n)%J

                Me%OpenChannelLinks(n)%CellWidth = (Me%ExtVar%DUX(i, j) + Me%ExtVar%DVY(i, j) ) / 2.0

                !Me%OpenChannelLinks(n)%FluxWidth = (Me%ExtVar%BasinPoints(i,j-1) + &
                !                                    Me%ExtVar%BasinPoints(i,j+1) + &
                !                                    Me%ExtVar%BasinPoints(i-1,j) + &
                !                                    Me%ExtVar%BasinPoints(i+1,j))* &
                !                                   Me%OpenChannelLinks(n)%CellWidth  !SOBRINHO
                
                Me%OpenChannelLinks(n)%FluxWidth = ((1-Me%ExtVar%BasinPoints(i,j-1)) + &
                                                    (1-Me%ExtVar%BasinPoints(i,j+1)) + &
                                                    (1-Me%ExtVar%BasinPoints(i-1,j)) + &
                                                    (1-Me%ExtVar%BasinPoints(i+1,j)))* &
                                                   Me%OpenChannelLinks(n)%CellWidth


                call GetData(Me%OpenChannelLinks(n)%LinkNodeName,                               &
                             OCLinksObjEnterData, iflag,                                        &
                             Keyword        = 'LINK_NAME',                                      &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR100'

                if(iflag == 0)then
                    write(*,*)"Please define LINK_NAME for open channel link"
                    write(*,*)"ID = ", n
                    write(*,*)"GRID_I = ", Me%OpenChannelLinks(n)%I
                    write(*,*)"GRID_J = ", Me%OpenChannelLinks(n)%J
                    stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR110'
                endif

                if(Me%OpenChannelLinks(n)%TypeOf == Weighted_)then

                    call GetData(Me%OpenChannelLinks(n)%SecondLinkNodeName,                     &
                                 OCLinksObjEnterData, iflag,                                    &
                                 Keyword        = 'SECOND_LINK_NAME',                           &
                                 SearchType     = FromBlock,                                    &
                                 ClientModule   ='ModuleRunoff',                                &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR120'

                    if(iflag == 0)then
                        write(*,*)"Please define SECOND_LINK_NAME for open channel link"
                        write(*,*)"ID = ", n
                        write(*,*)"GRID_I = ", Me%OpenChannelLinks(n)%I
                        write(*,*)"GRID_J = ", Me%OpenChannelLinks(n)%J
                        stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR130'
                    endif

                endif

                if(Me%OpenChannelLinks(n)%TypeOf == Weighted_ .or. Me%OpenChannelLinks(n)%TypeOf == OutfallLink_)then
                    
                    call GetData(Me%OpenChannelLinks(n)%Weight,                                    &
                                    OCLinksObjEnterData, iflag,                                    &
                                    Keyword        = 'WEIGHT',                                     &
                                    SearchType     = FromBlock,                                    &
                                    ClientModule   ='ModuleRunoff',                                &
                                    Default        = 1.0,                                          &
                                    STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR140'

                    if(iflag == 0)then
                        write(*,*)"Please define WEIGHT for open channel link"
                        write(*,*)"ID     = ", n
                        write(*,*)"NAME   = ", trim(adjustl(Me%OpenChannelLinks(n)%LinkNodeName))
                        write(*,*)"GRID_I = ", Me%OpenChannelLinks(n)%I
                        write(*,*)"GRID_J = ", Me%OpenChannelLinks(n)%J
                        stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR150'
                    endif

                    if(Me%OpenChannelLinks(n)%Weight < 0.0)then
                        write(*,*)"WEIGHT for open channel link is lower than 0.0"
                        write(*,*)"ID     = ", n
                        write(*,*)"NAME   = ", trim(adjustl(Me%OpenChannelLinks(n)%LinkNodeName))
                        write(*,*)"GRID_I = ", Me%OpenChannelLinks(n)%I
                        write(*,*)"GRID_J = ", Me%OpenChannelLinks(n)%J
                        stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR151'
                    endif

                    if(Me%OpenChannelLinks(n)%Weight > 1.0)then
                        write(*,*)"WEIGHT for open channel link is higher than 1.0"
                        write(*,*)"ID     = ", n
                        write(*,*)"NAME   = ", trim(adjustl(Me%OpenChannelLinks(n)%LinkNodeName))
                        write(*,*)"GRID_I = ", Me%OpenChannelLinks(n)%I
                        write(*,*)"GRID_J = ", Me%OpenChannelLinks(n)%J
                        stop 'ReadOpenChannelLinksFromFile - ModuleRunOff - ERR152'
                    endif

                endif

            else
                
                call Block_Unlock(OCLinksObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunoff - ERR200'                             
                
                exit do2   
                
            endif
            
        enddo do2


        call KillEnterData(OCLinksObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadOpenChannelLinksFromFile - ModuleRunoff - ERR900'  


    end subroutine ReadOpenChannelLinksFromFile

    !--------------------------------------------------------------------------
    
    subroutine ReadHeadwallsFromFile(Filename, RootSRT)
    
        !Arguments-------------------------------------------------------------
        character(len=PathLength), intent(in)       :: Filename
        character(len=PathLength), intent(in)       :: RootSRT

        !Local----------------------------------------------------------------
        integer                                     :: HeadwallsObjEnterData, ClientNumber, STAT_CALL
        integer                                     :: iflag, n
        logical                                     :: BlockFound
        character(len=PathLength)                   :: HeadwallOutputFilename

        !Begin----------------------------------------------------------------    
        
        HeadwallsObjEnterData  = 0
        Me%NumberOfHeadwalls   = 0

        call ConstructEnterData(HeadwallsObjEnterData, Filename, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunoff - ERR01'  
        

do1:    do         
            !Count number of headwalls
            call ExtractBlockFromBuffer(HeadwallsObjEnterData,                                  &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_headwall>',                   &
                                        block_end       = '<end_headwall>',                     &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
                
                Me%NumberOfHeadwalls = Me%NumberOfHeadwalls + 1
                
            else
                
                call Block_Unlock(HeadwallsObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunoff - ERR10'                             
                
                exit do1   
                
            endif
            
        enddo do1

        if(Me%NumberOfHeadwalls < 1)then
            write(*,*)
            write(*,*)"No headwalls were defined."
            write(*,*)"File: ", trim(adjustl(Filename))
            write(*,*)"ReadHeadwallsFromFile - ModuleRunoff - WRN01"
            write(*,*)
            return
        endif

        allocate(Me%Headwalls(1:Me%NumberOfHeadwalls))
        
        call RewindBuffer (HeadwallsObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunOff - ERR20'
        
        n = 0
        
do2:    do         
            !Count number of headwalls
            call ExtractBlockFromBuffer(HeadwallsObjEnterData,                                  &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_headwall>',                   &
                                        block_end       = '<end_headwall>',                     &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
                
                n = n + 1

                Me%Headwalls(n)%ID = n

                call GetData(Me%Headwalls(n)%Name,                                              &
                             HeadwallsObjEnterData, iflag,                                      &
                             Keyword        = 'NAME',                                           &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             STAT           = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunOff - ERR30'

                if(iflag == 0)then
                    write(*,*)"Please define NAME for headwall"
                    write(*,*)"ID = ", n
                    stop 'ReadHeadwallsFromFile - ModuleRunOff - ERR40'
                endif

                call GetData(Me%Headwalls(n)%OutputResults,                                     &
                             HeadwallsObjEnterData, iflag,                                      &
                             Keyword        = 'OUTPUT_RESULTS',                                 &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             Default        = .false.,                                          &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunOff - ERR40'

                if(Me%Headwalls(n)%OutputResults)then
                
                    call GetData(Me%Headwalls(n)%OutputTimeStep,                                &
                                 HeadwallsObjEnterData, iflag,                                  &
                                 Keyword        = 'DT_OUTPUT_TIME',                             &
                                 SearchType     = FromBlock,                                    &
                                 ClientModule   ='ModuleRunoff',                                &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunOff - ERR46'
                
                    if(iflag == 0)then
                        write(*,*)"Missing DT_OUTPUT_TIME in inlet ", trim(adjustl(Me%Headwalls(n)%Name))
                        write(*,*)"in file ", trim(adjustl(Filename))
                        stop 'ReadHeadwallsFromFile - ModuleRunOff - ERR47'
                    endif

                    Me%Headwalls(n)%OutputTime     = 0.0
                    Me%Headwalls(n)%NextOutputTime = Me%BeginTime
                
                    call UnitsManager(Me%Headwalls(n)%OutputUnit, OPEN_FILE, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunOff - ERR48'

                    HeadwallOutputFilename = trim(adjustl(RootSRT))//trim(adjustl(Me%Headwalls(n)%Name))//".srh"
                
                    open(Unit = Me%Headwalls(n)%OutputUnit,                                     &
                         File = trim(adjustl(HeadwallOutputFilename)),                          &
                         STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunOff - ERR49'

                    call WriteDataLine(Me%Headwalls(n)%OutputUnit, "NAME", Me%Headwalls(n)%Name)
                    call WriteDataLine(Me%Headwalls(n)%OutputUnit, 'SERIE_INITIAL_DATA', Me%BeginTime)
                    call WriteDataLine(Me%Headwalls(n)%OutputUnit, 'TIME_UNITS', 'SECONDS')
                    
                    !call WriteDataLine(Me%Headwalls(n)%OutputUnit, 'time WL2D_1 WL2D_2 WL1D_1 WL1D_2 Depth2D_1 Depth2D_2 Depth1D_1 Depth1D_2 flow flowEnteringCell')
                    call WriteDataLine(Me%Headwalls(n)%OutputUnit, 'time water_level water_depth flow_entering_cell')
                    call WriteDataLine(Me%Headwalls(n)%OutputUnit, '<BeginTimeSerie>')

                endif

            else
                
                call Block_Unlock(HeadwallsObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunoff - ERR150'                             
                
                exit do2   
                
            endif
            
        enddo do2


        call KillEnterData(HeadwallsObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadHeadwallsFromFile - ModuleRunoff - ERR200'  
    
    end subroutine ReadHeadwallsFromFile

    !--------------------------------------------------------------------------
    
    subroutine ReadInletsFromFile(Filename, RootSRT)
    
        !Arguments-------------------------------------------------------------

        character(len=PathLength), intent(in)       :: Filename
        character(len=PathLength), intent(in)       :: RootSRT

        !Local----------------------------------------------------------------
        integer                                     :: InletsObjEnterData, ClientNumber, STAT_CALL
        integer                                     :: iflag, n
        logical                                     :: BlockFound
        real                                        :: InletFlowBelowMinStage, InletFlowAboveMaxStage
        character(len=PathLength)                   :: InletOutputFilename

        !Begin----------------------------------------------------------------    
        
        InletsObjEnterData  = 0
        Me%NumberOfInlets   = 0

        call ConstructEnterData(InletsObjEnterData, Filename, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInletsFromFile - ModuleRunoff - ERR01'  
        

do1:    do         
            !Count number of inlets
            call ExtractBlockFromBuffer(InletsObjEnterData,                                     &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_inlet>',                      &
                                        block_end       = '<end_inlet>',                        &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
                
                Me%NumberOfInlets = Me%NumberOfInlets + 1
                
            else
                
                call Block_Unlock(InletsObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadInletsFromFile - ModuleRunoff - ERR10'                             
                
                exit do1   
                
            endif
            
        enddo do1

        if(Me%NumberOfInlets < 1)then
            write(*,*)
            write(*,*)"No inlets were defined."
            write(*,*)"File: ", trim(adjustl(Filename))
            write(*,*)"ReadInletsFromFile - ModuleRunoff - WRN01"
            write(*,*)
            return
        endif

        allocate(Me%Inlets(1:Me%NumberOfInlets))
        
        call RewindBuffer (InletsObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR20'
        
        n = 0
        
do2:    do         
            !Count number of inlets
            call ExtractBlockFromBuffer(InletsObjEnterData,                                     &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_inlet>',                      &
                                        block_end       = '<end_inlet>',                        &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
                
                n = n + 1

                Me%Inlets(n)%ID = n

                call GetData(Me%Inlets(n)%Name,                                                 &
                             InletsObjEnterData, iflag,                                         &
                             Keyword        = 'NAME',                                           &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             STAT           = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR30'

                if(iflag == 0)then
                    write(*,*)"Please define NAME for inlet"
                    write(*,*)"ID = ", n
                    stop 'ReadInletsFromFile - ModuleRunOff - ERR05'
                endif
                
                call GetData(Me%Inlets(n)%TypeOf,                                               &
                             InletsObjEnterData, iflag,                                         &
                             Keyword        = 'INLET_TYPE',                                     &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             Default        = Weir_,                                            &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR40'

                call GetData(Me%Inlets(n)%OutputResults,                                        &
                             InletsObjEnterData, iflag,                                         &
                             Keyword        = 'OUTPUT_RESULTS',                                 &
                             SearchType     = FromBlock,                                        &
                             ClientModule   ='ModuleRunoff',                                    &
                             Default        = .false.,                                          &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR45'

                if(Me%Inlets(n)%OutputResults)then
                
                    call GetData(Me%Inlets(n)%OutputTimeStep,                                   &
                                 InletsObjEnterData, iflag,                                     &
                                 Keyword        = 'DT_OUTPUT_TIME',                             &
                                 SearchType     = FromBlock,                                    &
                                 ClientModule   ='ModuleRunoff',                                &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR46'
                
                    if(iflag == 0)then
                        write(*,*)"Missing DT_OUTPUT_TIME in inlet ", trim(adjustl(Me%Inlets(n)%Name))
                        write(*,*)"in file ", trim(adjustl(Filename))
                        stop 'ReadInletsFromFile - ModuleRunOff - ERR47'
                    endif

                    Me%Inlets(n)%OutputTime     = 0.0
                    Me%Inlets(n)%NextOutputTime = Me%BeginTime
                
                    call UnitsManager(Me%Inlets(n)%OutputUnit, OPEN_FILE, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR48'

                    InletOutputFilename = trim(adjustl(RootSRT))//trim(adjustl(Me%Inlets(n)%Name))//".sri"
                
                    open(Unit = Me%Inlets(n)%OutputUnit,                                        &
                         File = trim(adjustl(InletOutputFilename)),                             &
                         STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR49'

                    call WriteDataLine(Me%Inlets(n)%OutputUnit, "NAME", Me%Inlets(n)%Name)
                    call WriteDataLine(Me%Inlets(n)%OutputUnit, 'SERIE_INITIAL_DATA', Me%BeginTime)
                    call WriteDataLine(Me%Inlets(n)%OutputUnit, 'TIME_UNITS', 'SECONDS')
                    
                    call WriteDataLine(Me%Inlets(n)%OutputUnit, 'time flow_entering_cell potential_flow effective_flow')
                    call WriteDataLine(Me%Inlets(n)%OutputUnit, '<BeginTimeSerie>')

                endif
                
                if(Me%Inlets(n)%TypeOf .eq. Weir_)then

                    call GetData(Me%Inlets(n)%Width,                                                &
                                 InletsObjEnterData, iflag,                                         &
                                 Keyword        = 'WIDTH',                                          &
                                 SearchType     = FromBlock,                                        &
                                 ClientModule   ='ModuleRunoff',                                    &
                                 STAT           = STAT_CALL)                                      
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR50'
                    
                    if(iflag == 0)then
                        write(*,*)"No WIDTH was specified for ", trim(adjustl(Me%Inlets(n)%Name))
                        stop 'ReadInletsFromFile - ModuleRunOff - ERR60'
                    endif 
                    
                    
                elseif(Me%Inlets(n)%TypeOf .eq. FlowCapture_)then
                    
                    call GetData(Me%Inlets(n)%CaptureFraction,                                      &
                                 InletsObjEnterData, iflag,                                         &
                                 Keyword        = 'FLOW_CAPTURE_FRACTION',                          &
                                 SearchType     = FromBlock,                                        &
                                 Default        = 1.0,                                              &
                                 ClientModule   ='ModuleRunoff',                                    &
                                 STAT           = STAT_CALL)                                      
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR70'

                    if(Me%Inlets(n)%CaptureFraction < 0.0 .or. Me%Inlets(n)%CaptureFraction > 1.0)then
                        write(*,*)"Inlet FLOW_CAPTURE_FRACTION is not valid: ", trim(adjustl(Me%Inlets(n)%Name))
                        stop 'ReadInletsFromFile - ModuleRunOff - ERR71'
                    endif

                elseif(Me%Inlets(n)%TypeOf .eq. DepthFlowRatingCurve_ .or. Me%Inlets(n)%TypeOf .eq. FlowFlowRatingCurve_)then

                    call GetData(Me%Inlets(n)%RatingCurveFileName,                                  &
                                 InletsObjEnterData, iflag,                                         &
                                 Keyword        = 'RATING_CURVE_FILENAME',                          &
                                 SearchType     = FromBlock,                                        &
                                 ClientModule   ='ModuleRunoff',                                    &
                                 STAT           = STAT_CALL)                                      
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR80'
                    
                    if(iflag == 0)then
                        write(*,*)"No RATING_CURVE_FILENAME was specified for ", trim(adjustl(Me%Inlets(n)%Name))
                        stop 'ReadInletsFromFile - ModuleRunOff - ERR90'
                    endif 

                    call ReadRatingCurveFile(Me%Inlets(n))

                    InletFlowBelowMinStage = Me%Inlets(n)%RatingCurveFlow(1)
                    InletFlowAboveMaxStage = Me%Inlets(n)%RatingCurveFlow(Me%Inlets(n)%RatingCurve_nValues)

                    call GetData(Me%Inlets(n)%RatingCurveBelowMin,                                  &
                                 InletsObjEnterData, iflag,                                         &
                                 Keyword        = 'FLOW_BELOW_MIN_STAGE',                           &
                                 SearchType     = FromBlock,                                        &
                                 ClientModule   = 'ModuleRunoff',                                   &
                                 Default        = InletFlowBelowMinStage,                           &
                                 STAT           = STAT_CALL)                                      
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadRatingCurvesFromFile - ModuleRunOff - ERR100'
                    
                    call GetData(Me%Inlets(n)%RatingCurveAboveMax,                                  &
                                 InletsObjEnterData, iflag,                                         &
                                 Keyword        = 'FLOW_ABOVE_MAX_STAGE',                           &
                                 SearchType     = FromBlock,                                        &
                                 ClientModule   = 'ModuleRunoff',                                   &
                                 Default        = InletFlowAboveMaxStage,                           &
                                 STAT           = STAT_CALL)                                      
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInletsFromFile - ModuleRunOff - ERR120'
                
                else
                    
                    write(*,*)"Invalid INLET_TYPE"
                    stop 'ReadInletsFromFile - ModuleRunOff - ERR140'
                    
                endif
               
            else
                
                call Block_Unlock(InletsObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadInletsFromFile - ModuleRunoff - ERR150'                             
                
                exit do2   
                
            endif
            
        enddo do2


        call KillEnterData(InletsObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInletsFromFile - ModuleRunoff - ERR200'  
    
    end subroutine ReadInletsFromFile
    
    !--------------------------------------------------------------------------

    subroutine ReadRatingCurveFile(Inlet)
    
        !Arguments-------------------------------------------------------------
        type(T_SewerGEMSInlet)                      :: Inlet

        !Local----------------------------------------------------------------
        integer                                     :: RatingCurveObjEnterData, ClientNumber, STAT_CALL
        integer                                     :: iflag, FirstLine, LastLine, iValue, iLine
        logical                                     :: BlockFound
        real, dimension(:), pointer                 :: BufferLine

        !Begin----------------------------------------------------------------

        call ConstructEnterData(RatingCurveObjEnterData, Inlet%RatingCurveFilename, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadRatingCurveFile - ModuleRunoff - ERR01'  
            
        !Get Block with rating curves values
        call ExtractBlockFromBuffer(RatingCurveObjEnterData, ClientNumber,                      &
                                    "<begin_rating_curve_table>", "<end_rating_curve_table>",   &
                                    BlockFound, FirstLine = FirstLine, LastLine = LastLine,     &
                                    STAT      = STAT_CALL)

        if (STAT_CALL .EQ. SUCCESS_  .and. BlockFound) then

            Inlet%RatingCurve_nValues = LastLine - FirstLine - 1
            
            allocate(Inlet%RatingCurveStage(1:Inlet%RatingCurve_nValues))
            allocate(Inlet%RatingCurveFlow (1:Inlet%RatingCurve_nValues))
                    
            allocate(BufferLine(2))

            iValue = 1;
            do  iLine = FirstLine+1, LastLine-1
                call GetData(BufferLine,                                &
                             RatingCurveObjEnterData,                   &
                             iflag, Buffer_Line = iLine,                &
                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .or. iflag /= 2)stop 'ReadRatingCurveFile - ModuleRunoff - ERR10'  
                    
                Inlet%RatingCurveStage(iValue) = BufferLine (1)
                Inlet%RatingCurveFlow (iValue) = BufferLine (2)
                iValue = iValue + 1
            end do

            deallocate (BufferLine)

            call Block_Unlock(RatingCurveObjEnterData, ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadRatingCurveFile - ModuleRunoff - ERR20'  
                    
            call KillEnterData(RatingCurveObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadRatingCurveFile - ModuleRunoff - ERR30'  

        else

            stop 'ReadRatingCurveFile - ModuleRunoff - ERR40'

        endif 

    end subroutine ReadRatingCurveFile

    !--------------------------------------------------------------------------

    ! This subroutine adds a new NodeGridPoint to the List  
    subroutine AddNodeGridPoint(NewNodeGridPoint)

        !Arguments-------------------------------------------------------------
        type(T_NodeGridPoint),              pointer     :: NewNodeGridPoint

        !----------------------------------------------------------------------

        ! Add to the NodeGridPoint List a new NodeGridPoint
        if (.not.associated(Me%FirstNodeGridPoint)) then
            Me%NodeGridPointNumber     = 1
            Me%FirstNodeGridPoint        => NewNodeGridPoint
            Me%LastNodeGridPoint         => NewNodeGridPoint
        else
            NewNodeGridPoint%Prev        => Me%LastNodeGridPoint
            Me%LastNodeGridPoint%Next    => NewNodeGridPoint
            Me%LastNodeGridPoint         => NewNodeGridPoint
            Me%NodeGridPointNumber     = Me%NodeGridPointNumber + 1
        end if 


    end subroutine AddNodeGridPoint 
    
    !--------------------------------------------------------------------------

    ! This subroutine adds a new BankGridPoint to the List  
    subroutine AddBankGridPoint(NewBankGridPoint)

        !Arguments-------------------------------------------------------------
        type(T_BankGridPoint),              pointer     :: NewBankGridPoint

        !----------------------------------------------------------------------

        ! Add to the BankGridPoint List a new BankGridPoint
        if (.not.associated(Me%FirstBankGridPoint)) then
            Me%BankGridPointNumber     = 1
            Me%FirstBankGridPoint        => NewBankGridPoint
            Me%LastBankGridPoint         => NewBankGridPoint
        else
            NewBankGridPoint%Prev        => Me%LastBankGridPoint
            Me%LastBankGridPoint%Next    => NewBankGridPoint
            Me%LastBankGridPoint         => NewBankGridPoint
            Me%BankGridPointNumber     = Me%BankGridPointNumber + 1
        end if 


    end subroutine AddBankGridPoint 
    
    !--------------------------------------------------------------------------

    ! This subroutine adds a new MarginGridPoint to the List  
    subroutine AddMarginGridPoint(NewMarginGridPoint)

        !Arguments-------------------------------------------------------------
        type(T_MarginGridPoint),              pointer     :: NewMarginGridPoint

        !----------------------------------------------------------------------

        ! Add to the MarginGridPoint List a new MarginGridPoint
        if (.not.associated(Me%FirstMarginGridPoint)) then
            Me%MarginGridPointNumber     = 1
            Me%FirstMarginGridPoint        => NewMarginGridPoint
            Me%LastMarginGridPoint         => NewMarginGridPoint
        else
            NewMarginGridPoint%Prev        => Me%LastMarginGridPoint
            Me%LastMarginGridPoint%Next    => NewMarginGridPoint
            Me%LastMarginGridPoint         => NewMarginGridPoint
            Me%MarginGridPointNumber     = Me%MarginGridPointNumber + 1
        end if 


    end subroutine AddMarginGridPoint   
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Allocates and fills river grid point array pointer structures for O(n) access
    !---------------------------------------------------------------------------
    subroutine AllocateRiverGridPointArrays()
    integer :: i, j
    logical :: found1, found2
    type(T_NodeGridPoint), pointer :: CurrNodeGridPoint
    type(T_MarginGridPoint), pointer :: CurrMarginGridPoint
    type(T_BankGridPoint), pointer :: CurrBankGridPoint

    allocate(Me%NodeGridPointArray(Me%NodeGridPointNumber))
    allocate(Me%MarginGridPointArray(Me%MarginGridPointNumber))
    allocate(Me%BankGridPointArray(Me%BankGridPointNumber))
    
    !print*, 'number of NodeGridPoint objects = ', size(Me%NodeGridPointArray)
    !print*, 'number of MarginGridPoint objects = ', size(Me%MarginGridPointArray)
    !print*, 'number of BankGridPoint objects = ', size(Me%BankGridPointArray)
    
    !filling the arrays
    CurrNodeGridPoint => Me%FirstNodeGridPoint
    do i=1, Me%NodeGridPointNumber
        Me%NodeGridPointArray(i)%ptr => CurrNodeGridPoint
        CurrNodeGridPoint => CurrNodeGridPoint%Next
    end do
    
    CurrMarginGridPoint => Me%FirstMarginGridPoint
    do i=1, Me%MarginGridPointNumber
        Me%MarginGridPointArray(i)%ptr => CurrMarginGridPoint
        CurrMarginGridPoint => CurrMarginGridPoint%Next
    end do
    
    CurrBankGridPoint => Me%FirstBankGridPoint
    do i=1, Me%BankGridPointNumber
        Me%BankGridPointArray(i)%ptr => CurrBankGridPoint
        CurrBankGridPoint => CurrBankGridPoint%Next
    end do
    
    !indexing the correct array position in the objects on the list
    CurrBankGridPoint => Me%FirstBankGridPoint
    do i=1, Me%BankGridPointNumber
        found1 = .false.
        do j=1, Me%NodeGridPointNumber
            if (CurrBankGridPoint%NGPId == Me%NodeGridPointArray(j)%ptr%ID) then
                CurrBankGridPoint%NGPIDidx = j
                found1 = .true.
                exit
            end if
        end do
        if (.not.found1) stop 'ModuleRunOff::AllocateRiverGridPointArrays - node grid point corresponding to bank grid point not found'
        CurrBankGridPoint => CurrBankGridPoint%Next
    end do
    
    CurrMarginGridPoint => Me%FirstMarginGridPoint
    do i=1, Me%MarginGridPointNumber
        found1 = .false.
        found2 = .false.
        do j=1, Me%BankGridPointNumber
            if (CurrMarginGridPoint%BGPUpId == Me%BankGridPointArray(j)%ptr%ID) then
                CurrMarginGridPoint%BGPUpIdidx = j
                found1 = .true.
            end if
            if (CurrMarginGridPoint%BGPDownId == Me%BankGridPointArray(j)%ptr%ID) then
                CurrMarginGridPoint%BGPDownIdidx = j
                found2 = .true.
            end if
        end do
        if (.not.found1) stop 'ModuleRunOff::AllocateRiverGridPointArrays - upper node grid point corresponding to bank grid point not found'
        if (.not.found2) stop 'ModuleRunOff::AllocateRiverGridPointArrays - lower node grid point corresponding to bank grid point not found'
        CurrMarginGridPoint => CurrMarginGridPoint%Next
    end do    
       
    end subroutine AllocateRiverGridPointArrays

    !--------------------------------------------------------------------------    
   
    subroutine StartOutputBoxFluxes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer :: STAT_CALL
        integer :: iflag
        logical :: exist, opened
        character(len=StringLength), dimension(:),  pointer :: FluxesOutputList
        character(len=StringLength), dimension(:),  pointer :: ScalarOutputList

        !Begin-----------------------------------------------------------------
       
        ! This keyword have two functions if exist fluxes between boxes are compute 
        ! and the value read is the name file where the boxes are defined
        call GetData(Me%Files%BoxesFile,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'BOXFLUXES',                                      &
                     SearchType     = FromFile,                                         &
                     ClientModule   ='ModuleRunoff',                                    &
                     STAT           = STAT_CALL)                                      

        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Subroutine StartOutputBoxFluxes - ModuleRunoff. ERR02.'

cd6 :   if (iflag .EQ. 1) then

            Me%Output%BoxFluxes = .true.
            
            inquire(FILE = Me%Files%BoxesFile, EXIST = exist)
cd4 :       if (exist) then
                
                inquire(FILE = Me%Files%BoxesFile, OPENED  = opened)
cd5 :           if (opened) then
                    write(*,*    ) 
                    write(*,'(A)') 'BoxFluxesFileName = ', Me%Files%BoxesFile
                    write(*,*    ) 'Already opened.'
                    stop           'Subroutine StartOutputBoxFluxes; ModuleRunoff. ERR04'    
                end if cd5

                allocate(FluxesOutputList(1), ScalarOutputList(1)) 

                FluxesOutputList = 'runoff_water'
                ScalarOutputList = 'runoff_water'

                call StartBoxDif(BoxDifID         = Me%ObjBoxDif,                    &
                                 TimeID           = Me%ObjTime,                      &
                                 HorizontalGridID = Me%ObjHorizontalGrid,            &
                                 BoxesFilePath    = Me%Files%BoxesFile,              &
                                 FluxesOutputList = FluxesOutputList,                &
                                 ScalarOutputList = ScalarOutputList,                &
                                 WaterPoints2D    = Me%ExtVar%BasinPoints,           &
                                 STAT             = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                         &
                    stop 'Subroutine StartOutputBoxFluxes - ModuleRunoff. ERR15.'

                deallocate(FluxesOutputList, ScalarOutputList) 
                nullify   (FluxesOutputList, ScalarOutputList)
                
            else
                write(*,*) 
                write(*,*)     'Error dont have the file box.'
                write(*,'(A)') 'BoxFileName = ', Me%Files%BoxesFile
                stop           'Subroutine StartOutputBoxFluxes; ModuleRunoff. ERR03'    
            end if cd4
        else
            Me%Output%BoxFluxes = .false.        
        end if cd6
        
    end subroutine StartOutputBoxFluxes

    !--------------------------------------------------------------------------
    
    subroutine ReadConvergenceParameters
    
        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL,               &
                                                       iflag,                   &
                                                       MIN_WATER_COLUMN_STABILIZE_flag    
                                                            
        real                                        :: dummy_real
        
        !----------------------------------------------------------------------    
        
        !----------------------------------------------------------------------
        !Find deprecated keywords in data file
        !----------------------------------------------------------------------
        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, MIN_WATER_COLUMN_STABILIZE_flag,          &
                     SearchType     = FromFile,                                 &
                     keyword        ='MIN_WATER_COLUMN_STABILIZE',              &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR010")

        if (MIN_WATER_COLUMN_STABILIZE_flag > 0) then
            
            write (*,*) '======================================================================='
            write (*,*) 'The following deprecated keywords were found in RunOff data file:'
            write (*,*) ''
            
            if (MIN_WATER_COLUMN_STABILIZE_flag > 0) &
                write(*,*) 'MIN_WATER_COLUMN_STABILIZE: Use STABILIZE_MIN_WATER_COLUMN instead.'
                
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR070")                              
        endif

        !----------------------------------------------------------------------
        !Read convergence options
        !----------------------------------------------------------------------  
        call GetData(Me%CV%Stabilize,                                           &
                     Me%ObjEnterData, iflag,                                    &  
                     keyword      = 'STABILIZE',                                &
                     ClientModule = 'ModuleRunOff',                             &
                     SearchType   = FromFile,                                   &
                     Default      = .false.,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) & 
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR080")        
        if (iflag <= 0) then 
            write(*,*) 'WARNING: Missing STABILIZE keyword in RunOff input data file.'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR081")        
            
        endif
        if (Me%CV%Stabilize) then                
            !Maximun change of water content (in %) allowed in one time step.
            call GetData(Me%CV%StabilizeFactor,                                     &
                         Me%ObjEnterData, iflag,                                    &  
                         keyword      = 'STABILIZE_FACTOR',                         &
                         ClientModule = 'ModuleRunOff',                             &
                         SearchType   = FromFile,                                   &
                         Default      = 0.1,                                        &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) & 
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR082")

            if (Me%CV%StabilizeFactor < 0.0 .or. Me%CV%StabilizeFactor > 1.0) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR083")
                
            call GetData(Me%CV%MinimumValueToStabilize,                     &
                         Me%ObjEnterData, iflag,                            &
                         SearchType   = FromFile,                           &
                         keyword      = 'STABILIZE_MIN_WATER_COLUMN',       &
                         default      = Me%MinimumWaterColumn,              &
                         ClientModule = 'ModuleRunOff',                     &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR084")
            if (Me%CV%MinimumValueToStabilize < Me%MinimumWaterColumn) then
                write (*,*)'Invalid Minimun Water Column to Stabilize value [STABILIZE_MIN_WATER_COLUMN]'
                write (*,*)'Value must be greater than MIN_WATER_COLUMN'            
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR085")
            endif      
            
            call GetData(dummy_real,                                            &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'STABILIZE_RESTART_FACTOR',             &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = 0.,                                     &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR086")
            if (dummy_real <= 0.) then
                Me%CV%MinToRestart = 0
            else
                call CountDomainPoints(dummy_real)
            endif  
            
            call GetData(Me%CV%CheckDecreaseOnly,                               &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'CHECK_DEC_ONLY',                       &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = .false.,                                &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR087")
            
            !Correcting user data can not be the default behaviour
            !User needs to specifically define that wants to correct so default is false
            call GetData(Me%CV%CorrectDischarge,                                &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'STABILIZE_CORRECT_DISCHARGE',          &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = .false.,                                &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR088")   
            
            !Bypass hyraulics correct by default (Paulo suggestion)
            call GetData(Me%CV%CorrectDischargeByPass,                          &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'STABILIZE_CORRECT_DISCHARGE_BYPASS',   &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = .true.,                                 &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR089")             
            
        endif        

       !Number of iterations threshold for starting to ask for a lower DT 
        call GetData(Me%CV%MinIterations,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='MIN_ITERATIONS',                          &
                     Default        = 1,                                        &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR090")
        if (Me%CV%MinIterations < 1) then
            write (*,*)'Invalid Minimun Iterations value [MIN_ITERATIONS]'
            write (*,*)'Value must be greater than 0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR091")
        endif                                 

        !Number of iterations threshold that causes the model to stop
        call GetData(Me%CV%MaxIterations,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='MAX_ITERATIONS',                          &
                     Default        = 1024,                                     &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR100")
        if (Me%CV%MaxIterations < Me%CV%MinIterations) then
            write (*,*)'Invalid Maximun Iterations value [MAX_ITERATIONS]'
            write (*,*)'Value must be greater than the value of MIN_ITERATIONS'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR101")              
        endif
                            
        !% of the maximun iterations that causes the DT to be cut to the value of one internal time step
        call GetData(dummy_real,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'DT_CUT_FACTOR',                    &
                     default      = 0.1,                                &
                     ClientModule = 'ModuleRunOff',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR110") 
        if (dummy_real <= 0.0 .or. dummy_real > 1.0) then
            write (*,*)'Invalid DT Cut Factor [DT_CUT_FACTOR]'
            write (*,*)'Value must be >= 0.0 and < 1.0'        
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR111") 
        endif
        Me%CV%StabilizeHardCutLimit = dummy_real * Me%CV%MaxIterations
        
       !Internal Time Step Split
        call GetData(Me%CV%DTSplitFactor,                                   &
                     Me%ObjEnterData, iflag,                                &
                     keyword      = 'DT_SPLIT_FACTOR',                      &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = 2.0,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadConvergenceParameters - ModuleRunOff - ERR120'        
        if (Me%CV%DTSplitFactor <= 1.0) then
            write (*,*)'Invalid DT Split Factor [DT_SPLIT_FACTOR]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR121")              
        endif            

        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR',                               &
                     Default        = 1.25,                                     &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR130")             
        if (dummy_real <= 1.0) then
            write (*,*)'Invalid DT Factor [DT_FACTOR]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR131")              
        endif            
        
        call GetData(Me%CV%DTFactorUp,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR_UP',                            &
                     Default        = dummy_real,                               &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR140")  
        if (Me%CV%DTFactorUp <= 1.0) then
            write (*,*)'Invalid DT Factor Up [DT_FACTOR_UP]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR141")              
        endif                  
                
        call GetData(Me%CV%DTFactorDown,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR_DOWN',                          &
                     Default        = dummy_real,                               &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR150")  
        if (Me%CV%DTFactorDown <= 1.0) then
            write (*,*)'Invalid DT Factor Down [DT_FACTOR_DOWN]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR151")
        endif                                           
        
        call GetData(Me%CV%LimitDTCourant,                                  &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'LIMIT_DT_COURANT',                     &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &                     
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR180") 
        if (iflag <= 0) then
            write(*,*) 'WARNING: Missing LIMIT_DT_COURANT keyword in RunOff input data file.'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR181")
        endif                    
        if (Me%CV%LimitDTCourant) then
            !Gets Maximum allowed Courant Number
            call GetData(Me%CV%MaxCourant,                                      &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'MAX_COURANT',                          &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = 1.0,                                    &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR181")        
        endif
        
        !----------------------------------------------------------------------
    
    end subroutine ReadConvergenceParameters
    
    !--------------------------------------------------------------------------    

    subroutine ConstructWaterLevelBoundaryConditions
        
        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                      :: CHUNK, i, j, di, dj
        real                                         :: Sum
        integer                                      :: STAT_CALL, n, line
        integer, dimension(:),   pointer             :: VectorI, VectorJ, VectorK
        !Begin-----------------------------------------------------------------

   
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        allocate(Me%WaterLevelBoundaryValue(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%WaterLevelBoundaryValue(:,:) = Me%BoundaryValue

        !$OMP PARALLEL PRIVATE(I,J,di,dj,Sum)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j)  == BasinPoint) then

                !Check if near a boundary point (no diagonal)
do3:            do dj = -1, 1
do4:            do di = -1, 1
                    Sum = dj + di
                    if ((Me%ExtVar%BasinPoints(i+di, j+dj) == 0) .and. (Sum .eq. -1 .or. Sum .eq. 1)) then
                        if(Me%ExtVar%Topography (i, j)  < Me%MaxDtmForBoundary)then
                            Me%BoundaryCells(i,j) = BasinPoint
                            exit do3 
                        endif
                    endif
                enddo do4
                enddo do3
                
            endif
        enddo do2
        enddo do1
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        if(Me%HasBoundaryLines)then

            do line = 1, Me%NumberOfBoundaryLines

                Me%BoundaryLines(line)%nCells = 0
        
                call GetCellZInterceptByLine(Me%ObjHorizontalGrid, Me%BoundaryLines(line)%Line, &
                                             Me%ExtVar%BasinPoints, VectorI, VectorJ, VectorK,  &
                                             Me%BoundaryLines(line)%nCells, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructWaterLevelBoundaryConditions - ModuleRunOff - ERR01' 
        
                if (Me%BoundaryLines(line)%nCells < 1) then
                    write(*,*) 'Boundary line does not intercept any cells'       
                    stop 'ConstructWaterLevelBoundaryConditions - ModuleRunOff - ERR02' 
                endif

                allocate(Me%BoundaryLines(line)%I(1:Me%BoundaryLines(line)%nCells))
                allocate(Me%BoundaryLines(line)%J(1:Me%BoundaryLines(line)%nCells))
        
                do n = 1, Me%BoundaryLines(line)%nCells

                    Me%BoundaryLines(line)%I(n) = VectorI(n)
                    Me%BoundaryLines(line)%J(n) = VectorJ(n)
            
                    Me%BoundaryCells(VectorI(n),VectorJ(n)) = 1

                    if(Me%ExtVar%BasinPoints(VectorI(n),VectorJ(n)) == 0)then
                        stop 'ConstructWaterLevelBoundaryConditions - ModuleRunOff - ERR03' 
                    endif


                    !if boundary line has fixed water level set it here and don't change it again
                    if(.not. Me%BoundaryLines(line)%Variable)then
                        Me%WaterLevelBoundaryValue(VectorI(n), VectorJ(n)) = Me%BoundaryLines(line)%WaterLevelValue
                    endif
            
                enddo

            enddo
        end if
        
        
    end subroutine ConstructWaterLevelBoundaryConditions
    
    !--------------------------------------------------------------------------   

    subroutine ReadBoundaryConditions
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        logical                                     :: BlockFound
        integer                                     :: ClientNumber, nBoundary

        !Begin-----------------------------------------------------------------


        call GetData(Me%MaxDtmForBoundary,                                  &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'MAX_DTM_FOR_BOUNDARY',                 &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadBoundaryConditions - ERR10'        

        if (iflag == 0) then
            write(*,*)'MAX_DTM_FOR_BOUNDARY must be defined in module Runoff'
            stop 'ReadBoundaryConditions - ModuleRunOff - ERR020'
        endif

        call GetData(Me%BoundaryMethod,                                     &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'BOUNDARY_METHOD',                      &
                     Default      = ComputeFlow_,                           &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR30'        

        if (Me%BoundaryMethod /= ComputeFlow_ .and. Me%BoundaryMethod /= InstantaneousFlow_) then
            write(*,*)'BOUNDARY_METHOD must be or 1 - Compute Flow or 2 - Instantaneous FlowOut'
            stop 'ReadBoundaryConditions - ModuleRunOff - ERR040'
        endif

        call GetData(Me%HasBoundaryLines,                                   &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'BOUNDARY_LINES',                       &
                     Default      = .false.,                                &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR50'

        !by default it's false
        Me%BoundaryImposedLevelInTime = .false.

        if(.not. Me%HasBoundaryLines)then !do it the old way, all open boundary cells are treated the same

            !Search for boundary block
            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR60'

            !The boundary block will only contain references to the boundary water level time series  
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_boundary>',                   &
                                        block_end       = '<end_boundary>',                     &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_ .and. BlockFound) then
        
                call ReadLevelTimeSerie(Me%ImposedLevelTS%TimeSerie)

                Me%BoundaryImposedLevelInTime = .true.
            
            else
      
                call GetData(Me%BoundaryValue,                                          &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType = FromFile,                                     &
                             keyword    = 'BOUNDARY_VALUE',                             &
                             ClientModule ='ModuleRunoff',                              &
                             STAT       = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR70' 

                if (iflag == 0) then
                    write(*,*)'if using water elevation boundary, BOUNDARY_VALUE must be defined in ModuleRunoff'
                    stop 'ReadBoundaryConditions - ModuleRunoff - ERR75'
                endif
        
            endif
        
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR076'


        else !use multiple lines to define different boundary conditions at different open boundary cells

            call GetData(Me%BoundaryValue,                                              &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType = FromFile,                                         &
                         keyword    = 'BOUNDARY_VALUE',                                 &
                         ClientModule ='ModuleRunoff',                                  &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR77' 

            if (iflag == 0) then
                write(*,*)'if using water elevation boundary, BOUNDARY_VALUE must be defined in ModuleRunoff'
                stop 'ReadBoundaryConditions - ModuleRunoff - ERR78'
            endif


            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR81'

            Me%NumberOfBoundaryLines = 0

do1:        do         
                call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                            ClientNumber    = ClientNumber,                 &
                                            block_begin     = '<begin_boundary_line>',      &
                                            block_end       = '<end_boundary_line>',        &
                                            BlockFound      = BlockFound,                   &   
                                            STAT            = STAT_CALL)
                if (STAT_CALL == SUCCESS_ .and. BlockFound) then

                    Me%NumberOfBoundaryLines = Me%NumberOfBoundaryLines + 1
                
                else
                
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR90'                             
                
                    exit do1   
                
                endif
            
            enddo do1
            
            !Stop the simluation if no boundary lines have been defined 
            if(Me%NumberOfBoundaryLines > 0)then
    
                allocate(Me%BoundaryLines(1:Me%NumberOfBoundaryLines))

            else

                write(*,*)
                write(*,*)"Option BOUNDARY_LINES is active, but no boundary"
                write(*,*)"lines have been defined in the RunOff input file"
                stop 'ReadBoundaryConditions - ModuleRunoff - ERR91'

            endif

            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR92'

            nBoundary = 0

do2:        do         
                call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                            ClientNumber    = ClientNumber,                     &
                                            block_begin     = '<begin_boundary_line>',          &
                                            block_end       = '<end_boundary_line>',            &
                                            BlockFound      = BlockFound,                       &   
                                            STAT            = STAT_CALL)
                if (STAT_CALL == SUCCESS_ .and. BlockFound) then

                    nBoundary = nBoundary + 1

                    call GetData(Me%BoundaryLines(nBoundary)%FileName,                          &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType = FromBlock,                                        &
                                 keyword    = 'LINE_FILENAME',                                  &
                                 ClientModule ='ModuleRunoff',                                  &
                                 STAT       = STAT_CALL)            
                    if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR95' 

                    if(iflag == 0)then
                        write(*,*)"Open boundary along a line is active "
                        write(*,*)"but line filepath was not defined."
                        write(*,*)"please set keyword LINE_FILENAME"
                        write(*,*)"in RunOff input data file."
                        stop 'ReadBoundaryConditions - ModuleRunoff - ERR100' 
                    endif

                    call New(Me%BoundaryLines(nBoundary)%Line, Me%BoundaryLines(nBoundary)%FileName)

                    call GetData(Me%BoundaryLines(nBoundary)%Variable,                          &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType = FromBlock,                                        &
                                 keyword    = 'VARIABLE_WATER_LEVEL',                           &
                                 ClientModule ='ModuleRunoff',                                  &
                                 STAT       = STAT_CALL)            
                    if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR110' 

                    if(Me%BoundaryLines(nBoundary)%Variable)then

                        !if there's at least one boundary line with variable water level
                        !set BoundaryImposedLevelInTime = .true.  
                        if(.not. Me%BoundaryImposedLevelInTime)then
                            Me%BoundaryImposedLevelInTime = .true.
                        endif

                        call ReadLevelTimeSerie(Me%BoundaryLines(nBoundary)%TimeSerie)

                    else

                        call GetData(Me%BoundaryLines(nBoundary)%WaterLevelValue,               &
                                     Me%ObjEnterData, iflag,                                    &
                                     SearchType = FromBlock,                                    &
                                     keyword    = 'DEFAULTVALUE',                               &
                                     ClientModule ='ModuleRunoff',                              &
                                     STAT       = STAT_CALL)            
                        if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR120' 

                    endif

                
                else
                
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModuleRunoff - ERR130'                             
                
                    exit do2  
                
                endif
            
            enddo do2


        endif
        
    end subroutine ReadBoundaryConditions

    !--------------------------------------------------------------------------

    subroutine ReadLevelTimeSerie(TimeSerie)   
    
        !Arguments-------------------------------------------------------------
        type(T_FromTimeSerieRunOff)                     :: TimeSerie
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(TimeSerie%FileName,                                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'FILENAME',                         &
                     ClientModule = 'FillMatrix',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevelTimeSerie - ModuleRunoff - ERR01'

        if(iflag == 0)then
            write(*,*)"Please set keyword FILENAMEto set the path to the"
            write(*,*)"open boundary time/elevation curve" 
            stop 'ReadLevelTimeSerie - ModuleRunoff - ERR02'
        endif


        call GetData(TimeSerie%DataColumn,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'DATA_COLUMN',                      &
                     ClientModule = 'FillMatrix',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevelTimeSerie - ModuleRunoff - ERR03'

        if(iflag == 0)then
            write(*,*)"Please set keyword DATA_COLUMN to set the time series column"
            write(*,*)"from where to read the open boundary water elevation"
            stop 'ReadLevelTimeSerie - ModuleRunoff - ERR04'
        endif

        call StartTimeSerieInput(TimeSerie%ObjTimeSerie,                &
                                 TimeSerie%FileName,                    &
                                 Me%ObjTime,                            &
                                 CheckDates = .true.,                   &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevelTimeSerie - ModuleRunoff - ERR05'

    end subroutine ReadLevelTimeSerie
    
    !------------------------------------------------------------------

    subroutine ModifyBoundaryLevel

        !Local-----------------------------------------------------------------
        real                                    :: WaterLevelValue
        integer                                 :: nLine, n
        !Begin-----------------------------------------------------------------

        if(Me%HasBoundaryLines)then

            do nLine = 1, Me%NumberOfBoundaryLines
                if(Me%BoundaryLines(nLine)%Variable)then

                    call UpDateLevelValue(Me%BoundaryLines(nLine)%TimeSerie, WaterLevelValue, Me%ExtVar%Now)

                    do n = 1, Me%BoundaryLines(nLine)%nCells
                        Me%WaterLevelBoundaryValue(Me%BoundaryLines(nLine)%I(n), Me%BoundaryLines(nLine)%J(n)) = WaterLevelValue
                    enddo

                endif
            enddo

        else
            
            !boundary values are given by the timeserie value in all boundary grid cells
            call UpDateLevelValue(Me%ImposedLevelTS%TimeSerie, Me%BoundaryValue, Me%ExtVar%Now)

            Me%WaterLevelBoundaryValue(:,:) = WaterLevelValue

        endif

    end subroutine ModifyBoundaryLevel

    !--------------------------------------------------------------------------


    subroutine UpDateLevelValue(TimeSerie, WaterLevelValue, CurrentTime)
        
        !Arguments-------------------------------------------------------------
        type(T_FromTimeSerieRunOff), intent(in)     :: TimeSerie
        real, intent(out)                           :: WaterLevelValue
        type(T_Time), intent(in)                    :: CurrentTime

        !Local-----------------------------------------------------------------
        logical                                     :: TimeCycle
        type (T_Time)                               :: Time1, Time2, InitialDate
        real                                        :: Value1, Value2
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        

        call GetTimeSerieInitialData(TimeSerie%ObjTimeSerie, InitialDate, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UpDateLeveTimeSerielValue - ModuleRunoff - ERR01'

        if (CurrentTime >= InitialDate) then

            !Gets Value for current Time
            call GetTimeSerieValue (TimeSerie%ObjTimeSerie,                             &
                                    CurrentTime,                                        &
                                    TimeSerie%DataColumn,                               &
                                    Time1, Value1, Time2, Value2, TimeCycle,            &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpDateLeveTimeSerielValue - ModuleRunoff - ERR10'

            if (TimeCycle) then
                WaterLevelValue     = Value1
            else
                
                !Interpolates Value for current instant
                call InterpolateValueInTime(CurrentTime, Time1, Value1, Time2, Value2, WaterLevelValue)
                                        
            endif

        else
            write(*,*) 'Water level time series does not have data' 
            write(*,*) 'for the start of the simulation'
            write(*,*) 'Water level time series name: ', Me%ImposedLevelTS%TimeSerie%FileName
            stop 'UpDateLeveTimeSerielValue - ModuleRunoff - ERR20'

        endif

    end subroutine UpDateLevelValue

    !------------------------------------------------------------------
    
    subroutine CountDomainPoints (percent)
    
        !Arguments-------------------------------------------------------------
        real                                        :: percent
        
        !Local----------------------------------------------------------------- 
        integer                                     :: i, j
        integer                                     :: count
        
        !Begin-----------------------------------------------------------------       
                
        count = 0
        
        !Initializes Water Column
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                count = count + 1
            endif

        enddo
        enddo
        
        Me%CV%MinToRestart = max(int(count * percent), 0)
    
    end subroutine CountDomainPoints
    
    !-------------------------------------------------------------------------

    subroutine InitializeVariables

        !Arguments-------------------------------------------------------------
        
        !Local----------------------------------------------------------------- 
        integer                                     :: i, j
        integer                                     :: di, dj
        real                                        :: lowestValue
        
        !Begin-----------------------------------------------------------------       
        
        if (Me%PresentInitialWaterColumn) then
            !Initializes Water Column
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    Me%myWaterLevel(i, j)       = Me%InitialWaterColumn(i,j) + Me%ExtVar%Topography(i, j)
                    Me%MyWaterColumn(i, j)      = Me%InitialWaterColumn(i,j)
                    Me%myWaterVolume(i, j)      = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
                    Me%MyWaterColumnOld(i, j)   = Me%MyWaterColumn(i,j)
                    Me%StabilityPoints(i, j)    = 1
                endif

            enddo
            enddo
        elseif (Me%PresentInitialWaterLevel) then            
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    if (Me%InitialWaterLevel(i,j) > Me%ExtVar%Topography(i, j)) then
                        Me%myWaterLevel(i, j)   = Me%InitialWaterLevel(i,j)
                        Me%MyWaterColumn(i, j)  = Me%InitialWaterLevel(i,j) - Me%ExtVar%Topography(i, j)
                    else
                        Me%myWaterLevel(i, j)   = Me%ExtVar%Topography(i, j)
                        Me%MyWaterColumn(i, j)  = 0.0
                    endif
                    Me%myWaterVolume(i, j)      = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
                    Me%MyWaterColumnOld(i, j)   = Me%MyWaterColumn(i,j)
                    Me%StabilityPoints(i, j)    = 1
                endif

            enddo
            enddo            
        endif
        
        !Finds lowest neighbor for from D8
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
            
                !Finds lowest neighbour
                lowestValue = Me%ExtVar%Topography(i, j)
                do dj = -1, 1
                do di = -1, 1
                    
                    if (dj /= 0 .and. di /= 0 .and. Me%ExtVar%BasinPoints(i+di, j+dj) == BasinPoint) then
                    
                        !Checks lowest neighbor
                        if (Me%ExtVar%Topography(i + di, j + dj) < lowestValue) then
                            
                            lowestValue = Me%ExtVar%Topography(i + di, j + dj)
                            Me%LowestNeighborI(i, j) = i + di
                            Me%LowestNeighborJ(i, j) = j + dj

                        endif
                        
                    endif
                    
                enddo
                enddo        

            endif
        
        enddo
        enddo

        !If drainage network module is associated and simple interaction, then don't apply stability
        !to river points
        if (Me%ObjDrainageNetwork /= 0 .and. Me%SimpleChannelInteraction) then
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%ExtVar%RiverPoints (i, j) == BasinPoint) then
                    Me%StabilityPoints(i, j)    =  0
                endif
                
            enddo
            enddo
        endif
        
        if (Me%RouteDFourPoints) then
            !Checks if a given point is a DFourSink Point -> No point in the four direction is lower then the current point
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                    if (((Me%ExtVar%BasinPoints(i+1, j) == BasinPoint                 .and. &
                          Me%ExtVar%Topography (i+1, j) >= Me%ExtVar%Topography(i, j)) .or.  &  
                          Me%ExtVar%BasinPoints(i+1, j) /= BasinPoint)                .and. &
                        ((Me%ExtVar%BasinPoints(i-1, j) == BasinPoint                 .and. &
                          Me%ExtVar%Topography (i-1, j) >= Me%ExtVar%Topography(i, j)) .or.  &
                          Me%ExtVar%BasinPoints(i-1, j) /= BasinPoint)                .and. &
                        ((Me%ExtVar%BasinPoints(i, j+1) == BasinPoint                 .and. &
                          Me%ExtVar%Topography (i, j+1) >= Me%ExtVar%Topography(i, j)) .or.  &
                          Me%ExtVar%BasinPoints(i, j+1) /= BasinPoint)                .and. &
                        ((Me%ExtVar%BasinPoints(i, j-1) == BasinPoint                 .and. &
                          Me%ExtVar%Topography (i, j-1) >= Me%ExtVar%Topography(i, j)) .or.  &
                          Me%ExtVar%BasinPoints(i, j-1) /= BasinPoint)) then
                        
                        if (Me%LowestNeighborI(i, j) /= i .or. Me%LowestNeighborJ(i, j) /= j) then

                            Me%DFourSinkPoint(i, j) = BasinPoint
                            
                            !D 4 Sink Points are not points where stability criteria is verified
                            Me%StabilityPoints(i, j)= 0
                        
                        endif
                        
                    endif

                endif
            
            enddo
            enddo
        
            !If drainage network modules is associated, then don't apply D4 on drainage network point
            if (Me%ObjDrainageNetwork /= 0) then
            
                if (.not. Me%RouteDFourPointsOnDN) then
                    do j = Me%Size%JLB, Me%Size%JUB
                    do i = Me%Size%ILB, Me%Size%IUB
                    
                        !Source Point is a DNet Point
                        if (Me%ExtVar%RiverPoints (i, j) == BasinPoint) then
                            Me%DFourSinkPoint(i, j) = 0
                        endif
                    
                    enddo
                    enddo
                endif       
            endif
        
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%DFourSinkPoint(i, j) == BasinPoint) then
                
                    if (Me%LowestNeighborI(i, j) /= null_int) then
                
                        !Neighbors of D 4 Sink Points are not points where stability criteria is verified
                        Me%StabilityPoints(Me%LowestNeighborI(i, j), Me%LowestNeighborJ(i, j)) = 0
                
                    endif

                endif

            enddo
            enddo        
        
        endif
                
    end subroutine InitializeVariables

    !--------------------------------------------------------------------------

    subroutine CheckRiverNetWorkConsistency

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real   , dimension(:, :), pointer           :: ChannelsNodeLength 


        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRiverNetWorkConsistency - ModuleRunOff - ERR01'


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
            
                if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) then
                
                    if (ChannelsNodeLength(i, j) < 0.0) then
                        write(*,*)'Inconsistent River Network', i, j
                        stop 'CheckRiverNetWorkConsistency - ModuleRunOff - ERR02'
                    endif
                
                else
                
                    if (ChannelsNodeLength(i, j) > 0.0) then
                        write(*,*)'Inconsistent River Network', i, j
                        stop 'CheckRiverNetWorkConsistency - ModuleRunOff - ERR03'
                    endif
                
                endif

            endif

        enddo
        enddo
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRiverNetWorkConsistency - ModuleRunOff - ERR04'
    
    
    end subroutine CheckRiverNetWorkConsistency

    !--------------------------------------------------------------------------
    
    subroutine ConstructDischarges

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        character(len=StringLength)                 :: DischargeName
        real                                        :: CoordinateX, CoordinateY
        logical                                     :: CoordinatesON, IgnoreOK
        integer                                     :: Id, Jd, dn, DischargesNumber
        integer                                     :: STAT_CALL
        type (T_Lines),   pointer                   :: LineX
        type (T_Polygon), pointer                   :: PolygonX
        integer, dimension(:),   pointer            :: VectorI, VectorJ, VectorK
        integer                                     :: SpatialEmission, nCells       

        call Construct_Discharges(Me%ObjDischarges, Me%ObjTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR01' 
                
        call GetDischargesNumber(Me%ObjDischarges, DischargesNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR02' 

        do dn = 1, DischargesNumber

            call GetDischargesGridLocalization(Me%ObjDischarges, dn,            &
                                               CoordinateX   = CoordinateX,     &
                                               CoordinateY   = CoordinateY,     & 
                                               CoordinatesON = CoordinatesON,   &
                                               STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR03' 
                    
            call GetDischargesIDName (Me%ObjDischarges, dn, DischargeName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR03' 

            if (CoordinatesON) then
                
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordinateX, CoordinateY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR04' 

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreDischarge(Me%ObjDischarges, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR05' 

                    if (IgnoreOK) then
                        write(*,*) 'Discharge outside the domain - ',trim(DischargeName),' - ',trim(Me%ModelName)
                        cycle
                    else
                        stop 'ModuleRunOff - ConstructDischarges - ERR06' 
                    endif

                endif

                call CorrectsCellsDischarges(Me%ObjDischarges, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR07' 
                    
            endif

            !ATTENTION - NEED TO VERIFY IF DISCHARGES ARE COLLINEAR.
            !Do not allow with properties since the flow used in PMP is not distributed by discharges
            !and will be accounted with flow duplicating
            call GetDischargeSpatialEmission(Me%ObjDischarges, dn, LineX, PolygonX, &
                                             SpatialEmission, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR08' 
                    
            if (SpatialEmission == DischPoint_) then
 
                call GetDischargesGridLocalization(Me%ObjDischarges, dn,            &
                                                   Igrid         = Id,              &
                                                   JGrid         = Jd,              &
                                                   STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR09' 

                if (Me%ExtVar%BasinPoints(Id,Jd) /= WaterPoint) then
                    call TryIgnoreDischarge(Me%ObjDischarges, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR10' 

                    write(*,*) 'Discharge outside the domain I=',Id,' J=',Jd,'Model name=',trim(Me%ModelName)

                    if (IgnoreOK) then
                        write(*,*) 'Discharge in a land cell - ',trim(DischargeName),' - ',trim(Me%ModelName)
                        cycle
                    else
                        stop 'ModuleRunOff - ConstructDischarges - ERR11' 
                    endif
                endif

                nCells    = 1
                allocate(VectorI(nCells), VectorJ(nCells))
                VectorJ(nCells) = Jd
                VectorI(nCells) = Id

            else

                if (SpatialEmission == DischLine_) then
                    call GetCellZInterceptByLine(Me%ObjHorizontalGrid, LineX,       &
                                                 Me%ExtVar%BasinPoints, VectorI, VectorJ, VectorK,   &
                                                 nCells, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR12' 

                    if (nCells < 1) then
                        write(*,*) 'Discharge line intercept 0 cells'       
                        stop 'ModuleRunOff - ConstructDischarges - ERR13' 
                    endif

                endif 


                if (SpatialEmission == DischPolygon_) then
                    call GetCellZInterceptByPolygon(Me%ObjHorizontalGrid, PolygonX, &
                                                 Me%ExtVar%BasinPoints, VectorI, VectorJ, VectorK,   &
                                                 nCells, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR14' 

                    if (nCells < 1) then
                        write(*,*) 'Discharge contains 0 center cells'       
                        write(*,*) 'Or the polygon is to small and is best to a discharge in a point or'
                        write(*,*) 'the polygon not define properly'
                        stop 'ModuleRunOff - ConstructDischarges - ERR15' 
                    endif

                endif


            endif
                        
            
            call SetLocationCellsZ (Me%ObjDischarges, dn, nCells, VectorI, VectorJ, VectorK, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR16' 
            

        enddo
        
        if (Me%OutPut%TimeSerieDischON) then
            call Construct_Time_Serie_Discharge
        endif            

   
    end subroutine ConstructDischarges
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine Construct_Time_Serie_Discharge

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: dis, i, j
        character(len=StringLength)                         :: Extension, DischargeName

        !Begin-----------------------------------------------------------------


        call GetDischargesNumber(Me%ObjDischarges, Me%OutPut%DischargesNumber, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_)stop 'Construct_Time_Serie_Discharge - ModuleRunOff - ERR10'

        allocate(Me%OutPut%TimeSerieDischID(Me%OutPut%DischargesNumber))
        
        Me%OutPut%TimeSerieDischID(:) = 0
        
        Me%OutPut%TS_Numb_DischProp = 6 !1 - flow; 2 - velocity; 3 - Area, 4 - water level Upstream ; 
                        !5 - water level Downstream; 6 - water flow without corrections
        
        allocate(Me%OutPut%TimeSerieDischProp(1:Me%OutPut%DischargesNumber,1:Me%OutPut%TS_Numb_DischProp))
        
        Me%OutPut%TimeSerieDischProp(:,:) = 0.

        !Allocates PropertyList
        allocate(PropertyList(Me%OutPut%TS_Numb_DischProp), STAT = STAT_CALL)

        if (STAT_CALL/=SUCCESS_)stop 'Construct_Time_Serie_Discharge - ModuleRunOff - ERR20'

        !Fills up PropertyList
        PropertyList(1) = "water_flux"
        PropertyList(2) = "velocity"
        PropertyList(3) = "area"        
        PropertyList(4) = "water_level_upstream"
        PropertyList(5) = "water_level_downstream"
        PropertyList(6) = "water_flow_no_correction"        

        do i=1,Me%OutPut%TS_Numb_DischProp
            do j=1,len_trim(PropertyList(i))
                if (PropertyList(i)(j:j)==' ') PropertyList(i)(j:j)='_'
            enddo
        enddo
  
        Extension = 'srd'

        do dis = 1, Me%OutPut%DischargesNumber
        
            call GetDischargesIDName (Me%ObjDischarges, dis, DischargeName, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_)stop 'Construct_Time_Serie_Discharge - ModuleRunOff - ERR60'
        
            call StartTimeSerie(TimeSerieID         = Me%OutPut%TimeSerieDischID(dis),      &
                                ObjTime             = Me%ObjTime,                           &
                                TimeSerieDataFile   = Me%Output%DiscTimeSerieLocationFile,  &
                                PropertyList        = PropertyList,                         &
                                Extension           = Extension,                            &
                                ResultFileName      = "hydro_"//trim(DischargeName),        &
                                STAT                = STAT_CALL)
            if (STAT_CALL/=SUCCESS_)stop 'Construct_Time_Serie_Discharge - ModuleRunOff - ERR70'
            
        enddo
        
        !----------------------------------------------------------------------
        
        
    end subroutine Construct_Time_Serie_Discharge

    !--------------------------------------------------------------------------    
    
    subroutine CheckHorizontalGridRotation
    
        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetCheckDistortion (Me%ObjHorizontalGrid,                                  &
                                 Me%ExtVar%Distortion,                                  &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckHorizontalGridRotation - ModuleRunOff - ERR01'
            

        if (Me%ExtVar%Distortion) then

            call GetGridRotation(Me%ObjHorizontalGrid,                                  &
                                 Me%ExtVar%RotationX,                                   &
                                 Me%ExtVar%RotationY,                                   &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckHorizontalGridRotation - ModuleRunOff - ERR02'





        else

            call GetGridAngle(Me%ObjHorizontalGrid,                                     &
                              Me%ExtVar%GridRotation,                                   &
                              STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'CheckHorizontalGridRotation - ModuleRunOff - ERR03'

            !Convert from degrees in to radians
            Me%ExtVar%GridRotation = Me%ExtVar%GridRotation * Pi / 180.

            Me%GridCosAngleX = cos(Me%ExtVar%GridRotation)
            Me%GridCosAngleY = cos(Pi/2. + Me%ExtVar%GridRotation)

            Me%GridSinAngleX = sin(Me%ExtVar%GridRotation)
            Me%GridSinAngleY = sin(Pi/2. + Me%ExtVar%GridRotation)

        endif
    
    
    end subroutine CheckHorizontalGridRotation

    !--------------------------------------------------------------------------    

    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------
        !Local----------------------------------------------------------------- 
        !Begin-----------------------------------------------------------------       

        allocate(Me%iFlowToChannels  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%lFlowToChannels  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%iFlowToChannels      = 0.0   !Sets values initially to zero, so 
        Me%lFlowToChannels      = 0.0   !model can run without DNet
        
        allocate(Me%lFlowBoundary    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%iFlowBoundary    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%lFlowBoundary        = 0.0   !Sets values initially to zero, so 
        Me%iFlowBoundary        = 0.0   !model can run without BC

        allocate(Me%lFlowDischarge    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%iFlowDischarge    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%lFlowDischarge        = 0.0   !Sets values initially to zero, so 
        Me%iFlowDischarge        = 0.0   !model can run without Dis

        allocate(Me%iFlowRouteDFour  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%iFlowRouteDFour       = 0.0   !Sets values initially to zero, so  
        
        allocate(Me%BoundaryCells     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%BoundaryCells = 0
        
        allocate(Me%myWaterLevel         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterColumn        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        
        
        allocate(Me%myWaterColumnAfterTransport (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterVolumePred   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%myWaterVolumePred = null_real
        allocate(Me%InitialWaterColumn   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%InitialWaterLevel    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterVolume        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterColumnOld     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterVolumeOld     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%MassError            (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%myWaterLevel            = null_real
        Me%myWaterColumn           = null_real
        Me%InitialWaterColumn      = null_real
        Me%InitialWaterLevel       = null_real
        Me%myWaterVolume           = 0.0        !For OpenMI
        Me%myWaterColumnOld        = null_real
        Me%myWaterVolumeOld        = null_real
        Me%MassError               = 0
        Me%myWaterColumnAfterTransport = null_real

        allocate(Me%iFlowX               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%iFlowY               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%lFlowX               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%lFlowY               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%FlowXOld             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%FlowYOld             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%InitialFlowX         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%InitialFlowY         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AreaU                (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AreaV                (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ComputeFaceU         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ComputeFaceV         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%OpenPoints           (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%NoAdvectionPoints    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ComputeAdvectionU    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ComputeAdvectionV    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%VelModFaceU          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%VelModFaceV          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))



        Me%iFlowX               = 0.0
        Me%iFlowY               = 0.0
        Me%lFlowX               = 0.0
        Me%lFlowY               = 0.0
        Me%FlowXOld             = 0.0
        Me%FlowYOld             = 0.0
        Me%InitialFlowX         = 0.0
        Me%InitialFlowY         = 0.0
        Me%AreaU                = AlmostZero
        Me%AreaV                = AlmostZero
        Me%ComputeFaceU         = 0
        Me%ComputeFaceV         = 0
        Me%OpenPoints           = 0
        
        Me%NoAdvectionPoints    = 0.0
        Me%ComputeAdvectionU    = 1
        Me%ComputeAdvectionV    = 1

        Me%VelModFaceU          = 0.0
        Me%VelModFaceV          = 0.0

        allocate(Me%OverLandCoefficient  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%OverLandCoefficientDelta (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%OverLandCoefficientX (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%OverLandCoefficientY (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%OverLandCoefficient       = null_real
        Me%OverLandCoefficientDelta  = null_real
        Me%OverLandCoefficientX      = null_real
        Me%OverLandCoefficientY      = null_real


       
        allocate (Me%CenterFlowX    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%CenterFlowY    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%FlowModulus    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%CenterVelocityX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%CenterVelocityY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%VelocityModulus(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        allocate (Me%Output%MaxFlowModulus (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%Output%MaxWaterColumn (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB)) 
        allocate (Me%Output%TimeOfMaxWaterColumn (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB)) 
        
        Me%Output%MaxFlowModulus = null_real
        Me%Output%MaxWaterColumn = Me%MinimumWaterColumn
        Me%Output%TimeOfMaxWaterColumn = -99.0
 
        allocate (Me%Output%VelocityAtMaxWaterColumn (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB)) 
        Me%Output%VelocityAtMaxWaterColumn = null_real

        allocate (Me%Output%MaxFloodRisk (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB)) 
        Me%Output%MaxFloodRisk = null_real
      
        allocate (Me%Output%FloodPeriod (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB)) 
        Me%Output%FloodPeriod = 0.

        allocate (Me%Output%FloodArrivalTime (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB)) 
        Me%Output%FloodArrivalTime = -99.0


        allocate (Me%LowestNeighborI (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%LowestNeighborJ (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%DFourSinkPoint  (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%StabilityPoints (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                    
        Me%LowestNeighborI = null_int
        Me%LowestNeighborJ = null_int
        Me%DFourSinkPoint  = 0
        Me%StabilityPoints = 0
        
        if (Me%ObjDrainageNetwork /= 0) then
            !allocate(Me%BoundaryRiverLevel   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            !Me%BoundaryRiverLevel      = null_real
            allocate(Me%NodeRiverLevel   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%NodeRiverLevel      = null_real  
            
            allocate(Me%MarginRiverLevel   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%MarginRiverLevel      = null_real        
            
            allocate(Me%MarginFlowToChannels   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%MarginFlowToChannels      = null_real              
        endif


    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine ConstructOverLandCoefficient

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB    
        integer                                             :: i, j

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !TODO: OpenMP - Missing implementation
        do j = JLB, JUB + 1
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j) + Me%ExtVar%BasinPoints(i, j-1) == 2) then !Two Basin Points
            
                Me%OverlandCoefficientX(i, j) = (Me%ExtVar%DUX(i, j  ) * Me%OverlandCoefficient(i, j-1  )  + &
                                                 Me%ExtVar%DUX(i, j-1) * Me%OverlandCoefficient(i, j)) / &
                                                 (Me%ExtVar%DUX(i, j) + Me%ExtVar%DUX(i, j-1))
            endif

        enddo
        enddo

        do j = JLB, JUB
        do i = ILB, IUB + 1

            if (Me%ExtVar%BasinPoints(i, j) + Me%ExtVar%BasinPoints(i-1, j) == 2) then !Two Basin Points
            
                Me%OverlandCoefficientY(i, j) =     (Me%ExtVar%DVY(i, j  ) * Me%OverlandCoefficient(i-1, j  )  + &
                                                     Me%ExtVar%DVY(i-1, j) * Me%OverlandCoefficient(i, j)) / &
                                                     (Me%ExtVar%DVY(i, j) + Me%ExtVar%DVY(i-1, j))
            endif

        enddo
        enddo


    end subroutine ConstructOverLandCoefficient
    
        !--------------------------------------------------------------------------

    subroutine ConstructAdvectionZones

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB    
        integer                                             :: i, j, STAT_CALL, ClientNumber
        logical                                             :: BlockFound
        
        if(Me%NoAdvectionZones)then
            
            call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructAdvectionZones - ModuleRunOff - ERR01'

            !Gets Block for OverLand Coef
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,           &
                                        '<BeginNoAdvectionZones>',               &
                                        '<EndNoAdvectionZones>', BlockFound,     &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructAdvectionZones - ModuleRunOff - ERR02'
            
            if (BlockFound) then
                call ConstructFillMatrix  ( PropertyID       = Me%NoAdvectionZonesID,        &
                                            EnterDataID      = Me%ObjEnterData,              &
                                            TimeID           = Me%ObjTime,                   &
                                            HorizontalGridID = Me%ObjHorizontalGrid,         &
                                            ExtractType      = FromBlock,                    &
                                            PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                            Matrix2D         = Me%NoAdvectionPoints,         &
                                            TypeZUV          = TypeZ_,                       &
                                            STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructAdvectionZones - ModuleRunOff - ERR03'

                call KillFillMatrix(Me%NoAdvectionZonesID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    write(*,*)"STAT_CALL = ", STAT_CALL
                    stop 'ConstructAdvectionZones - ModuleRunOff - ERR04'
                endif
                
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructAdvectionZones - ModuleRunOff - ERR04'

                
            endif
            
            !Check that advection zones are 0 or 1. 
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    if (int(Me%NoAdvectionPoints(i,j)) .lt. 0) then
                        write(*,*) 'No Advection Zones can only be 0 or 1'
                        write(*,*) 'in cell', i, j, Me%NoAdvectionPoints(i,j), int(Me%NoAdvectionPoints(i,j))
                        stop 'ConstructAdvectionZones - ModuleRunoff - ERR05'
                    endif
                    
                    if (int(Me%NoAdvectionPoints(i,j)) .gt. 1) then
                        write(*,*) 'No Advection Zones can only be 0 or 1'
                        write(*,*) 'in cell', i, j, Me%NoAdvectionPoints(i,j), int(Me%NoAdvectionPoints(i,j))
                        stop 'ConstructAdvectionZones - ModuleRunoff - ERR06'
                    endif
                
                endif
                
            enddo
            enddo

        else
            write(*,*)'Missing Block <BeginAdvectionZones> / <BeginAdvectionZones>' 
            stop      'ConstructAdvectionZones - ModuleRunoff - ERR07'
        endif
        

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !TODO: OpenMP - Missing implementation
        do j = JLB, JUB + 1
        do i = ILB, IUB
            if (Me%ExtVar%BasinPoints(i, j) + Me%ExtVar%BasinPoints(i, j-1) == 2) then !Two Basin Points
                if (int(Me%NoAdvectionPoints(i, j) + Me%NoAdvectionPoints(i, j-1)) == 2) then 
                    Me%ComputeAdvectionU(i, j) = 0
                endif
            endif
        enddo
        enddo

        do j = JLB, JUB
        do i = ILB, IUB + 1
            if (Me%ExtVar%BasinPoints(i, j) + Me%ExtVar%BasinPoints(i-1, j) == 2) then !Two Basin Points
                if ((Me%NoAdvectionPoints(i, j) + Me%NoAdvectionPoints(i-1, j)) == 2) then 
                    Me%ComputeAdvectionV(i, j) = 0
                endif
            endif
        enddo
        enddo

        deallocate(Me%NoAdvectionPoints)
        nullify(Me%NoAdvectionPoints)

    end subroutine ConstructAdvectionZones

    !--------------------------------------------------------------------------

    subroutine ConstructSewerStormWater
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
#ifdef _SEWERGEMSENGINECOUPLER_

        integer                                             :: ILB, IUB, JLB, JUB    

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
                
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        Me%TotalStormWaterVolume    = 0.0
        Me%TotalBoundaryFlowVolume  = 0.0
        Me%TotalBoundaryInflowVolume= 0.0
        Me%TotalBoundaryOutflowVolume = 0.0
        Me%TotalInletsVolume        = 0.0
        Me%TotalManholesVolume      = 0.0
        Me%TotalHeadwallsVolume     = 0.0
        Me%TotalPondsVolume         = 0.0
        Me%TotalOpenChannelVolume   = 0.0
        Me%TotalOutfallsVolume      = 0.0
        Me%Total1DVolume            = 0.0
        Me%Total1D2DVolume          = 0.0
        Me%MaxTotal1DVolume         = 0.0
        Me%MaxTotal2DVolume         = 0.0
        Me%MaxTotal1D2DVolume       = 0.0
        Me%TimeOfMaxTotal1DVolume   = 0.0
        Me%TimeOfMaxTotal2DVolume   = 0.0
        Me%TimeOfMaxTotal1D2DVolume = 0.0

        if (Me%StormWaterModel .and. Me%ObjDrainageNetwork /= 0) then
            write(*,*)'It is not possible to activate 1D Drainage Network and SWMM at the same time'
            stop 'ConstructSewerStormWater - ModuleRunOff - ERR01'            
        endif

        call ConstructSewerGEMS 

#else _SEWERGEMSENGINECOUPLER_    
    print*, 'Executable not compiled to use SewerGEMS Engine coupler, aborting'
    stop
#endif _SEWERGEMSENGINECOUPLER_
            
            
    end subroutine ConstructSewerStormWater
    
    !--------------------------------------------------------------------------
#ifdef _SEWERGEMSENGINECOUPLER_
    
    subroutine ConstructSewerGEMS
        
        !--------------------------------------------------------------------------
        integer                                         :: STAT_CALL, n, m, i, c, dpos, SaveResults, nodeType
        integer                                         :: ObjStormWaterEnterData = 0, iflag, xn
        character(len = :, kind = c_char), allocatable  :: inpFile, rptFile, outFile
        character(len = 99, kind = c_char)              :: nodeName
        character(len = 99, kind = c_char)              :: startDate, startTime, endDate, endTime
        character(len = :), allocatable                 :: year, month, day, hour, minute, second, fullDate
        type(T_Time)                                    :: SWMMBeginTime, SWMMEndTime
        integer(c_int)                                  :: isOpenChannel
        integer                                         :: nInlets, iNode, nManholes, nOutfalls, nCrossSections, nPonds
        logical                                         :: Exists
        integer                                         :: UncoupledElementsFileID, nHeadwalls
        !--------------------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Coupling SewerGEMS..."
        write(*,*)

        call ReadFileName('STORMWATER_DAT', Me%Files%SWMMdat, Message = "SWMM MOHID dat file", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR01'
        
        call ReadFileName('STORMWATER_HDF', Me%Files%SWMMHDF, Message = "SWMM HDF file", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR10'

        call ReadFileName('ROOT_SRT', Me%Files%SWMMTimeSeriesDir, Message = "SWMM Time Series folder", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR20'

        Me%Files%SWMMUncoupledElements = trim(adjustl(Me%Files%SWMMTimeSeriesDir))//"UncoupledElements.dat"

        call UnitsManager(UncoupledElementsFileID, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR21'

        open(UNIT = UncoupledElementsFileID, FILE = trim(Me%Files%SWMMUncoupledElements), STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR22'


        dpos = scan(trim(Me%Files%SWMMHDF),".", BACK= .true.)
    
        !name of output files is the same as the output hdf5
        if (dpos > 0) then
            Me%Files%SWMMrpt = Me%Files%SWMMHDF(1:dpos)//'rpt'
            Me%Files%SWMMout = Me%Files%SWMMHDF(1:dpos)//'swm'
        else
            write(*,*)'STORMWATER_HDF file extension in nomfich.dat not recognized'
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR30' 
        end if
        
        call ConstructEnterData (ObjStormWaterEnterData, Me%Files%SWMMdat, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR40' 

        call GetData(Me%Files%SWMMinp,                                                   &
                     ObjStormWaterEnterData, iflag,                                      &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MODEL_CONFIGURATION',                               &
                     ClientModule = 'ModuleStormWater',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR50' 

        if(iflag == 0)then
            write(*,*)"Please set keyword MODEL_CONFIGURATION."
            write(*,*)"File = " , trim(adjustl(Me%Files%SWMMdat))
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR60'
        endif

        inquire (FILE=trim(Me%Files%SWMMinp), EXIST = Exists)
        if(.not. Exists)then
            write(*,*)"Could not find SWMM stormwater file"
            write(*,*)"File = " , trim(adjustl(Me%Files%SWMMinp))
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR070'
        end if

        
        call GetData(Me%Files%SWMMTimeSeries,                                            &
                     ObjStormWaterEnterData, iflag,                                      &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TIME_SERIE_LOCATION',                               &
                     ClientModule = 'ModuleStormWater',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR80' 

        if(iflag == 0)then
            write(*,*)
            write(*,*)"Could not find path to SWMM time series config file"
            write(*,*)"using keyword TIME_SERIE_LOCATION in"
            write(*,*)"File = " , trim(adjustl(Me%Files%SWMMdat))
            write(*,*)'ConstructSewerGEMS - ModuleRunOff - WRNG08a'
            write(*,*)"SWMM timeseries results will not be converted to MOHID format" 
            Me%Files%SWMMTimeSeries     = ""
            Me%Files%SWMMTimeSeriesDir  = ""
        endif 

        
        call KillEnterData (ObjStormWaterEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR90' 

        
        inpFile = trim(ADJUSTL(Me%Files%SWMMinp))//C_NULL_CHAR
        rptFile = trim(ADJUSTL(Me%Files%SWMMrpt))//C_NULL_CHAR
        outFile = trim(ADJUSTL(Me%Files%SWMMout))//C_NULL_CHAR
        saveResults = 1

        STAT_CALL = SewerGEMSEngine_open(inpFile, rptFile, outFile)
        if (STAT_CALL /= SUCCESS_) then
       
            write(*,*)
            write(*,*)"Error initializing SWMM"
            write(*,*)"SWMM error code : ", STAT_CALL
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR100'
        endif
        
        STAT_CALL = SewerGEMSEngine_start(SaveResults)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR110'

        STAT_CALL = SewerGEMSEngine_getDates(startDate, startTime, endDate, endTime)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR120'
        
        month   = startDate(1:2)
        day     = startDate(4:5)
        year    = startDate(7:10)
        hour    = startTime(1:2)
        minute  = startTime(4:5)
        second  = startTime(7:8)
        fullDate = year//" "//month//" "//day//" "//hour//" "//minute//" "//second
        SWMMBeginTime = ConvertStringtToDate(fullDate, .true.)
        !extracting final date components from string, quick and dirty
        month   = endDate(1:2)
        day     = endDate(4:5)
        year    = endDate(7:10)
        hour    = endTime(1:2)
        minute  = endTime(4:5)
        second  = endTime(7:8)
        fullDate= year//" "//month//" "//day//" "//hour//" "//minute//" "//second
        SWMMEndTime = ConvertStringtToDate(fullDate, .true.)
    
        write(*,*)
        if (SWMMBeginTime .ne. Me%BeginTime .or. SWMMEndTime .ne. Me%EndTime) then
            print*, 'SewerGEMSEngine dates are not consistent with MOHID'
            print*, 'Please edit, run the coupling tools and try again'
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR125'
        else
            print*, 'Start end dates consistent, continuing'
        end if
        
        STAT_CALL = SewerGEMSEngine_getNumberOfNodes(Me%NumberOfNodes)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR130'
        
        write(*,*)
        write(*,*)"SewerGEMS - Total Number of Nodes :", Me%NumberOfNodes
        
        allocate(Me%NodesX      (1:Me%NumberOfNodes))
        allocate(Me%NodesY      (1:Me%NumberOfNodes))
        allocate(Me%NodesI      (1:Me%NumberOfNodes))
        allocate(Me%NodesJ      (1:Me%NumberOfNodes))
        allocate(Me%NodesID     (1:Me%NumberOfNodes))
        allocate(Me%NodesCellID (1:Me%NumberOfNodes))
        allocate(Me%NodesNames  (1:Me%NumberOfNodes))
        
        !0 = Inactive           -> Outside domain or not connected to 2D RunOff e.g. bolted manhole
        !1 = Manhole            -> can flow from SewerGEMS to 2D RunOff; cannot flow from 2D to SewerGEMS
        !2 = Inlet/Catch Basin  -> can flow from SewerGEMS to 2D RunOff and from 2D to SewerGEMS
        !3 = Cross Section      -> can flow from SewerGEMS to 2D RunOff and from 2D to SewerGEMS
        !4 = Outfall            -> can flow from SewerGEMS to 2D RunOff and from 2D to SewerGEMS
        allocate(Me%NodesType   (1:Me%NumberOfNodes))
        
        Me%ActiveNodes          = 0 
        Me%InactiveNodes        = 0 
        Me%NumberOfManholes     = 0 
        Me%NumberOfCrossSections= 0 
        Me%NumberOfOutfalls     = 0
        nInlets                 = 0 
        nPonds                  = 0 
        nHeadwalls              = 0 

        do n = 1, Me%NumberOfNodes
            
            !use (n-1) as SWMM node indexes start with zero
            Me%NodesID(n) = n-1
            
            nodeName = " "//C_NULL_CHAR
            STAT_CALL = SewerGEMSEngine_getNodeName(Me%NodesID(n), nodeName)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR140'
            
            do c = 1, 99
                i = iachar(nodeName(c:c))
                if(i == 0)then
                    nodeName(c:c) = space
                endif
            enddo

            !do c = 1, 99
               !if (nodeName(c:c) == C_NULL_CHAR) nodeName(c:c) = space
            !end do

            Me%NodesNames(n) = nodeName

            if(IgnoreSWMMNode(Me%NodesNames(n))) then

                Me%NodesType(n)     = NotCoupled_

                Me%InactiveNodes = Me%InactiveNodes + 1 !add to list of inactive nodes

                write(UncoupledElementsFileID, *) Me%NodesID(n)

                call SetError(WARNING_, INTERNAL_, Me%NodesNames(n), OFF)

                cycle !go directly to next node
            endif 

            !Check the type of node (0 = junction; 1 = outfall; 2 = storage)
            STAT_CALL = SewerGEMSEngine_getNodeType(Me%NodesID(n), nodeType)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR150'
            
            !Ponds - check if pond in SWMM file is valid and listed in the Ponds file; 
            !2D engine does not check for (I,J) location of pond as with other nodes
            if(nodeType == SWMMStorage_)then
                if(ValidatePond(Me%NodesNames(n), Me%NodesID(n)) > 0)then !if not pond returns null_int (-999999)
                    Me%NodesType(n) = Pond_
                    nPonds          = nPonds + 1
                    cycle !go directly to next node
                else
                    write(*,*)"Node is defined as pond/storage area in SWMM inp file"
                    write(*,*)"But it's not listed in the ponds file"
                    write(*,*)"Node name = ", trim(adjustl(Me%NodesNames(n)))
                    stop 'ConstructSewerGEMS - ModuleRunOff - ERR210'
                endif

            endif

            !Get node coordinates
            STAT_CALL = SewerGEMSEngine_getNodeXY(Me%NodesID(n), Me%NodesX(n), Me%NodesY(n))
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR160'
            
            !Check if node is located inside the computational grid
insidegrid: if (GetXYInsideDomain(Me%ObjHorizontalGrid, Me%NodesX(n), Me%NodesY(n))) then
                
                !Get the nodes I and J grid cell indexes from node's X and Y
                call GetXYCellZ(Me%ObjHorizontalGrid, Me%NodesX(n), Me%NodesY(n),         &
                                                      Me%NodesI(n), Me%NodesJ(n))
                

                !Check if grid cell is an active compute point
ifactivepoint:  if(Me%ExtVar%BasinPoints(Me%NodesI(n), Me%NodesJ(n)) == 1) then
                
                    if(nodeType == SWMMJunction_)then
                        
                        !Check the node is an open channel 
                        STAT_CALL = SewerGEMSEngine_getIsNodeOpenChannel(Me%NodesID(n), IsOpenChannel)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR170'

                        if(IsOpenChannel == 1)then
                            write(*,*)
                            write(*,*)"Cross section nodes cannot be located in active grid cell"
                            write(*,*)"Node name = ", trim(adjustl(Me%NodesNames(n)))
                            stop 'ConstructSewerGEMS - ModuleRunOff - ERR180' 
                        else
                            !Check inlet name and set inlet node index
                            if(ValidateInlet(Me%NodesNames(n), Me%NodesID(n)) > 0)then !if not inlet returns null_int (-999999)
                                Me%NodesType(n) = Inlet_
                                nInlets         = nInlets + 1
                            elseif(ValidateHeadwall(Me%NodesNames(n), Me%NodesID(n)) > 0)then
                                Me%NodesType(n) = Headwall_
                                nHeadwalls      = nHeadwalls + 1
                            else
                                !if it's not an inlet then it's a manhole
                                Me%NodesType(n)     = Manhole_ 
                                Me%NumberOfManholes = Me%NumberOfManholes + 1
                            endif
                        endif
                   
                    elseif(nodeType == SWMMOutfall_)then

                        Me%NodesType(n)     = Outfall_ 
                        Me%NumberOfOutfalls = Me%NumberOfOutfalls + 1

                    elseif(nodeType == SWMMStorage_)then

                        write(*,*)
                        write(*,*)"Ponds/Storage nodes cannot be located in active grid cell"
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR190' 

                    endif
                    !Convert I and J to cell index
                    call GetCellIDfromIJ(Me%ObjHorizontalGrid, Me%NodesI(n), Me%NodesJ(n), Me%NodesCellID(n))

                else ifactivepoint

                    !Check the node is an open channelF
                    !Cross section nodes are located in non-compute points
                    STAT_CALL = SewerGEMSEngine_getIsNodeOpenChannel(Me%NodesID(n), IsOpenChannel)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR200'

                    !0 = junction
                    !1 = outfall
                    !2 = storage
                    if(nodeType == SWMMJunction_)then

                        if(IsOpenChannel == 1)then
                            Me%NodesType(n)          = CrossSection_
                            Me%NumberOfCrossSections = Me%NumberOfCrossSections + 1 
                        else
                            !if node is not a cross section and is in a non-compute grid cell then it's ignored. 
                            Me%NodesType(n)     = NotCoupled_
                            Me%NodesCellID(n)   = null_int
                            Me%InactiveNodes    = Me%InactiveNodes + 1

                            write(UncoupledElementsFileID, *) Me%NodesID(n)
                            call SetError(WARNING_, INTERNAL_, Me%NodesNames(n), OFF)

                        endif

                    else
                        !if node is not a cross section nor pond and is in a non-compute grid cell then it's ignored. 
                        Me%NodesType(n)     = NotCoupled_
                        Me%NodesCellID(n)   = null_int
                        Me%InactiveNodes    = Me%InactiveNodes + 1

                        write(UncoupledElementsFileID, *) Me%NodesID(n)
                        call SetError(WARNING_, INTERNAL_, Me%NodesNames(n), OFF)

                    endif

                endif ifactivepoint
            
            else insidegrid
                !if node is in outside the computational grid then it's ignored. 
                Me%NodesI(n)        = null_int
                Me%NodesJ(n)        = null_int
                Me%NodesType(n)     = NotCoupled_
                Me%NodesCellID(n)   = null_int
                Me%InactiveNodes    = Me%InactiveNodes + 1

                write(UncoupledElementsFileID, *) Me%NodesID(n)
                call SetError(WARNING_, INTERNAL_, Me%NodesNames(n), OFF)

            endif insidegrid

        enddo

        call UnitsManager(UncoupledElementsFileID, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR211'


        !Double check inlets
        if(nInlets .ne. Me%NumberOfInlets)then
            write(*,*)"Number of Inlets in SWMM input file differs from Inlets config file"
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR220'
        endif

        !Double check headwalls
        if(nHeadwalls .ne. Me%NumberOfHeadwalls)then
            write(*,*)"Number of Headwalls in SWMM input file differs from Headwalls config file"
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR225'
        endif


        !Double check ponds
        if(nPonds .ne. Me%NumberOfPonds)then
            write(*,*)"Number of Ponds in SWMM input file differs from Ponds config file"
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR230'
        endif
        
        Me%ActiveNodes = Me%NumberOfManholes      + Me%NumberOfInlets    + &
                         Me%NumberOfCrossSections + Me%NumberOfHeadwalls + &
                         Me%NumberOfOutfalls      + Me%NumberOfPonds
        !Double check
        if(Me%ActiveNodes + Me%InactiveNodes .ne. Me%NumberOfNodes)then
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR250'
        endif

        write(*,*)
        write(*,*)"SewerGEMS - Nodes Inside Grid     : ", Me%ActiveNodes
        write(*,*)"SewerGEMS - Junctions             : ", Me%NumberOfManholes + Me%NumberOfInlets
        write(*,*)"SewerGEMS - Manholes              : ", Me%NumberOfManholes
        write(*,*)"SewerGEMS - Inlets                : ", Me%NumberOfInlets
        write(*,*)"SewerGEMS - Headwalls             : ", Me%NumberOfHeadwalls
        write(*,*)"SewerGEMS - Cross Sections        : ", Me%NumberOfCrossSections
        write(*,*)"SewerGEMS - Outfalls              : ", Me%NumberOfOutfalls
        write(*,*)"SewerGEMS - Ponds                 : ", Me%NumberOfPonds
        write(*,*)"SewerGEMS - Open Channel Links    : ", Me%NumberOfOpenChannelLinks
        write(*,*)"SewerGEMS - Inactive/Outside Grid : ", Me%InactiveNodes
        write(*,*)"SewerGEMS - Ignored Nodes         : ", Me%NumberOfIgnoredNodes
        write(*,*)

        do n = 1, Me%NumberOfInlets
            !iNode is inlet node index + 1
            !because MOHID nodes list starts with 1 and 
            !SWMM node list starts with 0
            iNode                           = Me%Inlets(n)%SWMM_ID + 1 
            Me%Inlets(n)%I                  = Me%NodesI(iNode)
            Me%Inlets(n)%J                  = Me%NodesJ(iNode)
            Me%Inlets(n)%CellID             = Me%NodesCellID(iNode)
            Me%Inlets(n)%EffectiveFlow      = 0.0
            Me%Inlets(n)%FlowEnteringCell   = 0.0
        enddo

        !Check for more than one inlet per grid cell
        do n = 1, Me%NumberOfInlets
            Me%Inlets(n)%nInletsInGridCell = 1
            do m = 1, Me%NumberOfInlets
                if(n .ne. m)then
                    if(Me%Inlets(n)%CellID == Me%Inlets(n)%CellID)then
                        Me%Inlets(n)%nInletsInGridCell = Me%Inlets(n)%nInletsInGridCell + 1
                        Me%Inlets(m)%nInletsInGridCell = Me%Inlets(m)%nInletsInGridCell + 1
                    endif
                endif
            enddo
        enddo
        
        nManholes = 0

        if(Me%NumberOfManholes > 0)then
            
            allocate(Me%Manholes(1:Me%NumberOfManholes))

            do n = 1, Me%NumberOfNodes
                if(Me%NodesType(n) == Manhole_)then
                    nManholes                        = nManholes + 1
                    Me%Manholes(nManholes)%SWMM_ID   = n - 1
                    Me%Manholes(nManholes)%I         = Me%NodesI(n)
                    Me%Manholes(nManholes)%J         = Me%NodesJ(n)
                    Me%Manholes(nManholes)%Name      = Me%NodesNames(n)
                    Me%Manholes(nManholes)%Outflow   = 0.0 !initialize to zero
                endif
            enddo

        endif

        if(nManholes .ne. Me%NumberOfManholes)then
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR260'
        endif

        nHeadwalls = 0

        if(Me%NumberOfHeadwalls > 0)then
            
            do n = 1, Me%NumberOfNodes
                if(Me%NodesType(n) == Headwall_)then
                    nHeadwalls                         = nHeadwalls + 1
                    Me%Headwalls(nHeadwalls)%SWMM_ID   = n - 1
                    Me%Headwalls(nHeadwalls)%I         = Me%NodesI(n)
                    Me%Headwalls(nHeadwalls)%J         = Me%NodesJ(n)
                    Me%Headwalls(nHeadwalls)%Flow      = 0.0
                    Me%Headwalls(nHeadwalls)%FlowEnteringCell = 0.0


                    STAT_CALL = SewerGEMSEngine_getNodeDownstreamLinkID(Me%Headwalls(nHeadwalls)%SWMM_ID, &
                                                                        Me%Headwalls(nHeadwalls)%SWMM_DownstreamLinkID)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR262'

                endif
            enddo

        endif

        if(nHeadwalls .ne. Me%NumberOfHeadwalls)then
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR265'
        endif

        nOutfalls = 0

        if(Me%NumberOfOutfalls > 0)then
            
            allocate(Me%Outfalls(1:Me%NumberOfOutfalls))
            
            do n = 1, Me%NumberOfNodes
                if(Me%NodesType(n) == Outfall_)then
                    nOutfalls                        = nOutfalls + 1
                    Me%Outfalls(nOutfalls)%SWMM_ID   = n - 1
                    Me%Outfalls(nOutfalls)%I         = Me%NodesI(n)
                    Me%Outfalls(nOutfalls)%J         = Me%NodesJ(n)
                    Me%Outfalls(nOutfalls)%Name      = Me%NodesNames(n)
                    Me%Outfalls(nOutfalls)%Flow      = 0.0
                endif
            enddo

        endif

        if(nOutfalls .ne. Me%NumberOfOutfalls)then
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR270'
        endif

        nCrossSections = 0

        if(Me%NumberOfCrossSections > 0)then
            
            allocate(Me%CrossSections(1:Me%NumberOfCrossSections))
            
            do n = 1, Me%NumberOfNodes
                if(Me%NodesType(n) == CrossSection_)then
                    nCrossSections                             = nCrossSections + 1
                    Me%CrossSections(nCrossSections)%SWMM_ID   = n - 1
                    Me%CrossSections(nCrossSections)%I         = Me%NodesI(n)
                    Me%CrossSections(nCrossSections)%J         = Me%NodesJ(n)
                    Me%CrossSections(nCrossSections)%Name      = Me%NodesNames(n)
                    Me%CrossSections(nCrossSections)%Flow      = 0.0
                    Me%CrossSections(nCrossSections)%WaterLevel= 0.0
                endif
            enddo

        endif

        if(nCrossSections .ne. Me%NumberOfCrossSections)then
            stop 'ConstructSewerGEMS - ModuleRunOff - ERR280'
        endif

        do n = 1, Me%NumberOfOpenChannelLinks

            Me%OpenChannelLinks(n)%LinkID = NodeIDFromName(Me%OpenChannelLinks(n)%LinkNodeName)

            if(Me%OpenChannelLinks(n)%LinkID  == null_int)then
                write(*,*)"Could not find node for open channel link"
                write(*,*)"Open channel link GRID_I : ", Me%OpenChannelLinks(n)%I
                write(*,*)"Open channel link GRID_J : ", Me%OpenChannelLinks(n)%J
                write(*,*)"SWMM node                : ", trim(adjustl(Me%OpenChannelLinks(n)%LinkNodeName))
                stop 'ConstructSewerGEMS - ModuleRunOff - ERR281'
            endif 

            if(Me%OpenChannelLinks(n)%TypeOf == Weighted_)then

                Me%OpenChannelLinks(n)%SecondLinkID = NodeIDFromName(Me%OpenChannelLinks(n)%SecondLinkNodeName)

                if(Me%OpenChannelLinks(n)%LinkID  == null_int)then
                    write(*,*)"Could not find secondary node for open channel link"
                    write(*,*)"Open channel link GRID_I : ", Me%OpenChannelLinks(n)%I
                    write(*,*)"Open channel link GRID_J : ", Me%OpenChannelLinks(n)%J
                    write(*,*)"SWMM node                : ", trim(adjustl(Me%OpenChannelLinks(n)%SecondLinkNodeName))
                    stop 'ConstructSewerGEMS - ModuleRunOff - ERR282'
                endif 

            endif

            if(Me%OpenChannelLinks(n)%TypeOf == OutfallLink_)then

                do xn = 1,  Me%NumberOfOutfalls
                    if(Me%Outfalls(xn)%SWMM_ID == Me%OpenChannelLinks(n)%LinkID)then
                        Me%OpenChannelLinks(n)%OutfallID = xn
                        exit
                    endif
                enddo

                if(Me%OpenChannelLinks(n)%OutfallID == null_int)then
                    write(*,*)"Could not find outfall link"
                    write(*,*)"Open channel link GRID_I : ", Me%OpenChannelLinks(n)%I
                    write(*,*)"Open channel link GRID_J : ", Me%OpenChannelLinks(n)%J
                    write(*,*)"SWMM node                : ", trim(adjustl(Me%OpenChannelLinks(n)%LinkNodeName))
                    write(*,*)"SWMM node is not an outfall"
                    stop 'ConstructSewerGEMS - ModuleRunOff - ERR283'
                endif

            elseif(Me%OpenChannelLinks(n)%TypeOf == PondLink_)then

                do xn = 1,  Me%NumberOfPonds
                    if(Me%Ponds(xn)%SWMM_ID == Me%OpenChannelLinks(n)%LinkID)then
                        Me%OpenChannelLinks(n)%PondID = xn
                        exit
                    endif
                enddo

                if(Me%OpenChannelLinks(n)%PondID == null_int)then
                    write(*,*)"Could not find pond link"
                    write(*,*)"Open channel link GRID_I : ", Me%OpenChannelLinks(n)%I
                    write(*,*)"Open channel link GRID_J : ", Me%OpenChannelLinks(n)%J
                    write(*,*)"SWMM node                : ", trim(adjustl(Me%OpenChannelLinks(n)%LinkNodeName))
                    write(*,*)"SWMM node is not a pond"
                    stop 'ConstructSewerGEMS - ModuleRunOff - ERR284'
                endif

            else

                do xn = 1,  Me%NumberOfCrossSections
                    if(Me%CrossSections(xn)%SWMM_ID == Me%OpenChannelLinks(n)%LinkID)then
                        Me%OpenChannelLinks(n)%CrossSectionID = xn
                        exit 
                    endif
                enddo

                if(Me%OpenChannelLinks(n)%CrossSectionID == null_int)then
                    write(*,*)"Could not find cross-section open channel link"
                    write(*,*)"Open channel link GRID_I : ", Me%OpenChannelLinks(n)%I
                    write(*,*)"Open channel link GRID_J : ", Me%OpenChannelLinks(n)%J
                    write(*,*)"SWMM node                : ", trim(adjustl(Me%OpenChannelLinks(n)%LinkNodeName))
                    write(*,*)"SWMM node is not a cross-section"
                    stop 'ConstructSewerGEMS - ModuleRunOff - ERR285'
                endif

            endif

        enddo

        !Get SewerGEMS SWMM current time step 
        STAT_CALL = SewerGEMSEngine_getTotalVolume(Me%Total1DVolume)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR289'

        Me%MaxTotal1DVolume = Me%Total1DVolume

        STAT_CALL = SewerGEMSEngine_getdt(Me%StormWaterModelDT)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSewerGEMS - ModuleRunOff - ERR290'


        write(*,*)
        write(*,*)"SewerGEMS - Successfully constructed"
        write(*,*)

        deallocate(Me%NodesX      )
        deallocate(Me%NodesY      )
        deallocate(Me%NodesI      )
        deallocate(Me%NodesJ      )
        deallocate(Me%NodesID     )
        deallocate(Me%NodesCellID )
        deallocate(Me%NodesNames  )
        deallocate(Me%NodesType   )
    
    end subroutine ConstructSewerGEMS

    !--------------------------------------------------------------------------

    logical function IgnoreSWMMNode(Name)
    
        character(len = 99), intent(in) :: Name
        integer                         :: n
        
        IgnoreSWMMNode = .false.
        do n = 1, Me%NumberOfIgnoredNodes
            if(trim(adjustl(Name)) == trim(adjustl(Me%IgnoredNodes(n))))then
                IgnoreSWMMNode           = .true.
                exit 
            end if 
        enddo
    
    end function IgnoreSWMMNode

    !--------------------------------------------------------------------------
    
    integer function ValidateInlet(Name, NodeID)
    
        character(len = 99), intent(in) :: Name
        integer,             intent(in) :: NodeID
        integer                         :: n
        
        ValidateInlet = null_int
        do n = 1, Me%NumberOfInlets
            if(trim(adjustl(Name)) == trim(adjustl(Me%Inlets(n)%Name)))then
                ValidateInlet           = Me%Inlets(n)%ID
                Me%Inlets(n)%SWMM_ID    = NodeID
                exit 
            end if 
        enddo
    
    end function ValidateInlet

    !--------------------------------------------------------------------------
    
    integer function ValidateHeadwall(Name, NodeID)
    
        character(len = 99), intent(in) :: Name
        integer,             intent(in) :: NodeID
        integer                         :: n
        
        ValidateHeadwall = null_int
        do n = 1, Me%NumberOfHeadwalls
            if(trim(adjustl(Name)) == trim(adjustl(Me%Headwalls(n)%Name)))then
                ValidateHeadwall        = Me%Headwalls(n)%ID
                Me%Headwalls(n)%SWMM_ID = NodeID
                exit 
            end if 
        enddo
    
    end function ValidateHeadwall


    !--------------------------------------------------------------------------

    integer function NodeIDFromName(Name)
    
        character(len = 99), intent(in) :: Name
        integer                         :: n

        NodeIDFromName = null_int

        do n = 1, Me%NumberOfNodes
            if(trim(adjustl(Name)) == trim(adjustl(Me%NodesNames(n))))then
                NodeIDFromName  = Me%NodesID(n)
                exit 
            end if 
        enddo
    
    end function NodeIDFromName

    !--------------------------------------------------------------------------
    
    integer function ValidatePond(Name, NodeID)
    
        character(len = 99), intent(in) :: Name
        integer,             intent(in) :: NodeID
        integer                         :: n
        
        ValidatePond = null_int
        do n = 1, Me%NumberOfPonds
            if(trim(adjustl(Name)) == trim(adjustl(Me%Ponds(n)%Name)))then
                ValidatePond           = Me%Ponds(n)%ID
                Me%Ponds(n)%SWMM_ID    = NodeID
                exit 
            end if 
        enddo
    
    end function ValidatePond


#endif _SEWERGEMSENGINECOUPLER_

    !--------------------------------------------------------------------------

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
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR04'

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR05'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",                   &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR06'

        !Writes the River Points
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "RiverPoints", "-",                   &
                              Array2D = Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR07'


        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR08'


    end subroutine ConstructHDF5Output

    !--------------------------------------------------------------------------
    
    subroutine ConstructTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        integer                                             :: nProperties
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=PathLength)                           :: TimeSerieLocationFile
        integer                                             :: TimeSerieNumber, dn, Id, Jd
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        character(len=StringLength)                         :: TimeSerieName
        integer                                             :: iProperty
        
        !Begin------------------------------------------------------------------

        nProperties = 8
        iProperty = 8
        !if(Me%StormWaterModel)then
        !    nProperties = nProperties + 4
        !endif      
        if (Me%Use1D2DInteractionMapping) then
            nProperties = nProperties + 4
        endif

        !Allocates PropertyList
        allocate(PropertyList(nProperties))
        
        !Property names
        PropertyList(1) = trim(GetPropertyName (WaterLevel_)) 
        PropertyList(2) = trim(GetPropertyName (WaterColumn_)) 
        PropertyList(3) = "flow X" 
        PropertyList(4) = "flow_Y" 
        PropertyList(5) = trim(GetPropertyName (FlowModulus_)) 
        PropertyList(6) = trim(GetPropertyName (VelocityU_))
        PropertyList(7) = trim(GetPropertyName (VelocityV_))
        PropertyList(8) = trim(GetPropertyName (VelocityModulus_))
      
        !if(Me%StormWaterModel)then
        !    iProperty = iProperty + 1
        !    PropertyList(iProperty)  = "storm water potential flow"
        !    iProperty = iProperty + 1
        !    PropertyList(iProperty) = "storm water effective flow"
        !    iProperty = iProperty + 1
        !    PropertyList(iProperty) = "street gutter potential flow"
        !    iProperty = iProperty + 1
        !    PropertyList(iProperty) = "street gutter effective flow"
        !endif
        
        if (Me%Use1D2DInteractionMapping) then
            iProperty = iProperty + 1
            PropertyList(iProperty)  = "node river level"
            iProperty = iProperty + 1
            PropertyList(iProperty) = "margin river level"
            iProperty = iProperty + 1
            PropertyList(iProperty) = "margin flow to channels"
            iProperty = iProperty + 1
            PropertyList(iProperty) = "node flow to channels"            
        endif
        
        call GetData (TimeSerieLocationFile,                  &
                      Me%ObjEnterData, iflag,                 &
                      SearchType   = FromFile,                &
                      keyword      = 'TIME_SERIE_LOCATION',   &
                      ClientModule = 'ModuleRunoff',          &
                      Default      = Me%Files%DataFile,       &
                      STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR010' 

        if (iflag == 1) then
            Me%OutPut%TimeSeries = .true.
        else
            Me%OutPut%TimeSeries = .false.
        endif
        
        !Constructs TimeSerie
        call StartTimeSerie (Me%ObjTimeSerie, Me%ObjTime,           &
                            TimeSerieLocationFile,                 &
                            PropertyList, "srr",                   &
                            WaterPoints2D = Me%ExtVar%BasinPoints, &
                            STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR030' 

        !Deallocates PropertyList
        deallocate(PropertyList)

        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR050'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation (Me%ObjTimeSerie, dn, &  
                                    CoordX   = CoordX,   &
                                    CoordY   = CoordY,   & 
                                    CoordON  = CoordON,  &
                                    STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR060'
            
            call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR070'
            
    i1:     if (CoordON) then
    
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR080'

                if (Id < 0 .or. Jd < 0) then
                
                   call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR090'

                   if (IgnoreOK) then
                      write(*,*) 'Time Serie outside the domain - ',trim(TimeSerieName)
                      cycle
                   else
                      stop 'ConstructTimeSeries - ModuleRunoff - ERR100'
                   endif

                endif

                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR110'

            endif i1

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,    &  
                                      LocalizationI   = Id,   &
                                      LocalizationJ   = Jd,   & 
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleRunoff - ERR120'

            if (Me%ExtVar%BasinPoints(Id, Jd) /= WaterPoint) then
                write(*,*) 'Time Serie in a cell outside basin - ',trim(TimeSerieName)
            endif

        enddo
       
    end subroutine ConstructTimeSeries

    !--------------------------------------------------------------------------
    
    subroutine ReadInitialFile_Bin

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: InitialFile
        type (T_Time)                               :: BeginTime, EndTimeFile, EndTime
        real                                        :: DT_error
        integer                                     :: STAT_CALL, i, j

        !----------------------------------------------------------------------

        call UnitsManager(InitialFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFileOld - ModuleRunoff - ERR01'

        open(Unit = InitialFile, File = Me%Files%InitialFile, Form = 'UNFORMATTED',     &
             status = 'OLD', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFileOld - ModuleRunoff - ERR02'

        !Reads Date
        read(InitialFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
        call SetDate(EndTimeFile, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFileOld - ModuleRunoff - ERR03'
        
        DT_error = EndTimeFile - BeginTime

        !Avoid rounding erros - Frank 08-2001
        !All runs are limited to second definition - David 10-2015     
        !if (abs(DT_error) >= 0.01) then
        if (abs(DT_error) >= 1) then
            write(*,*) 'The end time of the previous run is different from the start time of this run'
            write(*,*) 'Date in the file'
            write(*,*) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
            write(*,*) 'DT_error', DT_error
            if (Me%StopOnWrongDate) stop 'ReadInitialFileOld - ModuleRunoff - ERR04'   

        endif

        read(InitialFile)Me%myWaterColumn
        
   10   continue

        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFileOld - ModuleRunoff - ERR05'  
        
        !Updates Volume & Level from Column        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
            if (Me%ExtVar%BasinPoints(i, j) == 1) then
            
                Me%myWaterLevel(i, j)   = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
                Me%myWaterVolume(i, j)  = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
                
            endif

        enddo
        enddo      
      
    end subroutine ReadInitialFile_Bin
    
    !--------------------------------------------------------------------------
    
    subroutine ReadInitialFile_Hdf()


        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: EXIST
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ
        integer                                     :: i, j
        type (T_Time)                               :: BeginTime, EndTimeFile, EndTime
        real, dimension(:), pointer                 :: TimePointer
        real                                        :: DT_error        

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        !----------------------------------------------------------------------

        inquire (FILE=trim(Me%Files%InitialFile), EXIST = Exist)

cd0:    if (Exist) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%InitialFile),                              &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialFile - ModuleRunoff - ERR01'

            
            !Get Time
            call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialFile - ModuleRunoff - ERR010'             
            
            allocate(TimePointer(1:6))
            call HDF5ReadData   (ObjHDF5, "/Time",                                       &
                                    "Time",                                              &
                                    Array1D = TimePointer,                               &
                                    STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialFile - ModuleRunoff - ERR020'             
            
            
            call SetDate(EndTimeFile, TimePointer(1), TimePointer(2),    &
                                      TimePointer(3), TimePointer(4),    &
                                      TimePointer(5), TimePointer(6))
            
            
            call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, EndTime = EndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleRunoff - ERR030'
        
            DT_error = EndTimeFile - BeginTime

            !Avoid rounding erros - Frank 08-2001
            !All runs are limited to second definition - David 10-2015
            !if (abs(DT_error) >= 0.01) then
            if (abs(DT_error) >= 1) then
            
                write(*,*) 'The end time of the previous run is different from the start time of this run'
                write(*,*) 'Date in the file'
                write(*,*) TimePointer(1), TimePointer(2), TimePointer(3), TimePointer(4), TimePointer(5), TimePointer(6)
                write(*,*) 'DT_error', DT_error
                if (Me%StopOnWrongDate) stop 'ReadInitialFile - ModuleRunoff - ERR040'   

            endif                  
            deallocate(TimePointer)
            

            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB,                                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialFile - ModuleRunoff - ERR050'

            call HDF5ReadData   (ObjHDF5, "/Results/water column",                       &
                                 "water column",                                         &
                                 Array2D = Me%myWaterColumn,                             &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialFile - ModuleRunoff - ERR060'
            
            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialFile - ModuleRunoff - ERR080'
            
            
            
            !Updates Volume & Level from Column        
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                if (Me%ExtVar%BasinPoints(i, j) == 1) then
            
                    Me%myWaterLevel(i, j)   = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
                    Me%myWaterVolume(i, j)  = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
                
                endif

            enddo
            enddo               
            

        else
            
            write(*,*)
            stop 'ReadInitialFile - ModuleRunoff - ERR090'

        end if cd0

    end subroutine ReadInitialFile_Hdf

    !--------------------------------------------------------------------------    
 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
    subroutine GetOverLandFlow (ObjRunOffID, FlowX, FlowY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: FlowX, FlowY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowX => Me%iFlowX

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowY => Me%iFlowY

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetOverLandFlow

    !--------------------------------------------------------------------------

    subroutine GetFlowToChannels (ObjRunOffID, FlowToChannels, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: FlowToChannels
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowToChannels => Me%iFlowToChannels

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetFlowToChannels

    !--------------------------------------------------------------------------

    subroutine GetBoundaryImposed (ObjRunOffID, BoundaryOpen, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        logical, intent(OUT)                            :: BoundaryOpen
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BoundaryOpen = Me%ImposeBoundaryValue

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetBoundaryImposed
    
    !--------------------------------------------------------------------------

    subroutine GetRouteDFour (ObjRunOffID, RouteD4, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        logical, intent(OUT)                            :: RouteD4
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            RouteD4 = Me%RouteDFourPoints

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetRouteDFour

    !--------------------------------------------------------------------------   

    subroutine GetRouteDFourCells (ObjRunOffID, RouteD4Cells, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        integer, pointer, dimension (:,:)               :: RouteD4Cells
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            RouteD4Cells => Me%DFourSinkPoint

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetRouteDFourCells

    !--------------------------------------------------------------------------   

    subroutine GetRouteDFourNeighbours (ObjRunOffID, RouteD4LowerI, RouteD4LowerJ, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        integer, pointer, dimension (:,:)               :: RouteD4LowerI
        integer, pointer, dimension (:,:)               :: RouteD4LowerJ
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            RouteD4LowerI => Me%LowestNeighborI
            RouteD4LowerJ => Me%LowestNeighborJ

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetRouteDFourNeighbours

    !-------------------------------------------------------------------------- 

    subroutine GetRouteDFourFlux (ObjRunOffID, DFourFlow, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: DFourFlow
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            DFourFlow => Me%iFlowRouteDFour

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRouteDFourFlux
    
    !--------------------------------------------------------------------------

    subroutine GetBoundaryFlux (ObjRunOffID, FlowAtBoundary, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: FlowAtBoundary
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowAtBoundary => Me%iFlowBoundary

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoundaryFlux
    
    !--------------------------------------------------------------------------

    subroutine GetBoundaryCells (ObjRunOffID, BoundaryCells, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        integer,   pointer, dimension(:,:)              :: BoundaryCells
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BoundaryCells => Me%BoundaryCells

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoundaryCells

    !--------------------------------------------------------------------------

    subroutine GetFlowDischarge (ObjRunOffID, FlowDischarge, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: FlowDischarge
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowDischarge => Me%iFlowDischarge

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetFlowDischarge
    
    !--------------------------------------------------------------------------

    subroutine GetRunOffTotalDischargeFlowVolume (ObjRunOffID, Volume, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8)                                         :: Volume
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Volume = Me%TotalDischargeFlowVolume

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
    
    end subroutine GetRunOffTotalDischargeFlowVolume
    
    !--------------------------------------------------------------------------

        subroutine GetRunoffWaterLevel (ObjRunOffID, Waterlevel, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterLevel
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            WaterLevel => Me%MyWaterLevel

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffWaterLevel

    !--------------------------------------------------------------------------
    
    subroutine GetRunoffWaterColumn (ObjRunOffID, WaterColumn, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterColumn
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            WaterColumn => Me%MyWaterColumn

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffWaterColumn

    !--------------------------------------------------------------------------

    subroutine GetRunoffWaterColumnOld (ObjRunOffID, WaterColumnOld, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterColumnOld
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            WaterColumnOld => Me%MyWaterColumnOld

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffWaterColumnOld

    !--------------------------------------------------------------------------
    
    subroutine GetRunoffWaterColumnAT (ObjRunOffID, WaterColumn, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterColumn
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            !because runoff water may go all to the river in one time step and Runoff Conc would be zero
            !but not flow, transport has to be separated from drainage network interaction
            !and explicit/implicit transport is only evaluated with water column after transport            
            WaterColumn => Me%MyWaterColumnAfterTransport

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffWaterColumnAT    

    !--------------------------------------------------------------------------

    subroutine GetRunoffCenterVelocity (ObjRunOffID, VelX, VelY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: VelX, VelY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            VelX => Me%CenterVelocityX

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            VelY => Me%CenterVelocityY

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffCenterVelocity

    !--------------------------------------------------------------------------
    
    subroutine GetManning (ObjRunOffID, Manning, ManningX, ManningY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer, optional        :: Manning, ManningX, ManningY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (present(Manning)) then
                call Read_Lock(mRUNOFF_, Me%InstanceID)
                Manning => Me%OverlandCoefficient
            endif

            if (present(ManningX)) then
                call Read_Lock(mRUNOFF_, Me%InstanceID)
                ManningX => Me%OverlandCoefficientX
            endif

            if (present(ManningY)) then
                call Read_Lock(mRUNOFF_, Me%InstanceID)
                ManningY => Me%OverlandCoefficientY
            endif
            
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetManning

    !--------------------------------------------------------------------------    


    subroutine GetManningDelta (ObjRunOffID, ManningDelta, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: ManningDelta
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mRUNOFF_, Me%InstanceID)
            ManningDelta => Me%OverlandCoefficientDelta

           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetManningDelta

    !--------------------------------------------------------------------------   

    subroutine GetMassError (ObjRunOffID, MassError, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: MassError
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            MassError => Me%MassError

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetMassError
    
    !--------------------------------------------------------------------------
    
    subroutine GetRunoffTotalStoredVolume (ObjRunoffID, TotalStoredVolume, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunoffID
        real(8)                                         :: TotalStoredVolume
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunoffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            TotalStoredVolume = Me%TotalStoredVolume

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffTotalStoredVolume

    !-------------------------------------------------------------------------
    
    subroutine GetRunOffStoredVolumes (ID, Surface, StormSystem, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: ID
        real(8), intent(OUT), optional                  :: Surface
        real(8), intent(OUT), optional                  :: StormSystem
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        !Begin-----------------------------------------------------------------
        call Ready(ID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Surface))     Surface     = Me%VolumeStoredInSurface
            if (present(StormSystem)) StormSystem = Me%VolumeStoredInStormSystem

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_            
        !----------------------------------------------------------------------
    
    end subroutine GetRunOffStoredVolumes
    
    !-------------------------------------------------------------------------
    
    subroutine GetRunOffBoundaryFlowVolume (ID, BoundaryFlowVolume, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: ID
        real(8), intent(OUT)                            :: BoundaryFlowVolume
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        !Begin-----------------------------------------------------------------
        call Ready(ID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BoundaryFlowVolume = Me%BoundaryFlowVolume

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_            
        !----------------------------------------------------------------------    
    
    end subroutine GetRunOffBoundaryFlowVolume
    
    !-------------------------------------------------------------------------
        
    subroutine GetNextRunOffDT (ObjRunOffID, DT, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, intent(OUT)                               :: DT
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(ObjRunOffID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DT        = Me%CV%NextDT

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetNextRunOffDT

    !--------------------------------------------------------------------------

    subroutine SetBasinColumnToRunoff(ObjRunOffID, WaterColumnOld, WaterColumn, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterColumnOld
        real(8), dimension(:, :), pointer               :: WaterColumn
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL, ready_
        integer                                         :: i, j
        integer                                         :: ILB, IUB, JLB, JUB, CHUNK

        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        call Ready(ObjRunOffID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
        
            !Actualizes water column, water level and water volume
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR01'

            !Gets a pointer to Topography
            call GetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR10'
        
            call GetGridCellArea  (Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea,             &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR020'
        
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                    Me%myWaterColumnOld(i, j) = WaterColumnOld(i, j)
                    Me%myWaterColumn(i, j)    = WaterColumn(i, j)

                    Me%myWaterLevel (i, j) = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
                    !Here the water column is the uniformly distributed one. Inside 
                    Me%myWaterVolume(i, j) = WaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)

!                    Me%myWaterVolume(i, j) = WaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
!                    Me%myWaterColumn(i, j) = Me%myWaterVolume(i, j) / Me%FreeGridCellArea(i, j)
!                    Me%myWaterLevel (i, j) = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
                   
                endif
            enddo
            enddo            
            !$OMP END DO
            !$OMP END PARALLEL

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR030'

            !Ungets the Topography
            call UngetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR40'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR050'

            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
                    
    end subroutine SetBasinColumnToRunoff

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine SetBasinStatsToRunoff(ObjRunOffID, TotalCumulativeRainfallVolume,        &
                                                  TotalCumulativeInfiltrationVolume,    &
                                                  AverageCumulativeInfiltrationVolume,  &
                                                  AverageCumulativeInfiltrationDepth,   &
                                                  STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8)                                         :: TotalCumulativeRainfallVolume
        real(8)                                         :: TotalCumulativeInfiltrationVolume
        real(8)                                         :: AverageCumulativeInfiltrationVolume
        real(8)                                         :: AverageCumulativeInfiltrationDepth
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_


        call Ready(ObjRunOffID, ready_)

        if ((ready_ .EQ. IDLE_ERR_)) then
        
            Me%TotalRainfallVolume          = TotalCumulativeRainfallVolume
            Me%TotalInfiltrationVolume      = TotalCumulativeInfiltrationVolume
            Me%AvrgAccInfiltrationVolume    = AverageCumulativeInfiltrationVolume
            Me%AvrgAccInfiltrationDepth     = AverageCumulativeInfiltrationDepth
            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
                    
    end subroutine SetBasinStatsToRunoff

    !--------------------------------------------------------------------------

    subroutine UnGetRunOff2D_R4(ObjRunOffID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunOffID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRUNOFF_, Me%InstanceID, "UnGetRunOff2D_R4")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunOff2D_R4

    !--------------------------------------------------------------------------

    subroutine UnGetRunOff2D_R8(ObjRunOffID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunOffID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRUNOFF_, Me%InstanceID, "UnGetRunOff2D_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunOff2D_R8
        
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine ModifyRunOff(RunOffID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: RunOffID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL
        real                                        :: SumDT
        logical                                     :: Restart
        integer                                     :: Niter, iter
        integer                                     :: n_restart
        logical                                     :: IsFinalFile
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunOffID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ModifyRunOff")

            !Time Stuff
            call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyRunOff - ModuleRunOff - ERR10'

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyRunOff - ModuleRunOff - ERR11'

            
            !Stores initial values = from basin
            call SetMatrixValue(Me%myWaterColumnOld, Me%Size, Me%myWaterColumn)
            call SetMatrixValue(Me%InitialFlowX,     Me%Size, Me%iFlowX)
            call SetMatrixValue(Me%InitialFlowY,     Me%Size, Me%iFlowY)            
          
            
            !Set 1D River level in river boundary cells
            !From External model or DN
            if (Me%Use1D2DInteractionMapping) then
                call InterpolateRiverLevelToCells
            endif
            
            Restart = .true.
            n_restart = 0


            
            if (Me%CV%NextNiteration > 1 .and. Me%ExtVar%DT < (Me%CV%CurrentDT * Me%CV%NextNiteration)) then
                Me%CV%NextNiteration = max(aint(Me%ExtVar%DT / Me%CV%CurrentDT), 1.0)
            endif            
            
            do while (Restart)
            
                !Calculates local Watercolumn
                call ReadLockExternalVar   (StaticOnly = .true.)
                call LocalWaterColumn      (Me%myWaterColumnOld)
                call ReadUnLockExternalVar (StaticOnly = .true.)

                SumDT        = 0.0
                Restart      = .false.                                 
                iter         = 1
                Niter        = Me%CV%NextNiteration    !DB
                Me%CV%CurrentDT = Me%ExtVar%DT / Niter

                if (Niter > 1) then                
                    call WriteDTLog_ML ('ModuleRunOff', Niter, Me%CV%CurrentDT)
                endif
                
                call SetMatrixValue(Me%iFlowX, Me%Size, dble(0.0))
                call SetMatrixValue(Me%iFlowY, Me%Size, dble(0.0))
                call SetMatrixValue(Me%lFlowX, Me%Size, Me%InitialFlowX)
                call SetMatrixValue(Me%lFlowY, Me%Size, Me%InitialFlowY)
                call SetMatrixValue(Me%iFlowToChannels, Me%Size, 0.0)
                call SetMatrixValue(Me%iFlowBoundary, Me%Size, 0.0)
                call SetMatrixValue(Me%iFlowRouteDFour, Me%Size, 0.0)
                
doIter:         do while (iter <= Niter)

                    !Gets ExternalVars
                    call ReadLockExternalVar (StaticOnly = .false.)

                    !Stores WaterVolume for convergence test
                    call SetMatrixValue(Me%myWaterVolumeOld, Me%Size, Me%myWaterVolume)

                    call SetMatrixValue(Me%FlowXOld,         Me%Size, Me%lFlowX)
                    call SetMatrixValue(Me%FlowYOld,         Me%Size, Me%lFlowY)

                    !Updates Geometry
                    call ModifyGeometryAndMapping
                    
                    !save most recent water volume to predict if negative occur. in that case flux will be
                    !limited to water volume and next fluxes will be zero
                    if (.not. Me%LimitToCriticalFlow) then
                        call SetMatrixValue (Me%myWaterVolumePred, Me%Size, Me%myWaterVolume)
                    endif        
                    
                    select case (Me%HydrodynamicApproximation)
                        case (KinematicWave_)
                            call KinematicWave  ()            !Slope based on topography
                        case (DiffusionWave_)
                            call KinematicWave  ()            !Slope based on surface
                        case (DynamicWave_)
                            call ComputeFaceVelocityModulus
                            call DynamicWaveXX    (Me%CV%CurrentDT)   !Consider Advection, Friction and Pressure
                            call DynamicWaveYY    (Me%CV%CurrentDT)
                    end select

                    !Updates waterlevels, based on fluxes
                    call UpdateWaterLevels(Me%CV%CurrentDT)
                    
                    !Interaction with channels
                    if (.not. Me%Use1D2DInteractionMapping .and. Me%ObjDrainageNetwork /= 0 .and. .not. Me%SimpleChannelInteraction) then
                        call FlowIntoChannels       (Me%CV%CurrentDT)
                    endif

                    !Boundary Condition
!                    if (Me%ImposeBoundaryValue) then
!                        call ImposeBoundaryValue    (Me%CV%CurrentDT)
!                    endif

                    !Inputs Water from discharges
                    if (Me%Discharges) then
                        call ModifyWaterDischarges  (Me%CV%CurrentDT)                
                    endif

                    call CheckStability(Restart) 
                    
                    call ReadUnLockExternalVar (StaticOnly = .false.)
                    
                    if (Restart) then
                        exit doIter
                    endif

                    call IntegrateFlow     (Me%CV%CurrentDT, SumDT)  
                        
                    SumDT = SumDT + Me%CV%CurrentDT
                    iter  = iter  + 1                                        
                    
                enddo doIter
                                            
            enddo
            
!            !DB
!            if (Niter <= Me%LastGoodNiter) then
!                Me%CV%NextNiteration = max (min(int(Niter / Me%InternalTimeStepSplit), NIter - 1), 1)
!            else
!                Me%CV%NextNiteration = Niter
!            endif     
            
            !save water column before removes from next processes
            !important for property transport if river cells get out of water, conc has to be computed
            !after transport and not zero because there was no water left
            call SetMatrixValue(Me%myWaterColumnAfterTransport, Me%Size, Me%myWaterColumn)
            
            !Gets ExternalVars
            call ReadLockExternalVar (StaticOnly = .false.)

            if (Me%Use1D2DInteractionMapping) then
                !it will use mapping for any model (DN or SWMM or other)
                call OverLandChannelInteraction_6_NewMapping            
            else            
                if (Me%ObjDrainageNetwork /= 0) then
                
                    if (Me%SimpleChannelInteraction) then
                        !There were many methods to do this calculation that were implemented and abandoned
                        !Many of the abandoned methods were commented:
                        !OverLandChannelInteraction, OverLandChannelInteraction, 3, 4, 5, New
                        !In June/July 2022 they were deleted in a code clean up.
                        !Anyone trying to improved this code should look at code versions from around June 2022
                        !and check those routines as reference 
                    
                        !TODO: Remove this and have only one method!! This is a workaround
                        !Try to use the in _6 the code from _2 where water exits the river (celerity limited)
                        if (Me%ChannelHasTwoGridPoints) then
                            !new method adapted to hydraulic simulation with two grid points per river node (e.g. UKBenchmark Test7)
                            call OverLandChannelInteraction_6
                        else
                            !old method that is stable on normal cases only one cell per river node (e.g. Trancao sample generates "rendilhado"
                            !when water exits the river if _6 method is used but is stable with this one _2)
                            call OverLandChannelInteraction_2
                        endif
                    else
                        !Calculates flow from channels to land -> First ever implement approach
                        call FlowFromChannels
                    endif
                endif
            endif
            
            if (Me%StormWaterModel) then
                
                call ComputeStormWaterModel

                call ModifyGeometryAndMapping
                
            endif

            !Routes Ponded levels which occour due to X/Y direction (Runoff does not route in D8)
            !the defaul method was celerity (it was corrected) but it ccould create high flow changes. Manning method is stabler
            !because of resistance. However in both methods the area used is not consistent (regular faces flow
            !already used all the cell vertical areas and the route D4 will overlapp areas - review this in the future
            if (Me%RouteDFourPoints) then
                if (Me%RouteDFourMethod == Manning_) then
                    call RouteDFourPoints
                elseif (Me%RouteDFourMethod == Celerity_) then
                    call RouteDFourPoints_v3
                endif
            endif

            !Boundary Condition
            !Only compute if case of waterlevel higher than boundary (overflow)
            !the default method was instantaneous flow (instantaneous go to boundary level)
            !but it was changed to compute flow (based on celerity) to be more consistent
            !with a free drop to boundary level (that can be much lower than topography)
            if (Me%ImposeBoundaryValue) then

                if (Me%BoundaryImposedLevelInTime)then
                    call ModifyBoundaryLevel
                endif
                
                if (Me%BoundaryMethod == ComputeFlow_) then
                    call ImposeBoundaryValue
                elseif (Me%BoundaryMethod == InstantaneousFlow_) then
                    call ImposeBoundaryValue_v2
                endif
            endif

            !Calculates center flow and velocities (for output and next DT)
            call ComputeCenterValues

            call ComputeNextDT (Niter)                                    
            
            !Output Results
            if (Me%OutPut%Yes) then                   
                call RunOffOutput
            endif
            
            if(Me%OutPut%TimeSeries) then
                call OutputTimeSeries
            endif
            
            if (Me%Output%BoxFluxes) then
                call ComputeBoxesWaterFluxes
            endif

            if (Me%Output%WriteMaxWaterColumn .or. Me%Output%WriteMaxFloodRisk) then
                call OutputFlooding
            endif
            
            if (Me%Output%WriteFloodPeriod) then            
                call OutputFloodPeriod            
            endif  

            if (Me%Output%WriteFloodArrivalTime) then            
                call OutputFloodArrivalTime            
            endif              
            
            call CalculateTotalStoredVolume

            !Restart Output
            if (Me%Output%WriteRestartFile .and. .not. (Me%ExtVar%Now == Me%EndTime)) then
                if(Me%ExtVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                    IsFinalFile = .false.
                    if (Me%OutPut%RestartFormat == BIN_) then
                        call WriteFinalFile_Bin(IsFinalFile)
                    else if (Me%OutPut%RestartFormat == HDF_) then
                        call WriteFinalFile_Hdf(IsFinalFile)
                    endif
                    Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
                endif
            endif

            !Ungets external variables
            call ReadUnLockExternalVar (StaticOnly = .false.)

            STAT_ = SUCCESS_
            if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ModifyRunOff")

        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
        
    end subroutine ModifyRunOff
    
    !---------------------------------------------------------------------------
    
    subroutine InterpolateRiverLevelToCells
    
        !Argumnets--------------------------------------------------------------
    
        !Local------------------------------------------------------------------
        type(T_NodeGridPoint), pointer                  :: NodeGridPoint
        type(T_BankGridPoint), pointer                  :: BankGridPoint, BankGridPointUp, BankGridPointDown
        type(T_MarginGridPoint), pointer                :: MarginGridPoint
        real, dimension(:,:), pointer                   :: ChannelsWaterLevel
        integer                                         :: STAT_CALL
        !logical                                         :: Found, FoundUp, FoundDown
    
        !Begin------------------------------------------------------------------
    
        !output
        call SetMatrixValue(Me%MarginRiverLevel, Me%Size, null_real)
        
        !if using DN get water level from module and fill matrix. If using SWMM this will be filled by OpenMI
        !in case DN level appears in river points (in case of two banks are the bank grid points) and do not need
        !the node grid points
        if (Me%ObjDrainageNetwork /= 0) then
            
            call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InterpolateRiverLevelToCells - ModuleRunOff - ERR01'     

            !Set the matrix from DN
            call SetMatrixValue(Me%NodeRiverLevel, Me%Size, ChannelsWaterLevel)
            
            call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InterpolateRiverLevelToCells - ModuleRunOff - ERR05'            
            
            !2. Go directly to bank grid points and skip node point
            BankGridPoint => Me%FirstBankGridPoint            
            do while (associated(BankGridPoint))                
                BankGridPoint%RiverLevel = Me%NodeRiverLevel (BankGridPoint%GridI, BankGridPoint%GridJ)                    
                BankGridPoint => BankGridPoint%Next
            enddo 
            
        else            
            !1.        
            !Node Level from the matrix
            NodeGridPoint => Me%FirstNodeGridPoint
            do while (associated(NodeGridPoint))
                NodeGridPoint%RiverLevel = Me%NodeRiverLevel (NodeGridPoint%GridI, NodeGridPoint%GridJ)
                NodeGridPoint => NodeGridPoint%Next
            enddo
            
            !2.    
            !Bank river level from node
            BankGridPoint => Me%FirstBankGridPoint            
            do while (associated(BankGridPoint))                
                BankGridPoint%RiverLevel = Me%NodeGridPointArray(BankGridPoint%NGPIdidx)%ptr%RiverLevel
                    
                BankGridPoint => BankGridPoint%Next                
            enddo              
        endif    
        
        !3.
        !Margin river level from bank (interpolation)
        MarginGridPoint => Me%FirstMarginGridPoint            
        do while (associated(MarginGridPoint))
            BankGridPointUp => Me%BankGridPointArray(MarginGridPoint%BGPUpIdidx)%ptr
            BankGridPointDown => Me%BankGridPointArray(MarginGridPoint%BGPDownIdidx)%ptr
            !points that do not have interaction with river
            if (BankGridPointUp%RiverLevel < null_real / 2.0 .or. BankGridPointDown%RiverLevel < null_real / 2.0) then
                MarginGridPoint%RiverLevel = null_real
            else
                MarginGridPoint%RiverLevel = BankGridPointUp%RiverLevel - (BankGridPointUp%RiverLevel - BankGridPointDown%RiverLevel) * MarginGridPoint%InterpolationFraction                   
            endif
            !output
            Me%MarginRiverLevel(MarginGridPoint%GridI, MarginGridPoint%GridJ) = MarginGridPoint%RiverLevel
                    
            MarginGridPoint => MarginGridPoint%Next                
        enddo
        
    end subroutine InterpolateRiverLevelToCells
    
    !---------------------------------------------------------------------------
    
    subroutine FindNodeGridPoint (NodeID, NodeGridPoint, Found)

        !Arguments--------------------------------------------------------------
        integer, intent(IN)                             :: NodeID
        type (T_NodeGridPoint), pointer, intent(OUT)    :: NodeGridPoint
        logical, intent(OUT)                            :: Found
        !Local------------------------------------------------------------------


        Found = .FALSE.
        
        nullify(NodeGridPoint)
        NodeGridPoint => Me%FirstNodeGridPoint
        
        do while (associated(NodeGridPoint))
        
            if (NodeGridPoint%ID == NodeID) then
                Found = .TRUE.
                exit
            end if
            NodeGridPoint => NodeGridPoint%Next
        end do

    end subroutine FindNodeGridPoint

    !---------------------------------------------------------------------------      
    
    subroutine FindBankGridPoint (BankID, BankGridPoint, Found)

        !Arguments--------------------------------------------------------------
        integer, intent(IN)                             :: BankID
        type (T_BankGridPoint), pointer, intent(OUT)    :: BankGridPoint
        logical, intent(OUT)                            :: Found
        !Local------------------------------------------------------------------


        Found = .FALSE.
        
        nullify(BankGridPoint)
        BankGridPoint => Me%FirstBankGridPoint
        
        do while (associated(BankGridPoint))
        
            if (BankGridPoint%ID == BankID) then
                Found = .TRUE.
                exit
            end if
            BankGridPoint => BankGridPoint%Next
        end do

    end subroutine FindBankGridPoint

    !---------------------------------------------------------------------------       
    
    subroutine FindMarginGridPoint (i, j, Found)

        !Arguments--------------------------------------------------------------
        integer, intent(IN)                             :: i
        integer, intent(IN)                             :: j        
        logical, intent(OUT)                            :: Found
        !Local------------------------------------------------------------------
        type (T_MarginGridPoint), pointer               :: MarginGridPoint

        Found = .FALSE.
        
        nullify(MarginGridPoint)
        MarginGridPoint => Me%FirstMarginGridPoint
        
        do while (associated(MarginGridPoint))
        
            if (MarginGridPoint%GridI == i .and. MarginGridPoint%GridJ == j) then
                Found = .TRUE.
                exit
            end if
            MarginGridPoint => MarginGridPoint%Next
        end do

    end subroutine FindMarginGridPoint

    !---------------------------------------------------------------------------       
    

    subroutine ModifyWaterDischarges (LocalDT)

        !Arguments--------------------------------------------------------------
        real                                    :: LocalDT

        !Local------------------------------------------------------------------
        integer                                 :: iDis, nDischarges, nCells
        integer                                 :: i, j, k, ib, jb, n, FlowDistribution
        real                                    :: SurfaceElevation, SurfaceElevationByPass    
        real                                    :: MaxFlow, DischargeFlow, AuxFlowIJ, FlowArea
        real                                    :: MinVolume
        integer                                 :: STAT_CALL
        logical                                 :: ByPassON
        integer, dimension(:    ), pointer      :: VectorI, VectorJ
        real,    dimension(:    ), pointer      :: DistributionCoef
        real                                    :: CoordinateX, CoordinateY, XBypass, YBypass      
        logical                                 :: CoordinatesON
!        real                                    :: ByPassFlowCriticCenterCell, FlowCriticCenterCell
        real                                    :: variation, variation2, DV, StabilizeFactor, Vnew, Hold
        real                                    :: AuxFlow
       
        !Begin------------------------------------------------------------------        
       
         if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ModifyWaterDischarges")


        !Sets to 0
        call SetMatrixValue(Me%lFlowDischarge, Me%Size, 0.0)
        
        !The discharge flow is controled using two basic rules:
        ! 1 - when the flow is negative can not remove more than the volume present in the cell;
        ! 2 - the volume variation induce by the discharge can not be larger than a percentage of the volume present in the cell.
        !     This percentage is equal to 100 * Me%CV%StabilizeFactor. By default Me%CV%StabilizeFactor = 0.1  this means that by 
        !     default this percentage is 1000 %. The Me%CV%StabilizeFactor is used for estimate changes in the time step to 
        !     maintain the model stability  
        !David -> Correcting user defined values in discharge should not be the default behaviour. 
        !In "normal" discharges this is now controlled by a keyword and default is false. And only used when stabilize is ON. 
        !Paulo suggestion is that in by pass discharges it should 
        
        StabilizeFactor = Me%CV%StabilizeFactor * 100.


        !Gets the number of discharges
        call GetDischargesNumber(Me%ObjDischarges, nDischarges, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR10'

        do iDis = 1, nDischarges
        
        
            if (Me%OutPut%TimeSerieDischON) then
                Me%OutPut%TimeSerieDischProp(iDis,:) = 0.
            endif    

            call GetDischargesGridLocalization(Me%ObjDischarges,                        &
                                               DischargeIDNumber = iDis,                &
                                               Igrid         = i,                       &
                                               JGrid         = j,                       &
                                               KGrid         = k,                       &
                                               IByPass       = ib,                      &
                                               JByPass       = jb,                      & 
                                               CoordinateX   = CoordinateX,             &
                                               CoordinateY   = CoordinateY,             & 
                                               CoordinatesON = CoordinatesON,           &
                                               XBypass       = XBypass,                 &
                                               YBypass       = YBypass,                 &
                                               STAT          = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR20'
            
            if (k == 0) then

                !Check if this is a bypass discharge. If it is gives the water level of the bypass end cell
                call GetByPassON(Me%ObjDischarges, iDis, ByPassON, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR30'
                
                if (CoordinatesON) then
                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordinateX, CoordinateY, I, J, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ModuleHydrodynamic - ERR40'

                    call CorrectsCellsDischarges(Me%ObjDischarges, iDis, I, J, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ModuleHydrodynamic - ERR50'                
                    
                    if (ByPassON) then

                        call GetXYCellZ(Me%ObjHorizontalGrid, XBypass, YBypass, Ib, Jb, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ModuleHydrodynamic - ERR60'

                        call CorrectsBypassCellsDischarges(Me%ObjDischarges, iDis, Ib, Jb, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ModuleHydrodynamic - ERR70'

                    endif
                    
                endif                

                if (ByPassON) then
                    SurfaceElevationByPass = Me%myWaterLevel(ib, jb)
                else
                    SurfaceElevationByPass = FillValueReal
                endif

                !real(8) to real as expected in GetDischargeWaterFlow
                SurfaceElevation = Me%myWaterLevel(i, j)
                call GetDischargeWaterFlow(Me%ObjDischarges,                            &
                                        Me%ExtVar%Now, iDis,                            &
                                        SurfaceElevation,                               &
                                        DischargeFlow,                                  &
                                        SurfaceElevation2 = SurfaceElevationByPass,     &
                                        FlowArea = FlowArea,                            &
                                        STAT = STAT_CALL)
                if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR80'

                
                call GetDischargeFlowDistribuiton(Me%ObjDischarges, iDis, nCells, FlowDistribution, &
                                                  VectorI, VectorJ, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR90'
                
                if (ByPassON) then
                    if (nCells > 1) then
                        stop 'ModuleRunOff - ModifyWaterDischarges - ERR100'
                    endif                            
                endif                

                !Horizontal distribution
i1:             if (nCells .ge. 1) then
                    allocate(DistributionCoef(1:nCells))
i2:                 if      (FlowDistribution == DischByCell_ ) then
                    
                        DistributionCoef(1:nCells) = 1./float(nCells)

                    else i2
                    
                        stop 'ModuleRunOff - ModifyWaterDischarges - ERR110'

                    endif i2
                endif i1
                
                AuxFlowIJ = DischargeFlow
                
 dn:            do n=1, nCells
 
                    if (nCells .ge. 1) then
                        i         = VectorI(n)
                        j         = VectorJ(n)
                        
                        !For every cell get the total flow and multiply it by distribution coef
                        call GetDischargeWaterFlow(Me%ObjDischarges,                            &
                                                Me%ExtVar%Now, iDis,                            &
                                                SurfaceElevation,                               &
                                                AuxFlowIJ,                                      &
                                                FlowDistribution  = DistributionCoef(n),        &
                                                STAT = STAT_CALL)
                        if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR120'


                    endif
                    
                    !each additional flow can remove all water column left
                    if (AuxFlowIJ < 0.0 .or. ByPassON) then
                    
                        if (ByPassON .and. AuxFlowIJ > 0.) then
                            !m3 = m * m2
                            MinVolume = Me%MinimumWaterColumn * Me%ExtVar%GridCellArea(ib, jb)
                            
                            !m3/s = m3 /s
                            if (Me%myWaterVolume(ib, jb) > MinVolume) then
                                MaxFlow = (Me%myWaterVolume(ib, jb) - MinVolume) / LocalDT
                            else
                                MaxFlow = 0.
                            endif                                
                                             
                            if (abs(AuxFlowIJ) > abs(MaxFlow)) then
                                if (AuxFlowIJ > 0.) then
                                    AuxFlowIJ =   MaxFlow
                                else
                                    AuxFlowIJ = - MaxFlow
                                endif                                    
                            endif    
                                
                            !m3/s      = [m/s^2*m]^0.5*[m^2]^0.5 * [m] = [m/s] * [m] * [m]    
                            !ByPassFlowCriticCenterCell = sqrt(Gravity * Me%myWaterColumn (ib, jb)) * sqrt(Me%ExtVar%GridCellArea(ib, jb)) * &
                            !                             Me%myWaterColumn (ib, jb)   
                            
                            !ByPassFlowCriticCenterCell = Me%CV%MaxCourant / 2. * ByPassFlowCriticCenterCell
                            
                            !if (abs(AuxFlowIJ) > abs(ByPassFlowCriticCenterCell)) then
                            !    AuxFlowIJ = - ByPassFlowCriticCenterCell 
                            !endif                                                                            
                       
                        else

                            !m3 = m * m2
                            MinVolume = Me%MinimumWaterColumn * Me%ExtVar%GridCellArea(i, j)

                            !m3/s = m3 /s
                            if (Me%myWaterVolume(i, j) > MinVolume) then
                                MaxFlow = (Me%myWaterVolume(i, j) - MinVolume) / LocalDT
                            else
                                MaxFlow = 0.
                            endif                                       
                            
                            if (abs(AuxFlowIJ) > abs(MaxFlow)) then
                                if (AuxFlowIJ > 0.) then
                                    AuxFlowIJ =   MaxFlow
                                else
                                    AuxFlowIJ = - MaxFlow
                                endif                                    
                            endif                              
                                    
                            !m3/s      = [m/s^2*m]^0.5*[m^2]^0.5 * [m] = [m/s] * [m] * [m]    
                            !FlowCriticCenterCell = sqrt(Gravity * Me%myWaterColumn (i, j)) * sqrt(Me%ExtVar%GridCellArea(i, j)) * Me%myWaterColumn (i, j)        
                            
                            !FlowCriticCenterCell = Me%CV%MaxCourant / 2. * FlowCriticCenterCell

                            !if (abs(AuxFlowIJ) > abs(FlowCriticCenterCell)) then
                            !!    AuxFlowIJ = - FlowCriticCenterCell 
                            !endif                            

                       endif                            
                       
                    endif
                    
                    !if not bypass, look to normal discharge option
                    !if bypass, look to bypass discharge option
                    if ((.not. ByPassON .and. Me%CV%CorrectDischarge) .or. (ByPassON .and. Me%CV%CorrectDischargeByPass)) then      
                        
                        Vnew = Me%myWaterVolume(i, j) + AuxFlowIJ * LocalDT
                        Hold = Me%myWaterVolumeOld(i, j) / Me%ExtVar%GridCellArea(i, j)
                    

                        if ((.not. Me%CV%CheckDecreaseOnly) .or. Me%myWaterVolumeOld(i, j) > Vnew) then
                    
                            if (Hold >= Me%CV%MinimumValueToStabilize) then
                    
                                DV =  Me%myWaterVolume(i, j)  - Me%myWaterVolumeOld(i, j)
                            
                                variation = abs(DV + AuxFlowIJ * LocalDT) / Me%myWaterVolumeOld(i, j)
                            
                                if (variation > StabilizeFactor) then
                                    AuxFlow = AuxFlowIJ
                                    variation2 = abs(DV) / Me%myWaterVolumeOld(i, j)                    
                                    if (variation2 > StabilizeFactor) then
                                        AuxFlowIJ = 0.
                                    else
                                        if (AuxFlowIJ > 0.) then
                                            AuxFlowIJ =  (  StabilizeFactor * Me%myWaterVolumeOld(i, j) - DV) / LocalDT
                                        else
                                            AuxFlowIJ =  (- StabilizeFactor * Me%myWaterVolumeOld(i, j) - DV) / LocalDT
                                        endif                                
                                    endif              
                                    write(*,*) 'Flow in cell',i,j,'was corrected from ',AuxFlow,'to ',AuxFlowIJ
                                endif
                            endif                        
                        endif

                    endif

                    !correct if bypass correction ON (is the default)
                    if (ByPassON .and. Me%CV%CorrectDischargeByPass) then
                        
                        Vnew = Me%myWaterVolume   (ib, jb) - AuxFlowIJ * LocalDT                    
                        Hold = Me%myWaterVolumeOld(ib, jb) / Me%ExtVar%GridCellArea(ib, jb)


                        if ((.not. Me%CV%CheckDecreaseOnly) .or. Me%myWaterVolumeOld(ib, jb) > Vnew) then
                        
                            if (Hold >= Me%CV%MinimumValueToStabilize) then

                                DV =  Me%myWaterVolume(ib, jb)  - Me%myWaterVolumeOld(ib, jb)
                                
                                variation = abs(DV - AuxFlowIJ * LocalDT) / Me%myWaterVolumeOld(ib, jb)
                                
                                if (variation > StabilizeFactor) then
                                    
                                    AuxFlow = AuxFlowIJ
                                    
                                    variation2 = abs(DV) / Me%myWaterVolumeOld(ib, jb)                    
                                    if (variation2 > StabilizeFactor) then
                                        AuxFlowIJ = 0.
                                    else
                                        if (AuxFlowIJ < 0.) then
                                                AuxFlowIJ =  (- StabilizeFactor * Me%myWaterVolumeOld(ib, jb) + DV) / LocalDT
                                        else
                                                AuxFlowIJ =  (  StabilizeFactor * Me%myWaterVolumeOld(ib, jb) + DV) / LocalDT
                                        endif                                
                                    endif  
                                    write(*,*) 'Flow in cell',i,j,'was corrected from ',AuxFlow,'to ',AuxFlowIJ
                                endif
                            endif 
                        endif                  

                    
                    endif

                    Me%lFlowDischarge(i, j)     = Me%lFlowDischarge(i, j) + AuxFlowIJ

                    !Updates Water Volume
                    Me%myWaterVolume(i, j)      = Me%myWaterVolume(i, j) + AuxFlowIJ * LocalDT

                    !Updates Water Column
                    Me%myWaterColumn  (i, j)    = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)      = Me%myWaterColumn (i, j) + Me%ExtVar%Topography(i, j)
                    
                    if (ByPassON) then

                        Me%lFlowDischarge(ib, jb)     = Me%lFlowDischarge(ib, jb) - AuxFlowIJ

                        !Updates Water Volume
                        Me%myWaterVolume(ib, jb)      = Me%myWaterVolume(ib, jb) - AuxFlowIJ * LocalDT

                        !Updates Water Column
                        Me%myWaterColumn  (ib, jb)    = Me%myWaterVolume (ib, jb) / Me%ExtVar%GridCellArea(ib, jb)

                        !Updates Water Level
                        Me%myWaterLevel (ib, jb)      = Me%myWaterColumn (ib, jb) + Me%ExtVar%Topography(ib, jb)
                    endif

                    !if (Me%CheckMass) Me%TotalInputVolume = Me%TotalInputVolume + Me%DischargesFlow(iDis) * LocalDT                    


                enddo dn

                if (nCells .ge. 1) deallocate(DistributionCoef)
                
                
                if (Me%OutPut%TimeSerieDischON) then
                    if (ByPassON) then
                        !In the output is assumed the flow direction Cell i,j (upstream) -> Cell Bypass i,j (downstream) as positive 
                        Me%OutPut%TimeSerieDischProp(iDis,1) = - AuxFlowIJ 
                    else
                        Me%OutPut%TimeSerieDischProp(iDis,1) =   AuxFlowIJ 
                    endif
                    if (FlowArea > 0.) then
                        Me%OutPut%TimeSerieDischProp(iDis,2) = Me%OutPut%TimeSerieDischProp(iDis,1) / FlowArea
                    else                            
                        Me%OutPut%TimeSerieDischProp(iDis,2) = FillValueReal
                    endif              

                    Me%OutPut%TimeSerieDischProp(iDis,3) = FlowArea                                      

                    Me%OutPut%TimeSerieDischProp(iDis,4) = SurfaceElevation

                    Me%OutPut%TimeSerieDischProp(iDis,5) = SurfaceElevationByPass
                    
                    
                    if (ByPassON) then
                        !In the output is assumed the flow direction Cell i,j (upstream) -> Cell Bypass i,j (downstream) as positive 
                        Me%OutPut%TimeSerieDischProp(iDis,6) = - DischargeFlow 
                    else
                        Me%OutPut%TimeSerieDischProp(iDis,6) =   DischargeFlow 
                    endif                    
                endif                  

                call UnGetDischarges(Me%ObjDischarges, VectorI, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModuleRunOff - ModifyWaterDischarges - ERR130'

                call UnGetDischarges(Me%ObjDischarges, VectorJ, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModuleRunOff - ModifyWaterDischarges - ERR140'                               
 
            endif
           
        enddo

         if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ModifyWaterDischarges")


    end subroutine ModifyWaterDischarges  
    
    !--------------------------------------------------------------------------

    subroutine ModifyGeometryAndMapping

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: WCL, WCR, WCA, Bottom
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ModifyGeometryAndMapping")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, WCL, WCR, WCA, Bottom)


        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j-1) == BasinPoint .and. Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                if (Me%FaceWaterColumn == WCMaxBottom_) then
                    !Maximum Bottom Level
                    Bottom = max(Me%ExtVar%Topography(i, j-1), Me%ExtVar%Topography(i, j))
                elseif (Me%FaceWaterColumn == WCAverageBottom_) then
                    !Average Bottom Level
                    Bottom = (Me%ExtVar%Topography(i,j) + Me%ExtVar%Topography(i,j-1)) / 2.0                
                endif
                
                !Water Column Left (above MaxBottom)
                WCL       = max(Me%myWaterLevel(i, j-1) - Bottom, dble(AlmostZero))
            
                !Water Column Right (above MaxBottom)
                WCR       = max(Me%myWaterLevel(i, j  ) - Bottom, dble(AlmostZero))

                !In the case of kinematic wave, always consider the "upstream" area, otherwise the average above "max bottom"
                if (Me%HydrodynamicApproximation == KinematicWave_) then
                    if (Me%ExtVar%Topography(i, j-1) > Me%ExtVar%Topography(i, j)) then
                        WCA = WCL
                    else
                        WCA = WCR
                    endif
                else
                    !Average Water Column
                    WCA = (WCL + WCR) / 2.0
                    if (Me%myWaterLevel(i, j-1) > Me%myWaterLevel(i, j) ) then
                        WCA = WCL
                    else
                        WCA = WCR
                    endif
                endif
                
                !Area  = Water Column * Side lenght of cell
                Me%AreaU(i, j) = WCA * Me%ExtVar%DYY(i, j)
                
                if (WCA > Me%MinimumWaterColumn) then
                    Me%ComputeFaceU(i, j) = 1
                else
                    Me%ComputeFaceU(i, j) = 0
                endif

            endif
            
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        !$OMP PARALLEL PRIVATE(I,J, WCL, WCR, WCA, Bottom)
        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%BasinPoints(i-1, j) == BasinPoint .and. Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                if (Me%FaceWaterColumn == WCMaxBottom_) then
                    !Maximum Bottom Level
                    Bottom = max(Me%ExtVar%Topography(i-1, j), Me%ExtVar%Topography(i, j))
                elseif (Me%FaceWaterColumn == WCAverageBottom_) then
                    !Average Bottom Level
                    Bottom = (Me%ExtVar%Topography(i,j) +  Me%ExtVar%Topography(i-1,j)) / 2.0                
                endif
                
                !Water Column Left
                WCL       = max(Me%myWaterLevel(i-1, j) - Bottom, dble(AlmostZero))
            
                !Water Column Right
                WCR       = max(Me%myWaterLevel(i, j  ) - Bottom, dble(AlmostZero))
               
                !In the case of kinematic wave, always consider the "upstream" area, otherwise the average above "max bottom"
                if (Me%HydrodynamicApproximation == KinematicWave_) then
                    if (Me%ExtVar%Topography(i-1, j) > Me%ExtVar%Topography(i, j)) then
                        WCA = WCL
                    else
                        WCA = WCR
                    endif
                else
                    !Average Water Column
                    !WCA = (WCL + WCR) / 2.0
                    if (Me%myWaterLevel(i-1, j) > Me%myWaterLevel(i, j)) then
                        WCA = WCL
                    else
                        WCA = WCR
                    endif
                endif

                !Area  = Water Column * Side lenght of cell
                Me%AreaV(i, j) = WCA * Me%ExtVar%DXX(i, j)
                
                if (WCA > Me%MinimumWaterColumn) then
                    Me%ComputeFaceV(i, j) = 1
                else
                    Me%ComputeFaceV(i, j) = 0
                endif

            endif
            
        enddo
        enddo    
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL
        
        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%myWaterColumn(i, j) .gt. Me%MinimumWaterColumn) then
                Me%OpenPoints(i,j) = 1
            else
                Me%OpenPoints(i,j) = 0
            endif
        enddo
        enddo    
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ModifyGeometryAndMapping")

    
    end subroutine ModifyGeometryAndMapping

    !--------------------------------------------------------------------------

    subroutine KinematicWave ()
    
        !Arguments-------------------------------------------------------------
        !real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Slope
        real                                        :: level_left, level_right
        real                                        :: level_bottom, level_top
        real                                        :: HydraulicRadius, WetPerimeter
        real                                        :: Margin1, Margin2
        real                                        :: WaterDepth, MaxBottom
        integer                                     :: CHUNK, di, dj
        real(8)                                     :: MaxFlow
        !character(len=StringLength)                 :: Direction

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, Slope, level_left, level_right, level_bottom, level_top, &
        !$OMP HydraulicRadius, MaxFlow, Margin1, Margin2, WaterDepth, MaxBottom, WetPerimeter, di, dj)
        
        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ComputeFaceU(i, j) == Compute) then
            
                !Adds to the final level the height of the buidings, if any
                if (Me%HydrodynamicApproximation == KinematicWave_) then
                
                    level_left  = Me%ExtVar%Topography(i, j-1)
                    level_right = Me%ExtVar%Topography(i, j)
                    
                elseif (Me%HydrodynamicApproximation == DiffusionWave_) then
                
                    level_left  = Me%myWaterLevel(i, j-1)
                    level_right = Me%myWaterLevel(i, j)

                else
                
                    write(*,*)'Internal error'
                
                endif
                    
                !Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_left - level_right) / Me%ExtVar%DZX(i, j-1))
                else
                    Slope           = (level_left - level_right) / Me%ExtVar%DZX(i, j-1)
                endif
                
                !Hydraulic Radius
!                Direction = "X"
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_left,level_right)
                !Wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DYY(i, j)
                
                !only compute in water column as MaxBottom (topography stairs descritization)
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !Water Depth consistent with AreaU computed (only water above max bottom)
                    WaterDepth = Me%AreaU(i,j) / Me%ExtVar%DYY(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))
                    
                    !to check wich cell to use to use since areaU depends on higher water level and max bottom
                    if (level_left .gt. level_right) then
                        dj = -1
                    else
                        dj = 0
                    endif
                   
                    !Bottom Difference to adjacent cells (to check existence of “margins” on the side)
                    Margin1 = Me%ExtVar%Topography(i+1, j + dj) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i-1, j + dj) - MaxBottom

                    !if positive than there is a “margin” on the side and friction occurs at wet length
                    !If not basin points than result will be negative.
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif
                
                HydraulicRadius = Me%AreaU(i, j) / WetPerimeter
                             
                
                !
                !MANNING'S EQUATION -  KINEMATIC WAVE
                !
                !m3.s-1 = m2 * m(2/3) / (s.m(-1/3)) = m(8/3) * m(1/3) / s = m3.s-1
                if (Slope >= 0.0) then
                    Me%lFlowX(i, j) = Me%AreaU(i, j) * HydraulicRadius**(2./3.) * sqrt(Slope)          &
                                      / Me%OverlandCoefficientX(i,j)
                else
                    Me%lFlowX(i, j) = - Me%AreaU(i, j) * HydraulicRadius**(2./3.) * sqrt(-1.0 * Slope) &
                                      / Me%OverlandCoefficientX(i,j)
                endif
                
                
                !Limits Velocity to celerity if a free drop exists
                if (Me%HydrodynamicApproximation == DiffusionWave_ .and. Me%LimitToCriticalFlow) then
                    if ((level_left .lt. Me%ExtVar%Topography(i,j)) .or. (level_right .lt. Me%ExtVar%Topography(i,j-1))) then
                        
                        !already defined in shorter
                        !WaterDepth = max (level_left, level_right) - max(Me%ExtVar%Topography(i, j-1), Me%ExtVar%Topography(i, j))
                        WaterDepth      = Me%AreaU(i, j)/Me%ExtVar%DYY(i,j)
                        MaxFlow         = Me%AreaU(i, j) * sqrt(Gravity * WaterDepth)
                        Me%lFlowX(i, j) = Min (MaxFlow, Me%lFlowX(i, j))       
                                    
                    endif
    
                endif
                
            else
                
                Me%lFlowX(i, j) = 0.0
            
            endif
                
        enddo
        enddo        
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ComputeFaceV(i, j) == Compute) then
            
                !Adds to the final level the height of the buidings, if any
                if (Me%HydrodynamicApproximation == KinematicWave_) then
                
                    level_bottom = Me%ExtVar%Topography(i-1, j)
                    level_top    = Me%ExtVar%Topography(i, j)
                    
                else if (Me%HydrodynamicApproximation == DiffusionWave_) then

                    level_bottom = Me%myWaterLevel(i-1, j)
                    level_top    = Me%myWaterLevel(i, j)

                else
                
                    write(*,*)'Internal error'
                
                endif

                
                !Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_bottom - level_top) / Me%ExtVar%DZY(i-1, j))
                else
                    Slope           = (level_bottom - level_top) / Me%ExtVar%DZY(i-1, j)
                endif
                
                !Hydraulic Radius
!                Direction = "Y"
!               !This function produced an overhead in openmp and the simulation took
!               !double the time so it was abandoned
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_bottom,level_top)                
                !Wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DXX(i, j)
                
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !Water Depth consistent with AreaV computed (only water above max bottom)
                    WaterDepth = Me%AreaV(i,j) / Me%ExtVar%DXX(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))

                    !to check wich cell to use since areaV depends on higher water level
                    if (level_bottom .gt. level_top) then
                        di = -1
                    else
                        di = 0
                    endif

                    !Bottom Difference to adjacent cells (to check existence of “margins” on the side)
                    Margin1 = Me%ExtVar%Topography(i + di,j+1) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i + di,j-1) - MaxBottom

                    !if positive than there is a “margin” on the side and friction occurs at wet length
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif
                
                !m = m2 / m
                HydraulicRadius = Me%AreaV(i, j) / WetPerimeter
                
                                
                !
                !MANNING'S EQUATION -  KINEMATIC WAVE
                !
                !m3.s-1 = m2 * m(2/3) / (s.m(-1/3)) = m(8/3) * m(1/3) / s = m3.s-1
                if (Slope >= 0.0) then
                    Me%lFlowY(i, j) = Me%AreaV(i, j) * HydraulicRadius**(2./3.) * sqrt(Slope)            &
                                      / Me%OverlandCoefficientY(i,j)
                else
                    Me%lFlowY(i, j) = - Me%AreaV(i, j) * HydraulicRadius**(2./3.) * sqrt(-1.0 * Slope)   &
                                      / Me%OverlandCoefficientY(i,j)
                endif
                
                !Limits Velocity to reasonable values
                if (Me%HydrodynamicApproximation == DiffusionWave_ .and. Me%LimitToCriticalFlow) then

                    if ((level_bottom .lt. Me%ExtVar%Topography(i,j)) .or. (level_top .lt. Me%ExtVar%Topography(i-1,j))) then
                        
                        !already defined in shorter
                        !WaterDepth = max (level_bottom, level_top) - max(Me%ExtVar%Topography(i-1,j), Me%ExtVar%Topography(i, j))
                        WaterDepth      = Me%AreaV(i, j)/Me%ExtVar%DXX(i,j)
                        MaxFlow         = Me%AreaV(i, j) * sqrt(Gravity * WaterDepth)
                        Me%lFlowY(i, j) = Min (MaxFlow, Me%lFlowY(i, j))
                    
                    endif

                
                
                endif
                
            else
            
                Me%lFlowY(i, j) = 0.0
            
            endif

        enddo
        enddo              
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    
    end subroutine KinematicWave
    
    !--------------------------------------------------------------------------
    
    subroutine ComputeFaceVelocityModulus
    
    !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB    
        integer                                             :: i, j
        real                                                :: U, V, Uaverage, Vaverage
        integer                                             :: CHUNK
        
        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ComputeFaceVelocityModulus")

        
        !$OMP PARALLEL PRIVATE(I,J, U, Vaverage)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNKJ)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ComputeFaceU(i, j) == Compute) then
                
                Vaverage = (Me%FlowYOld(i,  j  )/Me%AreaV(i,  j  ) + &
                            Me%FlowYOld(i+1,j  )/Me%AreaV(i+1,j  ) + &
                            Me%FlowYOld(i+1,j-1)/Me%AreaV(i+1,j-1) + &
                            Me%FlowYOld(i,  j-1)/Me%AreaV(i,  j-1)) /4.0
                
                U = Me%FlowXOld(i,j)/Me%AreaU(i,j)
                
                !Me%VelModFaceU(i, j) = sqrt(U**2.0 + Vaverage**2.0)
                Me%VelModFaceU(i, j) = abs(cmplx(U, Vaverage))
                
            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        !$OMP PARALLEL PRIVATE(I,J, V, Uaverage)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNKJ)
        do j = JLB, JUB 
        do i = ILB, IUB

            if (Me%ComputeFaceV(i, j) == Compute) then
                
                Uaverage = (Me%FlowXOld(i   ,j)/Me%AreaU(i  ,j   ) + &
                            Me%FlowXOld(i-1,j )/Me%AreaU(i-1,j   ) + &
                            Me%FlowXOld(i-1,j+1)/Me%AreaU(i-1,j+1) + &
                            Me%FlowXOld(i  ,j+1)/Me%AreaU(i  ,j+1)) /4.0
                
                V = Me%FlowYOld(i,j)/Me%AreaV(i,j)
                
                !Me%VelModFaceV(i, j) = sqrt(V**2.0 + Uaverage**2.0)
                Me%VelModFaceV(i, j) = abs(cmplx(Uaverage, V))
                
            endif

        enddo
        enddo
        
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    
        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ComputeFaceVelocityModulus")

    end subroutine ComputeFaceVelocityModulus
    
    !-------------------------------------------------------------------------
    
    subroutine DynamicWaveXX (LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Slope
        real                                        :: level_left, level_right
        real                                        :: HydraulicRadius
        real                                        :: Friction
        real                                        :: Pressure
        !real                                        :: upAdv, downAdv, 
        real                                        :: XLeftAdv, XRightAdv, YBottomAdv, YTopAdv
        real                                        :: Advection, Qf, WetPerimeter
        real(8)                                     :: CriticalFlow
        real                                        :: Margin1, Margin2
        integer                                     :: CHUNK, dj
        real                                        :: MaxBottom, WaterDepth
        !character(len=StringLength)                 :: Direction


        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "DynamicWaveXX")


        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, Slope, level_left, level_right, &
        !$OMP HydraulicRadius, Friction, Pressure, XLeftAdv, XRightAdv, YBottomAdv, YTopAdv, Advection, Qf, &
        !$OMP CriticalFlow, Margin1, Margin2, MaxBottom, WaterDepth, dj, WetPerimeter)

        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNKJ)
        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%ComputeFaceU(i, j) == Compute) then
            
                level_left  = Me%myWaterLevel(i, j-1)
                level_right = Me%myWaterLevel(i, j)
                    
                !!Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_left - level_right) / Me%ExtVar%DZX(i, j-1))
                else
                    Slope           = (level_left - level_right) / Me%ExtVar%DZX(i, j-1)
                endif
                 
                !!Hydraulic Radius
!                Direction = "X"
!                !This function produced an overhead in openmp and the simulation took 
!                !double the time so it was abandoned
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_left,level_right)
                !wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DYY(i, j)
                
                !only compute margins if water column method is MaxBottom (topography discretization by "stairs")
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !Then, is checked if "margins" occur on the cell of the highest water level
                    !water depth consistent with AreaU computed (only water above max bottom)
                    WaterDepth = Me%AreaU(i,j) / Me%ExtVar%DYY(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))
                    
                    !to check which cell to use since areaU depends on higher water level
                    if (level_left .gt. level_right) then
                        dj = -1
                    else
                        dj = 0
                    endif
                    
                    !bottom Difference to adjacent cells (to check existence of “margins” on the side)
                    Margin1 = Me%ExtVar%Topography(i+1, j + dj) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i-1, j + dj) - MaxBottom

                    !if positive, than there is a “margin” on the side and friction occurs at wet length
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif
                
                HydraulicRadius = Me%AreaU(i, j) / WetPerimeter
       
                !
                !Sant Venant
                !

                !Pressure
                !m3/s             = s  * m/s2    * m2   * m/m
                Pressure          = LocalDT * Gravity * Me%AreaU(i, j) * Slope



                !FRICTION - semi-implicit -----------------------------------------------
                !   -    =  (s * m.s-2  * m3.s-1 * (s.m(-1/3))^2) / (m2 * m(4/3)) = m(10/3) / m(10/3)
                !Friction = LocalDT * Gravity * abs(Me%FlowXOld(i, j)) * Me%OverlandCoefficientX(i,j)** 2. &
                !     / ( Me%AreaU(i, j) * HydraulicRadius ** (4./3.) ) 
                
                !Friction = LocalDT * Gravity * &
                !           sqrt(Me%FlowXOld(i, j)**2. + Me%FlowYOld(i, j)**2.) * Me%OverlandCoefficientX(i,j)** 2. / &
                !           ( Me%AreaU(i, j) * HydraulicRadius ** (4./3.) ) 
                
                Friction = LocalDT * Gravity * &
                           Me%VelModFaceU(i,j) * Me%OverlandCoefficientX(i,j)** 2. / &
                           (HydraulicRadius ** (4./3.)) 
                
                !Advection (may be limited to water column height)
                if ((Me%ComputeAdvectionU(i,j) == 1) .and. (Me%myWaterColumn(i,j) .gt. Me%MinimumWaterColumnAdvection)  .and.   &
                    (Me%myWaterColumn(i,j-1) .gt. Me%MinimumWaterColumnAdvection)) then
                    
                    
                    !Face XU(i,j+1). Z U Faces have to be open
                    if ((Me%ComputeFaceU(i, j) +  Me%ComputeFaceU(i, j+1) == 2)) then 

                        !OLD Version
                        !Theold formulation had a problem when flows in adjacent reaches
                        !had opposite directions. Flow was the average and velocity would be
                        !in opposite  direction of average flow.
                       
                        !New Version 
                        !The new formulation, in case of opposite directions, in adjacent reaches does not compute
                        !advection. In case of same direction, is hard-upwind meaning that it will use flow and 
                        !velocity from the upwind reach. This option may be more stable than soft-upwind 
                        !(average flow and velocity from upwind reach) or central differences (average flow 
                        !and average velocity).
                        !if flows in same direction, advection is computed                        
                        if ((Me%FlowXOld(i, j) * Me%FlowXOld(i, j+1)).ge. 0.0) then
                            
                            Qf = (Me%FlowXOld(i, j) + Me%FlowXOld(i, j+1)) / 2.0

                            if (Qf > 0.0) then
                                XRightAdv = Me%FlowXOld(i, j)   * Me%FlowXOld(i, j) / Me%AreaU(i, j)
                            else
                                XRightAdv = Me%FlowXOld(i, j+1) * Me%FlowXOld(i, j+1) / Me%AreaU(i, j+1)
                            endif
                        else
                            XRightAdv = 0.0
                        endif
                        
                    else
                        XRightAdv = 0.0
                    endif       
                    
                    !Face XU(i,j). Z U Faces have to be open
                    if ((Me%ComputeFaceU(i, j-1) + Me%ComputeFaceU(i, j) == 2)) then  
                        
                        !New Version
                        if ((Me%FlowXOld(i, j-1) * Me%FlowXOld(i, j)) .ge. 0.0) then
                            
                            Qf = (Me%FlowXOld(i, j-1) + Me%FlowXOld(i, j)) / 2.0

                            if (Qf > 0.0) then
                                XLeftAdv = Me%FlowXOld(i, j-1) * Me%FlowXOld(i, j-1) / Me%AreaU(i, j-1)
                            else
                                XLeftAdv = Me%FlowXOld(i, j) * Me%FlowXOld(i, j) / Me%AreaU(i, j)
                            endif
                        else
                            XLeftAdv = 0.0
                        endif
                        
                        
                    else
                        XLeftAdv = 0.0
                    endif       
                    
                    if (Me%ComputeFaceV(i+1, j-1) +  Me%ComputeFaceV(i+1, j) .gt. 0) then

                        !if flows in same direction, advection is computed                        
                        if ((Me%FlowYOld(i+1, j-1) * Me%FlowYOld(i+1, j)).ge. 0.0) then
                            
                            Qf = (Me%FlowYOld(i+1, j-1) + Me%FlowYOld(i+1, j)) / 2.0
                            
                            if (Qf > 0.0) then
                                YTopAdv = Qf   * Me%FlowXOld(i, j) / Me%AreaU(i, j)
                            elseif (Qf < 0.0) then
                                if(Me%ComputeFaceU(i+1,j) == Compute) then
                                    YTopAdv = Qf * Me%FlowXOld(i+1, j) / Me%AreaU(i+1, j)
                                    !YTopAdv =  (Me%FlowYOld(i+1, j-1) * (Me%FlowXOld(i+1, j-1) / Me%AreaU(i+1, j-1)        + &
                                    !                                    Me%FlowXOld(i+1,   j) / Me%AreaU(i+1,   j)) / 2.0 + & 
                                    !           Me%FlowYOld(i+1,   j) * (Me%FlowXOld(i+1,   j) / Me%AreaU(i+1,   j)        + &
                                    !                                    Me%FlowXOld(i+1, j+1) / Me%AreaU(i+1, j+1)) / 2.0) / 2.0
                                else
                                    if(Me%ComputeFaceU(i+1, j-1)== Compute)then
                                        YTopAdv = Qf * Me%FlowXOld(i+1, j-1) / Me%AreaU(i+1, j-1) 
                                    elseif(Me%ComputeFaceU(i+1, j+1) == Compute)then
                                        YTopAdv = Qf * Me%FlowXOld(i+1, j+1) / Me%AreaU(i+1, j+1)
                                    else
                                        YTopAdv = 0.0
                                    endif
                                endif
                            else
                                YTopAdv = 0.0
                            endif
                        else
                            YTopAdv = 0.0
                        endif
                        
                    else
                        YTopAdv = 0.0
                    endif       
                    

                    if (Me%ComputeFaceV(i, j-1) +  Me%ComputeFaceV(i, j) .gt. 0) then

                        !if flows in same direction, advection is computed                        
                        if ((Me%FlowYOld(i, j-1) * Me%FlowYOld(i, j)).ge. 0.0) then
                            
                            Qf = (Me%FlowYOld(i, j-1) + Me%FlowYOld(i, j)) / 2.0
                            
                            if (Qf > 0.0)then
                                if(Me%ComputeFaceU(i-1,j) == Compute) then
                                    YBottomAdv =  Qf   * Me%FlowXOld(i-1, j) / Me%AreaU(i-1, j)
                                    
                                    !YBottomAdv =  (Me%FlowYOld(i, j-1) * (Me%FlowXOld(i-1, j-1) / Me%AreaU(i-1, j-1)       + &
                                    !                                     Me%FlowXOld(i-1,   j) / Me%AreaU(i-1,   j)) / 2.0 + & 
                                    !              Me%FlowYOld(i,   j) * (Me%FlowXOld(i-1,   j) / Me%AreaU(i-1,   j)        + &
                                    !                                     Me%FlowXOld(i-1, j+1) / Me%AreaU(i-1, j+1)) / 2.0) / 2.0
                                    !
                                    !YBottomAdv =  (Me%FlowYOld(i, j-1) * (Me%FlowXOld(i-1, j-1) / Me%AreaU(i-1, j-1)       + &
                                    !                                     Me%FlowXOld(i-1,   j) / Me%AreaU(i-1,   j)) / 2.0 + & 
                                    !              Me%FlowYOld(i,   j) * (Me%FlowXOld(i-1,   j) / Me%AreaU(i-1,   j)        + &
                                    !                                     Me%FlowXOld(i-1, j+1) / Me%AreaU(i-1, j+1)) / 2.0) / 2.0
                                    
                                    
                                    
                                else
                                    if(Me%ComputeFaceU(i-1, j-1) == Compute)then
                                        YBottomAdv =  Qf * Me%FlowXOld(i-1, j-1) / Me%AreaU(i-1, j-1) 
                                    elseif(Me%ComputeFaceU(i-1, j+1) == Compute)then
                                        YBottomAdv =  Qf * Me%FlowXOld(i-1, j+1) / Me%AreaU(i-1, j+1)
                                    else
                                        YBottomAdv =  0.0
                                    endif
                                    
                                endif
                            elseif ((Qf < 0.0)) then
                                YBottomAdv = Qf   * Me%FlowXOld(i, j) / Me%AreaU(i, j)
                            else
                                YBottomAdv = 0.0
                            endif
                        else
                            YBottomAdv = 0.0
                        endif
                        
                    else
                        YBottomAdv = 0.0
                    endif       
                    
                    
                    Advection = (XLeftAdv - XRightAdv) * LocalDT / Me%ExtVar%DZX(i, j-1) + &
                                (YBottomAdv - YTopAdv) * LocalDT / Me%ExtVar%DYY(i, j)
                                
                else
                
                    Advection = 0.0
                    
                endif

                Me%lFlowX(i, j) = (Me%FlowXOld(i, j) + Pressure + Advection) / (1.0 + Friction)

                if (Me%LimitToCriticalFlow) then

                    !Limit to critical flow. Using the critical flow limitation in all cells assumes "slow" flow or
                    !subcritical that is consistent with the formulation used (flow depends on downstream height)
                    !because in supercritical flow it is only dependent on upstream and descritization to describe it would have
                    !to change. Supercritical flow usually exists on hydraulic infraestructures (high drops) and a 
                    !hydraulic jump exists between fast flow and slow flow.
                    
                    !Test Limitation only if free drop exists
!                    if ((level_left .lt. Me%ExtVar%Topography(i,j)) .or. (level_right .lt. Me%ExtVar%Topography(i,j-1))) then

                        !Waterdepth at the center of the face - depending on flow direction since flow
                        !can be in opposite direction of height gradient (AreaU uses the higher water level)              
                        !WaterDepth = Me%AreaU(i,j)/Me%ExtVar%DYY(i,j)
                        if (Me%FaceWaterColumn == WCMaxBottom_) then
                            MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))     
                                                            
                            if (Me%lFlowX(i, j) .gt. 0.0) then           
                                WaterDepth = max(Me%MyWaterLevel(i,j-1) - MaxBottom, 0.0)
                            else
                                WaterDepth = max(Me%MyWaterLevel(i,j) - MaxBottom, 0.0)
                            endif
                        elseif (Me%FaceWaterColumn == WCAverageBottom_) then
                            if (Me%lFlowX(i, j) .gt. 0.0) then           
                                WaterDepth = Me%MyWaterColumn(i,j-1)
                            else
                                WaterDepth = Me%MyWaterColumn(i,j)
                            endif                        
                        endif
                        
                        !Critical Flow
                        !CriticalFlow = Me%AreaU(i, j) * sqrt(Gravity * WaterDepth)
                        !m3/s = m * m * m/s
                        CriticalFlow = WaterDepth * Me%ExtVar%DYY(i,j) * sqrt(Gravity * WaterDepth)
                        
                        !only limit if flow higher
                        if (abs(Me%lFlowX(i, j)) > CriticalFlow) then
                            if (Me%lFlowX(i, j) > 0) then
                                Me%lFlowX(i, j) = CriticalFlow
                            else
                                Me%lFlowX(i, j) = -1.0 * CriticalFlow
                            endif
                        endif
 !                   endif
                
                else
                    !Predict water column to avoid negative volumes since 4 fluxes exist and the sum may be more than exists
                    if (Me%lFlowX(i, j) .lt. 0.0) then
                        if (abs(Me%lFlowX(i, j))* LocalDT .gt. Me%myWaterVolumePred(i,j)) then
                            Me%lFlowX(i, j) = - Me%myWaterVolumePred(i,j) / LocalDT
                        endif
                    elseif (Me%lFlowX(i, j) .gt. 0.0) then
                        if (Me%lFlowX(i, j)* LocalDT .gt. Me%myWaterVolumePred(i,j-1)) then
                            Me%lFlowX(i, j) =  Me%myWaterVolumePred(i,j-1) / LocalDT
                        endif
                    endif 
                    
                    !m3 = m3 + (-m3/s * s)
                    Me%myWaterVolumePred(i,j  ) = Me%myWaterVolumePred(i,j  ) + (Me%lFlowX(i, j) * LocalDT)
                    Me%myWaterVolumePred(i,j-1) = Me%myWaterVolumePred(i,j-1) - (Me%lFlowX(i, j) * LocalDT) 
                                   
                endif
                
            else
            
                Me%lFlowX(i, j) = 0.0

            endif

        enddo
        enddo        
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "DynamicWaveXX")
        
        
    end subroutine DynamicWaveXX
    
    !--------------------------------------------------------------------------

    subroutine DynamicWaveYY (LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Slope
        real                                        :: level_bottom, level_top
        real                                        :: HydraulicRadius
        real                                        :: Friction
        real                                        :: Pressure
        !real                                        :: upAdv, downAdv, 
        real                                        :: XLeftAdv, XRightAdv, YBottomAdv, YTopAdv
        real                                        :: Advection, Qf, WetPerimeter
        real(8)                                     :: CriticalFlow
        real                                        :: Margin1, Margin2
        integer                                     :: CHUNK, di
        real                                        :: MaxBottom, WaterDepth
        !character(len=StringLength)                 :: Direction

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "DynamicWaveYY")


        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, Slope, level_bottom, level_top, &
        !$OMP HydraulicRadius, Friction, Pressure, XLeftAdv, XRightAdv, YBottomAdv, YTopAdv, Advection, Qf, &
        !$OMP CriticalFlow, Margin1, Margin2, MaxBottom, WaterDepth, di, WetPerimeter)

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNKJ)
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ComputeFaceV(i, j) == Compute) then
            
                level_bottom = Me%myWaterLevel(i-1, j)
                level_top    = Me%myWaterLevel(i, j)
                
                !!Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_bottom - level_top) / Me%ExtVar%DZY(i-1, j))
                else
                    Slope           = (level_bottom - level_top) / Me%ExtVar%DZY(i-1, j)
                endif
                
                !!Hydraulic Radius
!                Direction = "Y"
!                !This function produced an overhead with openmp so it was abandoned
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_bottom,level_top)
                
                !wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DXX(i, j)
                
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !water Depth consistent with AreaV computed (only water above max bottom)
                    WaterDepth = Me%AreaV(i,j) / Me%ExtVar%DXX(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))

                    !to check wich cell to use since areaV depends on higher water level
                    if (level_bottom .gt. level_top) then
                        di = -1
                    else
                        di = 0
                    endif

                    !bottom Difference to adjacent cells (to check existence of “margins” on the side)
                    Margin1 = Me%ExtVar%Topography(i + di,j+1) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i + di,j-1) - MaxBottom

                    !if positive than there is a “margin” on the side and friction occurs at wet length
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif
                
                !m = m2 / m
                HydraulicRadius = Me%AreaV(i, j) / WetPerimeter
               
                !
                !Sant Venant
                !

                !m3/s             = s  * m/s2    * m2   * m/m
                Pressure          = LocalDT * Gravity * Me%AreaV(i, j) * Slope

                !FRICTION - semi-implicit -----------------------------------------------
                !   -    =  (s * m.s-2  * m3.s-1 * (s.m(-1/3))^2) / (m2 * m(4/3)) = m(10/3) / m(10/3)
                !Friction = LocalDT * Gravity * abs(Me%FlowYOld(i, j)) * Me%OverlandCoefficientY(i,j) ** 2. &
                !         / ( Me%AreaV(i, j) * HydraulicRadius ** (4./3.) ) 
                
                !Friction = LocalDT * Gravity * &
                !           sqrt(Me%FlowXOld(i, j)**2. + Me%FlowYOld(i, j)**2.) * Me%OverlandCoefficientY(i,j)** 2. / &
                !           ( Me%AreaV(i, j) * HydraulicRadius ** (4./3.) ) 
                
                Friction = LocalDT * Gravity * &
                           Me%VelModFaceV(i,j) * Me%OverlandCoefficientY(i,j)** 2. / &
                           (HydraulicRadius ** (4./3.)) 
        

                !Advection
                if ((Me%ComputeAdvectionV(i,j) == 1) .and. (Me%myWaterColumn(i,j) .gt. Me%MinimumWaterColumnAdvection)  &
                     .and. (Me%myWaterColumn(i-1,j) .gt. Me%MinimumWaterColumnAdvection)) then
                    
                    !Face YV(i+1,j)
                    if ((Me%ComputeFaceV(i, j) +  Me%ComputeFaceV(i+1, j) == 2)) then 
                        
                        if ((Me%FlowYOld(i, j) * Me%FlowYOld(i+1, j)) .ge. 0.0) then
                            
                            Qf = (Me%FlowYOld(i, j) + Me%FlowYOld(i+1, j)) / 2.0

                            if (Qf > 0.0) then
                                YTopAdv = Me%FlowYOld(i, j)   * Me%FlowYOld(i, j) / Me%AreaV(i, j)
                            else
                                YTopAdv = Me%FlowYOld(i+1, j) * Me%FlowYOld(i+1, j) / Me%AreaV(i+1, j)
                            endif
                        else
                            YTopAdv = 0.0
                        endif
                        
                    else
                        YTopAdv = 0.0
                    endif
                    
                    !Face YV(i,j)
                    if ((Me%ComputeFaceV(i-1, j) + Me%ComputeFaceV(i, j) == 2)) then 

                        if ((Me%FlowYOld(i-1, j) * Me%FlowYOld(i, j)) .ge. 0.0) then
                            
                            Qf = (Me%FlowYOld(i-1, j) + Me%FlowYOld(i, j)) / 2.0

                            if (Qf > 0.0) then
                                YBottomAdv = Me%FlowYOld(i-1, j)   * Me%FlowYOld(i-1, j) / Me%AreaV(i-1, j)
                            else
                                YBottomAdv = Me%FlowYOld(i, j) * Me%FlowYOld(i, j) / Me%AreaV(i, j)
                            endif
                        else
                            YBottomAdv = 0.0
                        endif
                        
                    else
                        YBottomAdv = 0.0
                    endif                
        
                    if (Me%ComputeFaceU(i, j+1) +  Me%ComputeFaceU(i-1, j+1) .gt. 0) then
                        
                        !if flows in same direction, advection is computed                        
                        if ((Me%FlowXOld(i, j+1) * Me%FlowXOld(i-1, j+1)).ge. 0.0) then
                            
                            Qf = (Me%FlowXOld(i, j+1) + Me%FlowXOld(i-1, j+1)) / 2.0
                            
                            if (Qf > 0.0) then
                                
                                XRightAdv = Qf   * Me%FlowYOld(i, j) / Me%AreaV(i, j)
                                
                            elseif (Qf < 0.0)then
                                
                                if(Me%ComputeFaceV(i,j+1) == Compute) then
                                    XRightAdv = Qf   * Me%FlowYOld(i, j+1) / Me%AreaV(i, j+1)
                                    !XRightAdv =  (Me%FlowXOld(i,   j) * (Me%FlowYOld(i+1, j-1) / Me%AreaV(i+1, j-1)        + &
                                    !                                    Me%FlowYOld(i,   j-1) / Me%AreaV(i,   j-1)) / 2.0 + & 
                                    !             Me%FlowXOld(i-1, j) * (Me%FlowYOld(i,   j-1) / Me%AreaV(i,   j-1)        + &
                                    !                                    Me%FlowYOld(i-1, j-1) / Me%AreaV(i-1, j-1)) / 2.0)/2.0
                                else
                                    if(Me%ComputeFaceV(i-1, j+1) == Compute)then
                                        XRightAdv = Qf   *  Me%FlowYOld(i-1, j+1) / Me%AreaV(i-1, j+1)
                                    elseif(Me%ComputeFaceV(i+1, j+1) == Compute)then
                                        XRightAdv = Qf   * Me%FlowYOld(i+1, j+1) / Me%AreaV(i+1, j+1) 
                                    else
                                        XRightAdv = 0.0
                                    endif
                                endif
                                
                            else 
                                XRightAdv = 0.0
                            endif
                        else
                            XRightAdv = 0.0
                        endif
                        
                    else
                        XRightAdv = 0.0
                    endif       

                    if (Me%ComputeFaceU(i, j) +  Me%ComputeFaceU(i-1, j) .gt. 0) then
                        
                        !if flows in same direction, advection is computed                        
                        if ((Me%FlowXOld(i, j) * Me%FlowXOld(i-1, j)).ge. 0.0) then
                            
                            Qf = (Me%FlowXOld(i, j) + Me%FlowXOld(i-1, j)) / 2.0
                            
                            if (Qf > 0.0)then
                                if(Me%ComputeFaceV(i,j-1) == Compute) then
                                    XLeftAdv = Qf   * Me%FlowYOld(i, j-1) / Me%AreaV(i, j-1)
                                    
                                    !XLeftAdv =  (Me%FlowXOld(i,   j) * (Me%FlowYOld(i+1, j-1) / Me%AreaV(i+1, j-1)        + &
                                    !                                   Me%FlowYOld(i,   j-1) / Me%AreaV(i,   j-1)) / 2.0 + & 
                                    !            Me%FlowXOld(i-1, j) * (Me%FlowYOld(i,   j-1) / Me%AreaV(i,   j-1)        + &
                                    !                                   Me%FlowYOld(i-1, j-1) / Me%AreaV(i-1, j-1)) / 2.0)/2.0
                                    
                                else
                                    
                                    if(Me%ComputeFaceV(i+1, j-1) == Compute)then
                                        XLeftAdv = Qf * Me%FlowYOld(i+1, j-1) / Me%AreaV(i+1, j-1) 
                                    elseif(Me%ComputeFaceV(i-1, j-1) == Compute)then
                                        XLeftAdv = Qf * Me%FlowYOld(i-1, j-1) / Me%AreaV(i-1, j-1)
                                    else
                                        XLeftAdv = 0.0
                                    endif
                                   
                                endif
                            elseif (Qf < 0.0) then
                                XLeftAdv = Qf   * Me%FlowYOld(i, j) / Me%AreaV(i, j)
                            else
                                XLeftAdv = 0.0
                            endif
                        else
                            XLeftAdv = 0.0
                        endif
                        
                    else
                        XLeftAdv = 0.0
                    endif       
                           
                    !Advection = (upAdv - downAdv) * LocalDT / Me%ExtVar%DVY(i, j)
                    Advection = (YBottomAdv - YTopAdv) * LocalDT / Me%ExtVar%DZY(i-1, j)     &
                                + (XLeftAdv - XRightAdv) * LocalDT / Me%ExtVar%DXX(i, j)
                    
                else
                
                    Advection = 0.0
                    
                endif
                
                Me%lFlowY(i, j) = (Me%FlowYOld(i, j) + Pressure + Advection) / (1.0 + Friction)
                
                
                if (Me%LimitToCriticalFlow) then
                    
!                    if ((level_bottom .lt. Me%ExtVar%Topography(i,j)) .or. (level_top .lt. Me%ExtVar%Topography(i-1,j))) then
                    
                        !Waterdepth at the center of the face - depending on flow direction since flow
                        !can be in opposite direction of height gradient (AreaU uses the higher)
                        !WaterDepth = Me%AreaV(i,j)/Me%ExtVar%DXX(i,j)
                        if (Me%FaceWaterColumn == WCMaxBottom_) then
                            MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))
                                                            
                            if (Me%lFlowY(i, j) .gt. 0.0) then           
                                WaterDepth = max(Me%MyWaterLevel(i-1,j) - MaxBottom, 0.0)
                            else
                                WaterDepth = max(Me%MyWaterLevel(i,j) - MaxBottom, 0.0)
                            endif                
                        elseif (Me%FaceWaterColumn == WCAverageBottom_) then
                            if (Me%lFlowY(i, j) .gt. 0.0) then           
                                WaterDepth = Me%MyWaterColumn(i-1,j)
                            else
                                WaterDepth = Me%MyWaterColumn(i,j)
                            endif                                        
                        endif
                        
                        !Critical Flow
                        !CriticalFlow = Me%AreaV(i, j) * sqrt(Gravity * WaterDepth)
                        !m3/s = m * m * m/s
                        CriticalFlow = WaterDepth * Me%ExtVar%DXX(i,j) * sqrt(Gravity * WaterDepth)
                        
                        !only limit if flow higher
                        if (abs(Me%lFlowY(i, j)) > CriticalFlow) then
                            if (Me%lFlowY(i, j) > 0) then
                                Me%lFlowY(i, j) = CriticalFlow
                            else
                                Me%lFlowY(i, j) = -1.0 * CriticalFlow
                            endif
                        endif
 !                   endif
                
                else
                    if (Me%lFlowY(i, j) .lt. 0.0) then
                        if ( abs(Me%lFlowY(i, j))* LocalDT  .gt. Me%myWaterVolumePred(i,j)) then
                            Me%lFlowY(i, j) = - Me%myWaterVolumePred(i,j)  / LocalDT
                        endif
                    elseif (Me%lFlowY(i, j) .gt. 0.0) then
                        if ( Me%lFlowY(i, j)* LocalDT .gt. Me%myWaterVolumePred(i-1,j)) then
                            Me%lFlowY(i, j) =  Me%myWaterVolumePred(i-1,j) / LocalDT
                        endif
                    endif                
                    
                    Me%myWaterVolumePred(i  ,j) = Me%myWaterVolumePred(i,  j) + (Me%lFlowY(i, j) * LocalDT)                        
                    Me%myWaterVolumePred(i-1,j) = Me%myWaterVolumePred(i-1,j) - (Me%lFlowY(i, j) * LocalDT) 
                    
                endif

            else
            
                Me%lFlowY(i, j) = 0.0
            
            endif
            
        enddo
        enddo         
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "DynamicWaveYY")
        
        
    end subroutine DynamicWaveYY
    
    !-------------------------------------------------------------------------
    
!    real function HydraulicRadius(i,j,Direction, level_before, level_after)
!        
!        !Arguments-------------------------------------------------------------
!        integer                                           :: i,j
!        character(len=StringLength)                       :: Direction
!        real                                              :: level_before, level_after
!        !Local-----------------------------------------------------------------
!        real                                              :: WetPerimeter, WaterDepth, MaxBottom
!        real                                              :: Margin1, Margin2
!        integer                                           :: di, dj
!       
!        
!        if(Direction == "X") then
!        
!            !Hydraulic Radius
!            !Wet perimeter, first is bottom
!            WetPerimeter = Me%ExtVar%DYY(i, j)
!
!            !Water Depth consistent with AreaU computed (only water above max bottom)
!            WaterDepth = Me%AreaU(i,j) / Me%ExtVar%DYY(i, j)
!            MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))
!            
!            !to check wich cell to use to use since areaU depends on higher water level
!            if (level_before .gt. level_after) then
!                dj = -1
!            else
!                dj = 0
!            endif
!            
!            !Bottom Difference to adjacent cells (to check existence of “margins” on the side)
!            Margin1 = Me%ExtVar%Topography(i+1, j + dj) - MaxBottom
!            Margin2 = Me%ExtVar%Topography(i-1, j + dj) - MaxBottom
!
!            !if positive than there is a “margin” on the side and friction occurs at wet length
!            !If not basin points than result will be negative.
!            if (Margin1 .gt. 0.0) then
!                WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
!            endif
!            if (Margin2 .gt. 0.0) then
!                WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
!            endif
!
!            HydraulicRadius = Me%AreaU(i, j) / WetPerimeter
!                
!        elseif (Direction == "Y") then
!        
!            !Hydraulic Radius
!            !Wet perimeter, first is bottom
!            WetPerimeter = Me%ExtVar%DXX(i, j)
!
!            !Water Depth consistent with AreaV computed (only water above max bottom)
!            WaterDepth = Me%AreaV(i,j) / Me%ExtVar%DXX(i, j)
!            MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))
!
!            !to check wich cell to use since areaV depends on higher water level
!            if (level_before .gt. level_after) then
!                di = -1
!            else
!                di = 0
!            endif
!
!            !Bottom Difference to adjacent cells (to check existence of “margins” on the side)
!            Margin1 = Me%ExtVar%Topography(i + di,j+1) - MaxBottom
!            Margin2 = Me%ExtVar%Topography(i + di,j-1) - MaxBottom
!
!            !if positive than there is a “margin” on the side and friction occurs at wet length
!            if (Margin1 .gt. 0.0) then
!                WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
!            endif
!            if (Margin2 .gt. 0.0) then
!                WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
!            endif
!            
!            !m = m2 / m
!            HydraulicRadius = Me%AreaV(i, j) / WetPerimeter        
!        endif
!        
!    end function HydraulicRadius
    
    !---------------------------------------------------------------------------
    
    subroutine UpdateWaterLevels(LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: dVol

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "UpdateWaterLevels")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


        !X
        !$OMP PARALLEL PRIVATE(I,J,dVol)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNKI)
        do i = ILB, IUB
        do j = JLB, JUB
            if (Me%ComputeFaceU(i, j) == BasinPoint) then
                
                !dVol
                dVol                       = Me%lFlowX(i, j) * LocalDT
                
                !Updates Water Volume
                Me%myWaterVolume (i, j-1)  = Me%myWaterVolume (i, j-1) - dVol 
                Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   + dVol 

                !Updates Water Column
                Me%myWaterColumn  (i, j-1) = Me%myWaterVolume (i, j-1) / Me%ExtVar%GridCellArea(i, j-1)
                Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                !Updates Water Level
                Me%myWaterLevel (i, j-1)   = Me%myWaterColumn (i, j-1) + Me%ExtVar%Topography(i, j-1)
                Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)

            endif
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL        

        !Y
        !$OMP PARALLEL PRIVATE(I,J,dVol)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNKJ)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ComputeFaceV(i, j) == BasinPoint) then
                
                !dVol
                dVol                      = Me%lFlowY(i, j) * LocalDT
                
                !Updates Water Volume
                Me%myWaterVolume (i-1, j) = Me%myWaterVolume (i-1, j) - dVol 
                Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   + dVol 

                !Updates Water Column
                Me%myWaterColumn (i-1, j) = Me%myWaterVolume (i-1, j) / Me%ExtVar%GridCellArea(i-1, j)
                Me%myWaterColumn (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                !Updates Water Level
                Me%myWaterLevel (i-1, j)  = Me%myWaterColumn (i-1, j) + Me%ExtVar%Topography(i-1, j)
                Me%myWaterLevel (i, j)    = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)

            endif
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "UpdateWaterLevels")

    end subroutine UpdateWaterLevels 

    !--------------------------------------------------------------------------
    
    !old routine where flux is not taken into account level difference
    !and had an error where max volume was compared to a flow and not volume (fixed)
    subroutine RouteDFourPoints_v2
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, it, jt
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: dVol, AverageCellLength, FlowMax
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%DFourSinkPoint(i, j) == BasinPoint .and. Me%LowestNeighborI(i, j) /= null_int .and. &
                Me%myWaterColumn(i,  j) > Me%MinimumWaterColumn) then
            

                it = Me%LowestNeighborI(i, j)
                jt = Me%LowestNeighborJ(i, j)
                               
                !Critical Flow                    
                AverageCellLength  = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0
                !FlowMax = Min(sqrt(Gravity * Me%myWaterColumn(i, j)) *  Me%myWaterColumn(i, j) * AverageCellLength, &
                !              0.1 * Me%myWaterColumn(i, j) * AverageCellLength)
                
                ![m3/s] = [m/s] * [m] * [m]
                FlowMax = sqrt(Gravity * Me%myWaterColumn(i, j)) *  Me%myWaterColumn(i, j) * AverageCellLength
                

                !dVol -> max Critical Flow & Avaliable Volume
                !there was an error in units Flowmax is m3/s and not m3
                !dVol = min(Me%myWaterVolume(i,j), FlowMax)
                ![m3] = [m3/s] * [s]
                dVol = min(Me%myWaterVolume(i,j), FlowMax * Me%ExtVar%DT)
                
                !Updates Water Volume
                Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   - dVol 
                Me%myWaterVolume (it, jt)  = Me%myWaterVolume (it, jt) + dVol 

                !Updates Water Column
                Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterColumn  (it, jt) = Me%myWaterVolume (it, jt) / Me%ExtVar%GridCellArea(it, jt)

                !Updates Water Level
                Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                Me%myWaterLevel (it, jt)   = Me%myWaterColumn (it, jt) + Me%ExtVar%Topography(it, jt)

            
            endif
            
        enddo
        enddo
    

    end subroutine RouteDFourPoints_v2

    !--------------------------------------------------------------------------
    
    !new routine where dh is used and only dh may move not all water column. 
    !and water moves in level gradient and not always doenstream
    subroutine RouteDFourPoints_v3
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, it, jt
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: AverageCellLength, Flow, MaxFlow
        real                                        :: WaveHeight, Celerity, dh
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%DFourSinkPoint(i, j) == BasinPoint .and. Me%LowestNeighborI(i, j) /= null_int)  then         

                it = Me%LowestNeighborI(i, j)
                jt = Me%LowestNeighborJ(i, j)
                
                !topography of cell i,j is always higher than it, jt (is the max bottom)
                WaveHeight =  max(Me%myWaterLevel(i, j), Me%myWaterLevel(it,jt)) - Me%ExtVar%Topography(i,j)
                Celerity   = sqrt(Gravity * WaveHeight)

                if (WaveHeight .gt. Me%MinimumWaterColumn) then
                
                    !Critical Flow                    
                    AverageCellLength  = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0                

                    !dh>0 flow removes water, dh<0 flow brings water
                    dh =  Me%myWaterLevel(i, j) - Me%myWaterLevel(it,jt)
                    
                    !m3/s = m/s * m * m. if dh negative minimum is dh
                    Flow = Celerity *  min(dh, WaveHeight) * AverageCellLength
                    
                    !Max flow is volume given by area * dh
                    !Since it jt has always lower topography if dh negative there is not the
                    !possibility of using an abs(dh) higher than Waveheight (more flux than exists)
                    !if positive dh minimum is positive, if dh negative, negative flux with dh
                    MaxFlow = min(dh, WaveHeight) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                    
                    if (abs(Flow) > abs(MaxFlow)) then
                        Flow = MaxFlow
                    endif
                    
                    Me%iFlowRouteDFour(i,j)    = Flow

                    !Updates Water Volume
                    Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   - Flow *  Me%ExtVar%DT
                    Me%myWaterVolume (it, jt)  = Me%myWaterVolume (it, jt) + Flow *  Me%ExtVar%DT 

                    !Updates Water Column
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)
                    Me%myWaterColumn  (it, jt) = Me%myWaterVolume (it, jt) / Me%ExtVar%GridCellArea(it, jt)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    Me%myWaterLevel (it, jt)   = Me%myWaterColumn (it, jt) + Me%ExtVar%Topography(it, jt)
                
                else
                    Me%iFlowRouteDFour(i,j)    = 0.0
                endif
                
            endif
            
        enddo
        enddo
    

    end subroutine RouteDFourPoints_v3

    !--------------------------------------------------------------------------

    !new routine where flow is computed from manning. 
    subroutine RouteDFourPoints
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, it, jt, di, dj
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Flow, MaxFlow, dx, dy
        real                                        :: AverageCellLengthSink, AverageCellLengthLower
        real                                        :: WaveHeight, sign !, Celerity
        real                                        :: level_up, level_down, CenterDistance
        real                                        :: Slope, VertArea !, WetPerimeter
        real                                        :: HydraulicRadius, OverlandCoef
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%DFourSinkPoint(i, j) == BasinPoint .and. Me%LowestNeighborI(i, j) /= null_int)  then         

                it = Me%LowestNeighborI(i, j)
                jt = Me%LowestNeighborJ(i, j)
                
                !topography of cell i,j is always higher than it, jt (is the max bottom)
                WaveHeight =  max(Me%myWaterLevel(i, j), Me%myWaterLevel(it,jt)) - Me%ExtVar%Topography(i,j)

                if (WaveHeight .gt. Me%MinimumWaterColumn) then                
                    
                    level_up   = Me%myWaterLevel(i, j)
                    level_down = Me%myWaterLevel(it, jt)
                    
                    !diagonal is sqrt of squared distances

                    di = it - i
                    !distance to right cell
                    if (di > 0) then
                        dy = Me%ExtVar%DZY(i, j)
                    else
                       !distance to left cell
                        dy = Me%ExtVar%DZY(i-1, j)
                    endif

                    dj = jt - j
                    if (dj > 0) then
                        dx = Me%ExtVar%DZX(i, j)
                    else
                        dx = Me%ExtVar%DZX(i, j-1)
                    endif
                    
                    CenterDistance = sqrt((dx)**2 + (dy)**2)
                    
                    !Slope
                    if (Me%AdjustSlope) then
                        Slope           = AdjustSlope((level_up - level_down) / CenterDistance)
                    else
                        Slope           = (level_up - level_down) / CenterDistance
                    endif

                    if (Slope.LT.0.0) then
                        sign = -1.0
                    else
                        sign = 1.0
                    end if

                    AverageCellLengthSink   = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0  
                    AverageCellLengthLower  = ( Me%ExtVar%DUX (it, jt) + Me%ExtVar%DVY (it, jt) ) / 2.0              
                    VertArea                = ((AverageCellLengthSink + AverageCellLengthLower) / 2.0) * WaveHeight
                    
                    !Wet perimeter approximation to bottom (no walls effect)
                    !WetPerimeter    = (AverageCellLengthSink + AverageCellLengthLower) / 2.0
                    
                    !Same as wave height. short circuit
                    !HydraulicRadius = VertArea / WetPerimeter
                    HydraulicRadius = WaveHeight
                                 
                    OverlandCoef    = (AverageCellLengthSink * Me%OverlandCoefficient(i, j) +      &
                                       AverageCellLengthLower * Me%OverlandCoefficient(it, jt)) /  &
                                       (AverageCellLengthSink + AverageCellLengthLower)               
                    !
                    !MANNING'S EQUATION -  KINEMATIC WAVE
                    !
                    !m3.s-1 = m2 * m(2/3) / (s.m(-1/3)) = m(8/3) * m(1/3) / s = m3.s-1
                    Flow = sign * VertArea * HydraulicRadius**(2./3.) * sqrt(sign * Slope)          &
                           / OverlandCoef

                    !MaxFlow  = sign * VertArea * sqrt(Gravity * WaveHeight)
                    
                    if (sign > 0.0) then
                        MaxFlow  = min(VertArea * sqrt(Gravity * WaveHeight) * Me%ExtVar%DT, Me%myWaterVolume (i, j)) / Me%ExtVar%DT
                    else
                        MaxFlow  = sign * min(VertArea * sqrt(Gravity * WaveHeight) * Me%ExtVar%DT, Me%myWaterVolume (it, jt)) / &
                            Me%ExtVar%DT
                    endif                    
                    
                    if (abs(Flow) > abs(MaxFlow)) then
                        Flow = MaxFlow
                    endif
                    
                    Me%iFlowRouteDFour(i,j)    = Flow

                    !Updates Water Volume
                    Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   - Flow *  Me%ExtVar%DT
                    Me%myWaterVolume (it, jt)  = Me%myWaterVolume (it, jt) + Flow *  Me%ExtVar%DT 

                    !Updates Water Column
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)
                    Me%myWaterColumn  (it, jt) = Me%myWaterVolume (it, jt) / Me%ExtVar%GridCellArea(it, jt)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    Me%myWaterLevel (it, jt)   = Me%myWaterColumn (it, jt) + Me%ExtVar%Topography(it, jt)
                
                else
                    Me%iFlowRouteDFour(i,j)    = 0.0
                endif
                
            endif
            
        enddo
        enddo
    

    end subroutine RouteDFourPoints

    !--------------------------------------------------------------------------
    
    subroutine ComputeStormWaterModel

#ifdef _SEWERGEMSENGINECOUPLER_
   
        !--------------------------------------------------------------------------
        real(c_double)              :: dt, elapsedTime
        integer                     :: STAT_CALL, n, i, j, xn
        real                        :: Flow, Flow_1, Flow_2, Flow_3, Flow_4, Flow_5, myWaterLevel, Flow_formulation_2, tF
        real                        :: SecondLinkWaterLevel, dh, area, WaterLevelSWMM, sign, HydraulicRadius

        !--------------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ComputeStormWaterModel")

        !Compute inlet potential flow
        call ComputeInletsPotentialFlow
        
        !Sets Inlet Potential Flow and Node Surcharge Depth (MOHIDLand WaterColumn) in SWMM model
        call setInlets_SewerGems
        
        !Sets Manholes Node Surcharge Depth (MOHIDLand WaterColumn) in SWMM model
        call setManholes_SewerGems
        
        !Sets Outfalls Water level (MOHIDLand WaterLevel) in SWMM model
        call setOutFalls_SewerGems

        !Sets HeadWalls Water level (MOHIDLand WaterLevel) in SWMM model and updates FlowEnteringCell
        call setHeadWalls_SewerGems

        !------------------Run SewerGEMS SWMM engine time step-----------------------------------------
        STAT_CALL = SewerGEMSEngine_step_imposed_dt(elapsedTime, Me%ExtVar%DT)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR60'
        !------------------Run SewerGEMS SWMM engine time step-----------------------------------------
        
        
        !Gets outflow coming from SWMM manhole and adds it to the MOHIDLand Cell
        !Updates WaterLevel, WaterColumn and volumes
        call FlowFromManholes
        
        !Gets flow from or to SWMM inlets
        !Computes Effective flow and updates WaterLevel, WaterColumn and volumes. Writes results to inlets output file
        call FlowFromToInlets
        
        !Gets flow to or from outfalls and updates volumes
        call FlowFromToOutfalls
        
        !Gets flow from or to SWMM HeadWalls
        !Computes Effective flow and updates WaterLevel, WaterColumn and volumes. Writes results to headwall output file
        call FlowFromHeadWalls
        
        !set the flow at each cross section node to zero before computing open channel links
        !each open channel link will contribute with flow for its main node link (a cross-section)
        do xn = 1, Me%NumberOfCrossSections
            Me%CrossSections(xn)%Flow = 0.0
        enddo

        !set the flow at each pond node to zero before computing open channel links
        !each open channel link will contribute with flow for its main node link (the pond)
        do xn = 1, Me%NumberOfPonds
            Me%Ponds(xn)%Flow = 0.0
        enddo

        do n = 1, Me%NumberOfOpenChannelLinks

            i   = Me%OpenChannelLinks(n)%I
            j   = Me%OpenChannelLinks(n)%J
            
            if(Me%OpenChannelLinks(n)%TypeOf == OutfallLink_)then

                xn  = Me%OpenChannelLinks(n)%OutfallID

                !Divide total outfall flow by the weight of each grid cell it intercepts (= 1/nCells_InterceptedByOutfall) 
                Flow = Me%Outfalls(xn)%Flow * Me%OpenChannelLinks(n)%Weight

                Me%OpenChannelLinks(n)%Flow = Flow

                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Flow * Me%ExtVar%DT)

                if(Me%myWaterVolume (i, j) < 0.0)then
                    Me%myWaterVolume (i, j) = 0.0
                    Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
                endif

                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            else

                if(Me%OpenChannelLinks(n)%TypeOf == PondLink_)then
                    xn  = Me%OpenChannelLinks(n)%PondID
                else
                    xn  = Me%OpenChannelLinks(n)%CrossSectionID
                endif


                !Get water level for main link node (if link type is Direct_ no more information is needed)
                STAT_CALL = SewerGEMSEngine_getNodeWaterLevel(Me%OpenChannelLinks(n)%LinkID, Me%OpenChannelLinks(n)%WaterLevel)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR140'

                if(Me%OpenChannelLinks(n)%TypeOf == Weighted_)then

                    !Get water level for secondary link nodes
                    STAT_CALL = SewerGEMSEngine_getNodeWaterLevel(Me%OpenChannelLinks(n)%SecondLinkID, SecondLinkWaterLevel)
                    if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR150'

                    !Interpolate between the 2 link nodes
                    Me%OpenChannelLinks(n)%WaterLevel = SecondLinkWaterLevel               + & 
                                                        (Me%OpenChannelLinks(n)%WaterLevel - & 
                                                         SecondLinkWaterLevel)             * &
                                                        Me%OpenChannelLinks(n)%Weight

                endif

                WaterLevelSWMM = Me%OpenChannelLinks(n)%WaterLevel

                Me%OpenChannelLinks(n)%Flow = 0.0
                
                !Set default flow to 0
                Flow = 0.0
                if (WaterLevelSWMM < Me%ExtVar%Topography(i, j)) then      
                    !Only consider flow if level above minimum
                    if (Me%myWaterColumn (i, j) > Me%MinimumWaterColumn) then

                        dh = Me%myWaterColumn (i, j)
                        !Weir equation with 0.4 as coeficient.
                        Flow_1 = 0.4 * Me%OpenChannelLinks(n)%FluxWidth  * sqrt(2.0 * Gravity) * dh ** 1.5
                        
                        Flow_2 = (Me%myWaterColumn (i, j) - Me%MinimumWaterColumn) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                        
                        Flow_3 = (Me%myWaterLevel (i, j) - WaterLevelSWMM) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                        
                        !Maximum empty cell or in between levels (if river level is close to topography)
                        Flow = min(Flow_1, Flow_2, Flow_3)
                        
                        !Check if WaterLevelSWMM is close to topography (using 5mm for no specific reason).
                        !If it is then the equation may soon change and a spike
                        !in flow will appear. So a transition factor between the 2 formulations is computed to smooth it.
                        !The use of mywatercolumn is because the transition can only occur until the maxwater level available
                        !Which is the myWaterColumn.
                        if (Me%ExtVar%Topography(i, j) - WaterLevelSWMM < min(Me%myWaterColumn(i, j),0.005)) then
                            !m2 = m * (m + m)/2
                            area  = Me%OpenChannelLinks(n)%FluxWidth * (((WaterLevelSWMM - Me%ExtVar%Topography(i, j)) +  &
                                    (Me%myWaterLevel (i, j) - Me%ExtVar%Topography(i, j))) / 2.0)
                            
                            !m = m + m
                            HydraulicRadius = Me%myWaterLevel(i, j) - WaterLevelSWMM
                            !m3/s = m2 * (m2/3) * [] / s/(m1/3)
                            Flow_4 = area * HydraulicRadius ** (2./3.) *                        &
                                sqrt(HydraulicRadius/Me%OpenChannelLinks(n)%CellWidth) / Me%OverlandCoefficient(i, j)
                            
                            Flow_5 = (Me%myWaterLevel (i, j) - WaterLevelSWMM) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                            ![] - this factor varies from 1 when the distance is > 0.005 (assumes formulation for flow
                            !from 2D to SWMM) and 0 when the distance is -0.005 (above topography)
                            tF = (Me%ExtVar%Topography(i, j) - WaterLevelSWMM) * 100 + 0.5
                            
                            Flow_formulation_2 = min(Flow_4,Flow_5)
                            
                            Flow = tF * Flow + (1-tF) * Flow_formulation_2
                        endif
                        
                    endif
                else
                    if (abs(WaterLevelSWMM - Me%myWaterLevel(i, j)) > Me%MinimumWaterColumn) then
                        
                        !m2 = m * (m + m)/2
                        area  = Me%OpenChannelLinks(n)%FluxWidth * (((WaterLevelSWMM - Me%ExtVar%Topography(i, j)) +  &
                                (Me%myWaterLevel (i, j) - Me%ExtVar%Topography(i, j))) / 2.0)
                        
                        if (Me%myWaterLevel(i, j) > WaterLevelSWMM) then !Water flows from 2D to channel/pond
                            
                            !m = m + m
                            HydraulicRadius = Me%myWaterLevel(i, j) - WaterLevelSWMM
                            !m3/s = m2 * (m2/3) * [] / s/(m1/3)
                            Flow_1 = area * HydraulicRadius ** (2./3.) *                        &
                                sqrt(HydraulicRadius/Me%OpenChannelLinks(n)%CellWidth) / Me%OverlandCoefficient(i, j)
                            
                            Flow_2 = (Me%myWaterLevel (i, j) - WaterLevelSWMM) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                            
                        else ! Water flows from channel/pond to 2D
                            
                            HydraulicRadius =  WaterLevelSWMM - Me%myWaterLevel(i, j)
                            Flow_1 = -1 * area * HydraulicRadius ** (2./3.) *                        &
                                sqrt(HydraulicRadius/Me%OpenChannelLinks(n)%CellWidth) / Me%OverlandCoefficient(i, j)
                            
                            Flow_2 = (WaterLevelSWMM - Me%myWaterLevel (i, j)) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT

                        end if
                        !m3/s
                        Flow = min (Flow_1, Flow_2)
                    else
                        Flow = 0.0
                    endif
                endif
                
               ! Aplicar aqui a transicao?
                
                Me%OpenChannelLinks(n)%Flow = Flow

                if(Me%OpenChannelLinks(n)%TypeOf == PondLink_)then
                    Me%Ponds(xn)%Flow = Me%Ponds(xn)%Flow + Flow
                else
                    Me%CrossSections(xn)%Flow = Me%CrossSections(xn)%Flow + Flow
                endif
                

                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
                
                if(Me%myWaterVolume (i, j) < 0.0)then
                    Me%myWaterVolume (i, j) = 0.0
                    Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
                endif

                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)
                
            endif
            
        enddo

        do n = 1, Me%NumberOfCrossSections
            
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume  - Me%CrossSections(n)%Flow * Me%ExtVar%DT
            Me%TotalOpenChannelVolume= Me%TotalOpenChannelVolume - Me%CrossSections(n)%Flow * Me%ExtVar%DT
            
            !Set 2D flow to/from cross section node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setSurfaceLinkFlow(Me%CrossSections(n)%SWMM_ID, Me%CrossSections(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR160'
                
        enddo

        do n = 1, Me%NumberOfPonds
                
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume - Me%Ponds(n)%Flow * Me%ExtVar%DT
            Me%TotalPondsVolume      = Me%TotalPondsVolume - Me%Ponds(n)%Flow * Me%ExtVar%DT

            STAT_CALL = SewerGEMSEngine_setSurfaceLinkFlow(Me%Ponds(n)%SWMM_ID, Me%Ponds(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR170'

        enddo
        
        !Get SewerGEMS SWMM current total volume
        STAT_CALL = SewerGEMSEngine_getTotalVolume(Me%Total1DVolume)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR180'

        if(Me%Total1DVolume > Me%MaxTotal1DVolume)then
            Me%MaxTotal1DVolume         = Me%Total1DVolume
            Me%TimeOfMaxTotal1DVolume   = Me%ExtVar%Now - Me%BeginTime            
        endif

        !Get SewerGEMS SWMM current time step 
        STAT_CALL = SewerGEMSEngine_getdt(dt)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR190'
        
        !Store SewerGEMS SWMM current time step
        !Check subroutine ComputeNextTimeStep where 
        Me%StormWaterModelDT = dt
!999 format(a20,1x,3f20.6)
!        write(99,999) TimeToString(Me%ExtVar%Now), elapsedTime *86400.0, Me%ExtVar%DT, Me%StormWaterModelDT

        if(Me%ExtVar%Now .ge. Me%EndTime)then
            if(Me%StormWaterModelDT > 0.0)then
                !Run SewerGEMS SWMM engine just to make sure last output is written in case of time step rounding error
                STAT_CALL = SewerGEMSEngine_step_imposed_dt(elapsedTime, Me%StormWaterModelDT)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR200'
            endif
        endif

        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ComputeStormWaterModel")


#endif _SEWERGEMSENGINECOUPLER_
 
    end subroutine ComputeStormWaterModel
    
    !---------------------------------------------------------------------------
    
    subroutine ComputeStormWaterModel_v0

#ifdef _SEWERGEMSENGINECOUPLER_
   
        !--------------------------------------------------------------------------
        real(c_double)              :: dt, elapsedTime
        integer                     :: STAT_CALL, n, i, j, xn, FlowType, Pond_OutType, Compute_Flow
        real                        :: Flow, Overflow, Flow_1, Flow_2, Flow_3
        real                        :: SecondLinkWaterLevel, dh, Area, WaterLevelSWMM, sign, HydraulicRadius

        !--------------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ComputeStormWaterModel")

        !Compute inlet potential flow
        call ComputeInletsPotentialFlow

        do n = 1, Me%NumberOfInlets

            i = Me%Inlets(n)%I
            j = Me%Inlets(n)%J

            
            !Set SewerGEMS SWMM potential flow into inlet. 
            !SWMM will decide if it can accomodate this flow and if not 
            !it will return the effective flow into the inlet in the next iteration
            STAT_CALL = SewerGEMSEngine_setInletPotentialFlow(Me%Inlets(n)%SWMM_ID, Me%Inlets(n)%PotentialFlow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR10'

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Inlets(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR20'

        enddo

        do n = 1, Me%NumberOfManholes

            i = Me%Manholes(n)%I
            j = Me%Manholes(n)%J

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Manholes(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR30'

        enddo

        do n = 1, Me%NumberOfOutfalls
            
            !I and J of the outfall node location - if outfall is connected to a cross-section 
            !the outfall may intercept multiple cells, but all have same terrain elevation and water elevation
            i = Me%Outfalls(n)%I
            j = Me%Outfalls(n)%J

            
            !Set 2D water level over outfall node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setOutfallWaterLevel(Me%Outfalls(n)%SWMM_ID, Me%myWaterLevel(i,j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR40'
            
        enddo

        do n = 1, Me%NumberOfHeadwalls
            
            !I and J of the headwall node location
            i = Me%Headwalls(n)%I
            j = Me%Headwalls(n)%J

            Me%Headwalls(n)%FlowEnteringCell = 0.0

            !Compute flow entering grid cell 
            if(Me%iFlowX(i,j)   > 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell + Me%iFlowX(i,  j  )
            if(Me%iFlowX(i,j+1) < 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell - Me%iFlowX(i,  j+1)
            if(Me%iFlowY(i,j)   > 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell + Me%iFlowY(i,  j  )
            if(Me%iFlowY(i+1,j) < 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell - Me%iFlowY(i+1,j  )

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Headwalls(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR30'

            !Set 2D water depth for headwall node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setHeadwallWaterDepth(Me%Headwalls(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR50'
            
        enddo


        !Run SewerGEMS SWMM engine time step
        STAT_CALL = SewerGEMSEngine_step_imposed_dt(elapsedTime, Me%ExtVar%DT)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR60'


            
        do n = 1, Me%NumberOfManholes

            i = Me%Manholes(n)%I
            j = Me%Manholes(n)%J

            !Get SewerGEMS SWMM flow out of manhole
            STAT_CALL = SewerGEMSEngine_getNodeOverflow(Me%Manholes(n)%SWMM_ID, Me%Manholes(n)%Outflow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR70'

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Manholes(n)%Outflow * Me%ExtVar%DT
            Me%TotalManholesVolume   = Me%TotalManholesVolume   + Me%Manholes(n)%Outflow * Me%ExtVar%DT
            
            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Me%Manholes(n)%Outflow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            STAT_CALL = SewerGEMSEngine_setNodeSurfaceDepth(Me%Manholes(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR80'

        enddo

        do n = 1, Me%NumberOfInlets

            i = Me%Inlets(n)%I
            j = Me%Inlets(n)%J

            !Get SewerGEMS SWMM flow in or out of inlet
            STAT_CALL = SewerGEMSEngine_getNodeOverflow(Me%Inlets(n)%SWMM_ID, Overflow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR90'

            if(Overflow > 0.0)then
                if(Me%Inlets(n)%PotentialFlow > 0)then
                    Me%Inlets(n)%EffectiveFlow = Overflow - Me%Inlets(n)%PotentialFlow 
                else
                    Me%Inlets(n)%EffectiveFlow = Overflow
                endif
            else
                Me%Inlets(n)%EffectiveFlow = Me%Inlets(n)%PotentialFlow * -1.0
            endif

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT
            Me%TotalInletsVolume     = Me%TotalInletsVolume     + Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT
                
            if(Me%Inlets(n)%OutputResults)then
                if(Me%ExtVar%Now >= Me%Inlets(n)%NextOutputTime)then 

                    write(Me%Inlets(n)%OutputUnit,500)Me%Inlets(n)%OutputTime,          &
                                                      Me%Inlets(n)%FlowEnteringCell,    &
                                                      Me%Inlets(n)%PotentialFlow,       &
                                                      Me%Inlets(n)%EffectiveFlow*-1.0

                    Me%Inlets(n)%NextOutputTime = Me%Inlets(n)%NextOutputTime + Me%Inlets(n)%OutputTimeStep
                    Me%Inlets(n)%OutputTime     = Me%Inlets(n)%OutputTime     + Me%Inlets(n)%OutputTimeStep

                end if
            end if    

            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            STAT_CALL = SewerGEMSEngine_setNodeSurfaceDepth(Me%Inlets(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR120'


        enddo

    500 format(1x, f13.2, 1x, 3(1x, e20.12e3))

        do n = 1, Me%NumberOfOutfalls
            
            !Get SewerGEMS SWMM flow in or out of outfall. 
            STAT_CALL = SewerGEMSEngine_getOutfallFlow(Me%Outfalls(n)%SWMM_ID, Me%Outfalls(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR130'
        
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Outfalls(n)%Flow * Me%ExtVar%DT
            Me%TotalOutfallsVolume   = Me%TotalOutfallsVolume   + Me%Outfalls(n)%Flow * Me%ExtVar%DT
        
        enddo

        do n = 1, Me%NumberOfHeadwalls

            i = Me%Headwalls(n)%I
            j = Me%Headwalls(n)%J

            !Get SewerGEMS Headwall downstream conduit flow
            STAT_CALL = SewerGEMSEngine_getHeadwallDownstreamFlow(Me%Headwalls(n)%SWMM_DownstreamLinkID, Me%Headwalls(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR131'

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume - Me%Headwalls(n)%Flow * Me%ExtVar%DT
            Me%TotalHeadwallsVolume  = Me%TotalHeadwallsVolume  - Me%Headwalls(n)%Flow * Me%ExtVar%DT
            
            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Me%Headwalls(n)%Flow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            if(Me%Headwalls(n)%OutputResults)then
                if(Me%ExtVar%Now >= Me%Headwalls(n)%NextOutputTime)then

                    write(Me%Headwalls(n)%OutputUnit,600)Me%Headwalls(n)%OutputTime,        &
                                                         Me%myWaterLevel  (i, j),           &
                                                         Me%myWaterColumn (i, j),           &
                                                         Me%Headwalls(n)%FlowEnteringCell

                    Me%Headwalls(n)%NextOutputTime = Me%Headwalls(n)%NextOutputTime + Me%Headwalls(n)%OutputTimeStep
                    Me%Headwalls(n)%OutputTime     = Me%Headwalls(n)%OutputTime     + Me%Headwalls(n)%OutputTimeStep

                end if
            end if

        enddo

    600 format(1x, f13.2, 1x, 3(1x, e20.12e3))

        !set the flow at each cross section node to zero before computing open channel links
        !each open channel link will contribute with flow for its main node link (a cross-section)
        do xn = 1, Me%NumberOfCrossSections
            Me%CrossSections(xn)%Flow = 0.0
        enddo

        !set the flow at each pond node to zero before computing open channel links
        !each open channel link will contribute with flow for its main node link (the pond)
        do xn = 1, Me%NumberOfPonds
            Me%Ponds(xn)%Flow = 0.0
        enddo
        
        Me%Counter = Me%Counter + 1

        do n = 1, Me%NumberOfOpenChannelLinks

            i   = Me%OpenChannelLinks(n)%I
            j   = Me%OpenChannelLinks(n)%J
            
            FlowType = 0
            Flow_1 = 0
            Flow_2 = 0
            Flow_3 = 0
            area = 0
            Compute_Flow = 0
            dh = 0
            if(Me%OpenChannelLinks(n)%TypeOf == OutfallLink_)then

                xn  = Me%OpenChannelLinks(n)%OutfallID

                !Divide total outfall flow by the weight of each grid cell it intercepts (= 1/nCells_InterceptedByOutfall) 
                Flow = Me%Outfalls(xn)%Flow * Me%OpenChannelLinks(n)%Weight

                Me%OpenChannelLinks(n)%Flow = Flow

                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Flow * Me%ExtVar%DT)

                if(Me%myWaterVolume (i, j) < 0.0)then
                    Me%myWaterVolume (i, j) = 0.0
                    Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
                endif

                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            else

                if(Me%OpenChannelLinks(n)%TypeOf == PondLink_)then
                    xn  = Me%OpenChannelLinks(n)%PondID
                else
                    xn  = Me%OpenChannelLinks(n)%CrossSectionID
                endif


                !Get water level for main link node (if link type is Direct_ no more information is needed)
                STAT_CALL = SewerGEMSEngine_getNodeWaterLevel(Me%OpenChannelLinks(n)%LinkID, Me%OpenChannelLinks(n)%WaterLevel)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR140'

                if(Me%OpenChannelLinks(n)%TypeOf == Weighted_)then

                    !Get water level for secondary link nodes
                    STAT_CALL = SewerGEMSEngine_getNodeWaterLevel(Me%OpenChannelLinks(n)%SecondLinkID, SecondLinkWaterLevel)
                    if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR150'

                    !Interpolate between the 2 link nodes
                    Me%OpenChannelLinks(n)%WaterLevel = SecondLinkWaterLevel               + & 
                                                        (Me%OpenChannelLinks(n)%WaterLevel - & 
                                                         SecondLinkWaterLevel)             * &
                                                        Me%OpenChannelLinks(n)%Weight

                endif

                WaterLevelSWMM = Me%OpenChannelLinks(n)%WaterLevel

                Me%OpenChannelLinks(n)%Flow = 0.0
                
                if (WaterLevelSWMM < Me%ExtVar%Topography(i, j)) then
                    Pond_OutType = 1        
                    !Only consider flow if level above minimum
                    if (Me%myWaterColumn (i, j) > Me%MinimumWaterColumn) then
                        Compute_Flow = 1

                        dh = Me%myWaterColumn (i, j)
                        !Weir equation with 0.4 as coeficient.
                        Flow_1 = 0.4 * Me%OpenChannelLinks(n)%FluxWidth  * sqrt(2.0 * Gravity) * dh ** 1.5
                        
                        Flow_2 = (Me%myWaterColumn (i, j) - Me%MinimumWaterColumn) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                        
                        Flow_3 = (Me%myWaterLevel (i, j) - WaterLevelSWMM) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                        
                        
                        !Maximum empty cell or in between levels (if river level is close to topography)
                        Flow = min(Flow_1, Flow_2, Flow_3)

                    else
                        Flow = 0.0
                    endif
                else
                    Pond_OutType = 2
                    dh = WaterLevelSWMM - Me%myWaterLevel(i, j)
                    if (abs(WaterLevelSWMM - Me%myWaterLevel(i, j)) > Me%MinimumWaterColumn) then
                        dh = Me%myWaterLevel(i, j) - WaterLevelSWMM
                        Compute_Flow = 1
                        if (dh.LT.0.0) then
                            sign = -1.0
                        else
                            sign = 1.0
                        end if

                        !correta
                        area  = Me%OpenChannelLinks(n)%FluxWidth * (((WaterLevelSWMM - Me%ExtVar%Topography(i, j)) +  &
                                (Me%myWaterLevel (i, j) - Me%ExtVar%Topography(i, j))) / 2.0)
                        
                        !Using the hydraulic radius as difference of heights
                        HydraulicRadius = abs(dh)
                        Flow_1 = area * HydraulicRadius ** (2./3.) * sign *                        &
                                sqrt(ABS(dh)/Me%OpenChannelLinks(n)%CellWidth) / Me%OverlandCoefficient(i, j)
                        
                        Flow_2 = sign * ABS(Me%myWaterLevel (i, j) - WaterLevelSWMM) / &
                               2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                        
                        !Maximum equal levels
                        Flow = sign * min( ABS(Flow_1), ABS(Flow_2))
                        
                        if (Flow == sign*abs(Flow_1)) FlowType = 1
                        if (Flow == sign*abs(Flow_2)) FlowType = 2
                        
                    else
                        Flow = 0.0
                    endif
                endif
                
                Me%OpenChannelLinks(n)%Flow = Flow

                if(Me%OpenChannelLinks(n)%TypeOf == PondLink_)then
                    Me%Ponds(xn)%Flow = Me%Ponds(xn)%Flow + Flow
                else
                    Me%CrossSections(xn)%Flow = Me%CrossSections(xn)%Flow + Flow
                endif
                

                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
                
                if(Me%myWaterVolume (i, j) < 0.0)then
                    Me%myWaterVolume (i, j) = 0.0
                    Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
                endif

                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)
                
            999 format(i5,1x,i3,1x,i3,1x,i1,1x,i1,4f20.6,1x,i1,1x,6f20.6)
                if ((Me%Counter > 15000) .and. (Me%Counter < 16000)) then
                    write(99,999) Me%Counter,i,j,Pond_OutType,Compute_Flow,dh,Me%myWaterLevel(i, j),WaterLevelSWMM,area,&
                        FlowType,Flow,Flow_1,Flow_2,Flow_3,Me%OpenChannelLinks(n)%FluxWidth,Me%ExtVar%Topography(i, j)
                end if
            endif
            
        enddo
        
        write(*,*) "Counter and flow = ", Me%Counter, Me%Ponds(1)%Flow

        do n = 1, Me%NumberOfCrossSections
            
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume  - Me%CrossSections(n)%Flow * Me%ExtVar%DT
            Me%TotalOpenChannelVolume= Me%TotalOpenChannelVolume - Me%CrossSections(n)%Flow * Me%ExtVar%DT
            
            !Set 2D flow to/from cross section node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setSurfaceLinkFlow(Me%CrossSections(n)%SWMM_ID, Me%CrossSections(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR160'
                
        enddo

        do n = 1, Me%NumberOfPonds
                
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume - Me%Ponds(n)%Flow * Me%ExtVar%DT
            Me%TotalPondsVolume      = Me%TotalPondsVolume - Me%Ponds(n)%Flow * Me%ExtVar%DT

            STAT_CALL = SewerGEMSEngine_setSurfaceLinkFlow(Me%Ponds(n)%SWMM_ID, Me%Ponds(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR170'

        enddo
        
        !Get SewerGEMS SWMM current total volume
        STAT_CALL = SewerGEMSEngine_getTotalVolume(Me%Total1DVolume)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR180'

        if(Me%Total1DVolume > Me%MaxTotal1DVolume)then
            Me%MaxTotal1DVolume         = Me%Total1DVolume
            Me%TimeOfMaxTotal1DVolume   = Me%ExtVar%Now - Me%BeginTime            
        endif

        !Get SewerGEMS SWMM current time step 
        STAT_CALL = SewerGEMSEngine_getdt(dt)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR190'
        
        !Store SewerGEMS SWMM current time step
        !Check subroutine ComputeNextTimeStep where 
        Me%StormWaterModelDT = dt
!999 format(a20,1x,3f20.6)
!        write(99,999) TimeToString(Me%ExtVar%Now), elapsedTime *86400.0, Me%ExtVar%DT, Me%StormWaterModelDT

        if(Me%ExtVar%Now .ge. Me%EndTime)then
            if(Me%StormWaterModelDT > 0.0)then
                !Run SewerGEMS SWMM engine just to make sure last output is written in case of time step rounding error
                STAT_CALL = SewerGEMSEngine_step_imposed_dt(elapsedTime, Me%StormWaterModelDT)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR200'
            endif
        endif

        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ComputeStormWaterModel")


#endif _SEWERGEMSENGINECOUPLER_
 
    end subroutine ComputeStormWaterModel_v0

    !---------------------------------------------------------------------------
    
    subroutine ComputeStormWaterModel_original

#ifdef _SEWERGEMSENGINECOUPLER_
   
        !--------------------------------------------------------------------------
        real(c_double)              :: dt, elapsedTime
        integer                     :: STAT_CALL, n, i, j, xn
        real                        :: Flow, Overflow
        real                        :: SecondLinkWaterLevel, dh, Area, WaterLevelSWMM, sign

        !--------------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ComputeStormWaterModel")

        !Compute inlet potential flow
        call ComputeInletsPotentialFlow

        do n = 1, Me%NumberOfInlets

            i = Me%Inlets(n)%I
            j = Me%Inlets(n)%J

            
            !Set SewerGEMS SWMM potential flow into inlet. 
            !SWMM will decide if it can accomodate this flow and if not 
            !it will return the effective flow into the inlet in the next iteration
            STAT_CALL = SewerGEMSEngine_setInletPotentialFlow(Me%Inlets(n)%SWMM_ID, Me%Inlets(n)%PotentialFlow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR10'

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Inlets(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR20'

        enddo

        do n = 1, Me%NumberOfManholes

            i = Me%Manholes(n)%I
            j = Me%Manholes(n)%J

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Manholes(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR30'

        enddo

        do n = 1, Me%NumberOfOutfalls
            
            !I and J of the outfall node location - if outfall is connected to a cross-section 
            !the outfall may intercept multiple cells, but all have same terrain elevation and water elevation
            i = Me%Outfalls(n)%I
            j = Me%Outfalls(n)%J

            
            !Set 2D water level over outfall node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setOutfallWaterLevel(Me%Outfalls(n)%SWMM_ID, Me%myWaterLevel(i,j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR40'
            
        enddo

        do n = 1, Me%NumberOfHeadwalls
            
            !I and J of the headwall node location
            i = Me%Headwalls(n)%I
            j = Me%Headwalls(n)%J

            Me%Headwalls(n)%FlowEnteringCell = 0.0

            !Compute flow entering grid cell 
            if(Me%iFlowX(i,j)   > 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell + Me%iFlowX(i,  j  )
            if(Me%iFlowX(i,j+1) < 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell - Me%iFlowX(i,  j+1)
            if(Me%iFlowY(i,j)   > 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell + Me%iFlowY(i,  j  )
            if(Me%iFlowY(i+1,j) < 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell - Me%iFlowY(i+1,j  )

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Headwalls(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR30'

            !Set 2D water depth for headwall node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setHeadwallWaterDepth(Me%Headwalls(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR50'
            
        enddo


        !Run SewerGEMS SWMM engine time step
        STAT_CALL = SewerGEMSEngine_step_imposed_dt(elapsedTime, Me%ExtVar%DT)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR60'


            
        do n = 1, Me%NumberOfManholes

            i = Me%Manholes(n)%I
            j = Me%Manholes(n)%J

            !Get SewerGEMS SWMM flow out of manhole
            STAT_CALL = SewerGEMSEngine_getNodeOverflow(Me%Manholes(n)%SWMM_ID, Me%Manholes(n)%Outflow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR70'

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Manholes(n)%Outflow * Me%ExtVar%DT
            Me%TotalManholesVolume   = Me%TotalManholesVolume   + Me%Manholes(n)%Outflow * Me%ExtVar%DT
            
            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Me%Manholes(n)%Outflow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            STAT_CALL = SewerGEMSEngine_setNodeSurfaceDepth(Me%Manholes(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR80'

        enddo

        do n = 1, Me%NumberOfInlets

            i = Me%Inlets(n)%I
            j = Me%Inlets(n)%J

            !Get SewerGEMS SWMM flow in or out of inlet
            STAT_CALL = SewerGEMSEngine_getNodeOverflow(Me%Inlets(n)%SWMM_ID, Overflow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR90'

            if(Overflow > 0.0)then
                if(Me%Inlets(n)%PotentialFlow > 0)then
                    Me%Inlets(n)%EffectiveFlow = Overflow - Me%Inlets(n)%PotentialFlow 
                else
                    Me%Inlets(n)%EffectiveFlow = Overflow
                endif
            else
                Me%Inlets(n)%EffectiveFlow = Me%Inlets(n)%PotentialFlow * -1.0
            endif

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT
            Me%TotalInletsVolume     = Me%TotalInletsVolume     + Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT
                
            if(Me%Inlets(n)%OutputResults)then
                if(Me%ExtVar%Now >= Me%Inlets(n)%NextOutputTime)then 

                    write(Me%Inlets(n)%OutputUnit,500)Me%Inlets(n)%OutputTime,          &
                                                      Me%Inlets(n)%FlowEnteringCell,    &
                                                      Me%Inlets(n)%PotentialFlow,       &
                                                      Me%Inlets(n)%EffectiveFlow*-1.0

                    Me%Inlets(n)%NextOutputTime = Me%Inlets(n)%NextOutputTime + Me%Inlets(n)%OutputTimeStep
                    Me%Inlets(n)%OutputTime     = Me%Inlets(n)%OutputTime     + Me%Inlets(n)%OutputTimeStep

                end if
            end if    

            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            STAT_CALL = SewerGEMSEngine_setNodeSurfaceDepth(Me%Inlets(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR120'


        enddo

    500 format(1x, f13.2, 1x, 3(1x, e20.12e3))

        do n = 1, Me%NumberOfOutfalls
            
            !Get SewerGEMS SWMM flow in or out of outfall. 
            STAT_CALL = SewerGEMSEngine_getOutfallFlow(Me%Outfalls(n)%SWMM_ID, Me%Outfalls(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR130'
        
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Outfalls(n)%Flow * Me%ExtVar%DT
            Me%TotalOutfallsVolume   = Me%TotalOutfallsVolume   + Me%Outfalls(n)%Flow * Me%ExtVar%DT
        
        enddo

        do n = 1, Me%NumberOfHeadwalls

            i = Me%Headwalls(n)%I
            j = Me%Headwalls(n)%J

            !Get SewerGEMS Headwall downstream conduit flow
            STAT_CALL = SewerGEMSEngine_getHeadwallDownstreamFlow(Me%Headwalls(n)%SWMM_DownstreamLinkID, Me%Headwalls(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR131'

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume - Me%Headwalls(n)%Flow * Me%ExtVar%DT
            Me%TotalHeadwallsVolume  = Me%TotalHeadwallsVolume  - Me%Headwalls(n)%Flow * Me%ExtVar%DT
            
            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Me%Headwalls(n)%Flow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            if(Me%Headwalls(n)%OutputResults)then
                if(Me%ExtVar%Now >= Me%Headwalls(n)%NextOutputTime)then 
                    !'time WL2D_1 WL2D_2 WL1D_1 WL1D_2 Depth2D_1 Depth2D_2 Depth1D_1 Depth1D_2
                    !write(Me%Headwalls(n)%OutputUnit,600)Me%Headwalls(n)%OutputTime,        &
                    !                                     Me%Headwalls(n)%Headwater2D_1,     &
                    !                                     Me%Headwalls(n)%Headwater2D_2,     &
                    !                                     Me%Headwalls(n)%Headwater1D_1,     &
                    !                                     Me%Headwalls(n)%Headwater2D_2,     &
                    !                                     Me%Headwalls(n)%WaterDepth2D_1,    &
                    !                                     Me%Headwalls(n)%WaterDepth2D_2,    &
                    !                                     Me%Headwalls(n)%WaterDepth1D_1,    &
                    !                                     Me%Headwalls(n)%WaterDepth1D_2,    &
                    !                                     Me%Headwalls(n)%Flow 

                    write(Me%Headwalls(n)%OutputUnit,600)Me%Headwalls(n)%OutputTime,        &
                                                         Me%myWaterLevel  (i, j),           &
                                                         Me%myWaterColumn (i, j),           &
                                                         Me%Headwalls(n)%FlowEnteringCell

                    Me%Headwalls(n)%NextOutputTime = Me%Headwalls(n)%NextOutputTime + Me%Headwalls(n)%OutputTimeStep
                    Me%Headwalls(n)%OutputTime     = Me%Headwalls(n)%OutputTime     + Me%Headwalls(n)%OutputTimeStep

                end if
            end if

        enddo

    600 format(1x, f13.2, 1x, 3(1x, e20.12e3))

        !set the flow at each cross section node to zero before computing open channel links
        !each open channel link will contribute with flow for its main node link (a cross-section)
        do xn = 1, Me%NumberOfCrossSections
            Me%CrossSections(xn)%Flow = 0.0
        enddo

        !set the flow at each pond node to zero before computing open channel links
        !each open channel link will contribute with flow for its main node link (the pond)
        do xn = 1, Me%NumberOfPonds
            Me%Ponds(xn)%Flow = 0.0
        enddo

        do n = 1, Me%NumberOfOpenChannelLinks

            i   = Me%OpenChannelLinks(n)%I
            j   = Me%OpenChannelLinks(n)%J

            if(Me%OpenChannelLinks(n)%TypeOf == OutfallLink_)then

                xn  = Me%OpenChannelLinks(n)%OutfallID

                !Divide total outfall flow by the weight of each grid cell it intercepts (= 1/nCells_InterceptedByOutfall) 
                Flow = Me%Outfalls(xn)%Flow * Me%OpenChannelLinks(n)%Weight

                Me%OpenChannelLinks(n)%Flow = Flow

                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Flow * Me%ExtVar%DT)

                if(Me%myWaterVolume (i, j) < 0.0)then
                    Me%myWaterVolume (i, j) = 0.0
                    Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
                endif

                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            else

                if(Me%OpenChannelLinks(n)%TypeOf == PondLink_)then
                    xn  = Me%OpenChannelLinks(n)%PondID
                else
                    xn  = Me%OpenChannelLinks(n)%CrossSectionID
                endif


                !Get water level for main link node (if link type is Direct_ no more information is needed)
                STAT_CALL = SewerGEMSEngine_getNodeWaterLevel(Me%OpenChannelLinks(n)%LinkID, Me%OpenChannelLinks(n)%WaterLevel)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR140'

                if(Me%OpenChannelLinks(n)%TypeOf == Weighted_)then

                    !Get water level for secondary link nodes
                    STAT_CALL = SewerGEMSEngine_getNodeWaterLevel(Me%OpenChannelLinks(n)%SecondLinkID, SecondLinkWaterLevel)
                    if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR150'

                    !Interpolate between the 2 link nodes
                    Me%OpenChannelLinks(n)%WaterLevel = SecondLinkWaterLevel               + & 
                                                        (Me%OpenChannelLinks(n)%WaterLevel - & 
                                                         SecondLinkWaterLevel)             * &
                                                        Me%OpenChannelLinks(n)%Weight

                endif

                WaterLevelSWMM = Me%OpenChannelLinks(n)%WaterLevel

                Me%OpenChannelLinks(n)%Flow = 0.0

                if (WaterLevelSWMM < Me%ExtVar%Topography(i, j)) then
                            
                    !Only consider flow if level above minimum
                    if (Me%myWaterColumn (i, j) > Me%MinimumWaterColumn) then                            
                                
                        dh = Me%myWaterColumn (i, j)
                        !Weir equation with 0.4 as coeficient.
                        Flow  = 0.4 * Me%OpenChannelLinks(n)%CellWidth  * sqrt(2.0 * Gravity) * dh ** 1.5

                        !Maximum empty cell or in between levels (if river level is close to topography)
                        Flow     = min(Flow, (Me%myWaterColumn (i, j) - Me%MinimumWaterColumn) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT,     &
                                        (Me%myWaterLevel (i, j) - WaterLevelSWMM) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT)
                    else
                        Flow = 0.0
                    endif
                else
                            
                    if (abs(WaterLevelSWMM - Me%myWaterLevel(i, j)) > Me%MinimumWaterColumn) then
                        
                        dh = Me%myWaterLevel(i, j) - WaterLevelSWMM
                        if (dh.LT.0.0) then
                            sign = -1.0
                        else
                            sign = 1.0
                        end if                        
                        area  = Me%OpenChannelLinks(n)%FluxWidth * (WaterLevelSWMM - Me%ExtVar%Topography(i, j)) +  &
                                (Me%myWaterLevel (i, j) - Me%ExtVar%Topography(i, j)) / 2.0
                        Flow  = Area *  Me%OpenChannelLinks(n)%FluxWidth ** (2./3.) * sign *                        &
                                sqrt(ABS(dh)/Me%OpenChannelLinks(n)%CellWidth) / Me%OverlandCoefficient(i, j)
                
                        !Maximum equal levels
                        Flow = sign * min( ABS(Flow), ABS(Me%myWaterLevel (i, j) - WaterLevelSWMM) / &
                               2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT)
                    else
                        Flow = 0.0
                    endif
                endif

                Me%OpenChannelLinks(n)%Flow = Flow

                if(Me%OpenChannelLinks(n)%TypeOf == PondLink_)then
                    Me%Ponds(xn)%Flow = Me%Ponds(xn)%Flow + Flow
                else
                    Me%CrossSections(xn)%Flow = Me%CrossSections(xn)%Flow + Flow
                endif
                

                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
                
                if(Me%myWaterVolume (i, j) < 0.0)then
                    Me%myWaterVolume (i, j) = 0.0
                    Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
                endif

                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            endif

        enddo 

        do n = 1, Me%NumberOfCrossSections
            
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume  - Me%CrossSections(n)%Flow * Me%ExtVar%DT
            Me%TotalOpenChannelVolume= Me%TotalOpenChannelVolume - Me%CrossSections(n)%Flow * Me%ExtVar%DT
            
            !Set 2D flow to/from cross section node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setSurfaceLinkFlow(Me%CrossSections(n)%SWMM_ID, Me%CrossSections(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR160'
                
        enddo

        do n = 1, Me%NumberOfPonds
                
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume - Me%Ponds(n)%Flow * Me%ExtVar%DT
            Me%TotalPondsVolume      = Me%TotalPondsVolume - Me%Ponds(n)%Flow * Me%ExtVar%DT

            STAT_CALL = SewerGEMSEngine_setSurfaceLinkFlow(Me%Ponds(n)%SWMM_ID, Me%Ponds(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR170'

        enddo
        
        !Get SewerGEMS SWMM current total volume
        STAT_CALL = SewerGEMSEngine_getTotalVolume(Me%Total1DVolume)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR180'

        if(Me%Total1DVolume > Me%MaxTotal1DVolume)then
            Me%MaxTotal1DVolume         = Me%Total1DVolume
            Me%TimeOfMaxTotal1DVolume   = Me%ExtVar%Now - Me%BeginTime            
        endif

        !Get SewerGEMS SWMM current time step 
        STAT_CALL = SewerGEMSEngine_getdt(dt)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR190'
        
        !Store SewerGEMS SWMM current time step
        !Check subroutine ComputeNextTimeStep where 
        Me%StormWaterModelDT = dt
!999 format(a20,1x,3f20.6)
!        write(99,999) TimeToString(Me%ExtVar%Now), elapsedTime *86400.0, Me%ExtVar%DT, Me%StormWaterModelDT

        if(Me%ExtVar%Now .ge. Me%EndTime)then
            if(Me%StormWaterModelDT > 0.0)then
                !Run SewerGEMS SWMM engine just to make sure last output is written in case of time step rounding error
                STAT_CALL = SewerGEMSEngine_step_imposed_dt(elapsedTime, Me%StormWaterModelDT)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeStormWaterModel - ModuleRunOff - ERR200'
            endif
        endif

        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ComputeStormWaterModel")


#endif _SEWERGEMSENGINECOUPLER_
 
    end subroutine ComputeStormWaterModel_original
    
     !--------------------------------------------------------------------------

    subroutine ComputeInletsPotentialFlow
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: InletInflow, FlowEnteringCell
        real                                        :: AverageCellLength, y0, dH1, dH2
        integer                                     :: n, iStage


        !Local-----------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do n = 1, Me%NumberOfInlets
                
            i = Me%Inlets(n)%I
            j = Me%Inlets(n)%J

            FlowEnteringCell = 0.0
            Me%Inlets(n)%FlowEnteringCell = 0.0

            !Compute flow entering grid cell 
            if(Me%iFlowX(i,j)   > 0.0) FlowEnteringCell = FlowEnteringCell + Me%iFlowX(i,  j  )
            if(Me%iFlowX(i,j+1) < 0.0) FlowEnteringCell = FlowEnteringCell - Me%iFlowX(i,  j+1)
            if(Me%iFlowY(i,j)   > 0.0) FlowEnteringCell = FlowEnteringCell + Me%iFlowY(i,  j  )
            if(Me%iFlowY(i+1,j) < 0.0) FlowEnteringCell = FlowEnteringCell - Me%iFlowY(i+1,j  )

                    
            if(Me%Inlets(n)%TypeOf == Weir_)then
                    
                AverageCellLength  = (Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j)) / 2.0
                    
                !Considering an average side slope of 5% (1/0.05 = 20) of the street
                y0 = sqrt(2.0*Me%myWaterColumn(i, j)*AverageCellLength / 20.0)
                    
                !When triangle of street is full, consider new head 
                if (y0 * 20.0 > AverageCellLength) then
                    y0 = AverageCellLength / 40.0 + Me%myWaterColumn(i, j)
                endif

                !Q  = L * K * y0^(3/2) * sqrt(g)
                !L  = inlet width = 0.5
                !K  = Coef = 0.2
                !y0 = downstream level
                InletInflow = Me%Inlets(n)%Width * 0.2 * y0**1.5 * sqrt(Gravity)       
                    
            elseif(Me%Inlets(n)%TypeOf == FlowCapture_)then
                    
                if(FlowEnteringCell > 0.0)then
                    InletInflow = FlowEnteringCell * Me%Inlets(n)%CaptureFraction
                else
                    InletInflow = 0.0
                end if
                
            elseif(Me%Inlets(n)%TypeOf == DepthFlowRatingCurve_)then

                if    (Me%myWaterColumn(i, j) < Me%Inlets(n)%RatingCurveStage(1))then
                    !if lower than minimum stage then set to minimum flow
                    InletInflow = Me%Inlets(n)%RatingCurveBelowMin
                elseif(Me%myWaterColumn(i, j) > Me%Inlets(n)%RatingCurveStage(Me%Inlets(n)%RatingCurve_nValues))then
                    !if higher than maximum stage then set to maximum flow
                    InletInflow = Me%Inlets(n)%RatingCurveAboveMax
                else
                    !if within stage levels average the flow
                    iStage = 2 !start at second stage level
                    do iStage = 2, Me%Inlets(n)%RatingCurve_nValues
                        if(Me%myWaterColumn(i, j) .le. Me%Inlets(n)%RatingCurveStage(iStage))then
                            exit
                        endif
                    enddo

                    dH1 = Me%myWaterColumn(i, j) - Me%Inlets(n)%RatingCurveStage(iStage-1) 
                    dH2 = Me%Inlets(n)%RatingCurveStage(iStage) - Me%myWaterColumn(i, j)

                    InletInflow = (dH1 * Me%Inlets(n)%RatingCurveFlow(iStage  )  + &
                                    dH2 * Me%Inlets(n)%RatingCurveFlow(iStage-1)) / &
                                    (Me%Inlets(n)%RatingCurveStage(iStage )        - &
                                    Me%Inlets(n)%RatingCurveStage(iStage-1))

                endif

            elseif(Me%Inlets(n)%TypeOf == FlowFlowRatingCurve_)then

                if    (FlowEnteringCell < Me%Inlets(n)%RatingCurveStage(1))then
                    !if lower than minimum stage then set to minimum flow
                    InletInflow = Me%Inlets(n)%RatingCurveBelowMin
                elseif(FlowEnteringCell > Me%Inlets(n)%RatingCurveStage(Me%Inlets(n)%RatingCurve_nValues))then
                    !if higher than maximum stage then set to maximum flow
                    InletInflow = Me%Inlets(n)%RatingCurveAboveMax
                else
                    !if within stage levels average the flow
                    iStage = 2 !start at second stage level
                    do iStage = 2, Me%Inlets(n)%RatingCurve_nValues
                        if(FlowEnteringCell .le. Me%Inlets(n)%RatingCurveStage(iStage))then
                            exit
                        endif
                    enddo

                    dH1 = FlowEnteringCell - Me%Inlets(n)%RatingCurveStage(iStage-1) 
                    dH2 = Me%Inlets(n)%RatingCurveStage(iStage) - FlowEnteringCell

                    InletInflow = (dH1 * Me%Inlets(n)%RatingCurveFlow(iStage  )  + &
                                    dH2 * Me%Inlets(n)%RatingCurveFlow(iStage-1)) / &
                                    (Me%Inlets(n)%RatingCurveStage(iStage )        - &
                                    Me%Inlets(n)%RatingCurveStage(iStage-1))

                endif
                
            endif

            Me%Inlets(n)%FlowEnteringCell= FlowEnteringCell
            Me%Inlets(n)%PotentialFlow = Min(InletInflow, Me%myWaterVolume(i, j) / Me%ExtVar%DT)

        enddo
            

    end subroutine ComputeInletsPotentialFlow
        
    !--------------------------------------------------------------------------
    
    subroutine setInlets_SewerGems
        !--------------------------------------------------------------------------
        integer                     :: STAT_CALL, n, i, j
        !Begin---------------------------------------------------------------------
        do n = 1, Me%NumberOfInlets
            i = Me%Inlets(n)%I
            j = Me%Inlets(n)%J

            !Set SewerGEMS SWMM potential flow into inlet. 
            !SWMM will decide if it can accomodate this flow and if not 
            !it will return the effective flow into the inlet in the next iteration
            STAT_CALL = SewerGEMSEngine_setInletPotentialFlow(Me%Inlets(n)%SWMM_ID, Me%Inlets(n)%PotentialFlow)
            if (STAT_CALL /= SUCCESS_) stop 'setInlets_SewerGems - ModuleRunOff - ERR10'

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Inlets(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'setInlets_SewerGems - ModuleRunOff - ERR20'
        enddo
    end subroutine setInlets_SewerGems
    
    !--------------------------------------------------------------------------
    
    subroutine setManholes_SewerGems
        !--------------------------------------------------------------------------
        integer                     :: STAT_CALL, n, i, j
        !Begin---------------------------------------------------------------------
        do n = 1, Me%NumberOfManholes
            i = Me%Manholes(n)%I
            j = Me%Manholes(n)%J

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Manholes(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'setManholes_SewerGems - ModuleRunOff - ERR30'
        enddo
    end subroutine setManholes_SewerGems
    
    !--------------------------------------------------------------------------
    
    subroutine setOutFalls_SewerGems
        !--------------------------------------------------------------------------
        integer                     :: STAT_CALL, n, i, j
        !Begin---------------------------------------------------------------------
        do n = 1, Me%NumberOfOutfalls
            !I and J of the outfall node location - if outfall is connected to a cross-section 
            !the outfall may intercept multiple cells, but all have same terrain elevation and water elevation
            i = Me%Outfalls(n)%I
            j = Me%Outfalls(n)%J

            !Set 2D water level over outfall node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setOutfallWaterLevel(Me%Outfalls(n)%SWMM_ID, Me%myWaterLevel(i,j))
            if (STAT_CALL /= SUCCESS_) stop 'setOutFalls_SewerGems - ModuleRunOff - ERR10'
        enddo
    end subroutine setOutFalls_SewerGems
    
    !--------------------------------------------------------------------------
    
    subroutine setHeadWalls_SewerGems
        !--------------------------------------------------------------------------
        integer                     :: STAT_CALL, n, i, j
        !Begin---------------------------------------------------------------------
        do n = 1, Me%NumberOfHeadwalls
            
            !I and J of the headwall node location
            i = Me%Headwalls(n)%I
            j = Me%Headwalls(n)%J

            Me%Headwalls(n)%FlowEnteringCell = 0.0

            !Compute flow entering grid cell 
            if(Me%iFlowX(i,j)   > 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell + Me%iFlowX(i,  j  )
            if(Me%iFlowX(i,j+1) < 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell - Me%iFlowX(i,  j+1)
            if(Me%iFlowY(i,j)   > 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell + Me%iFlowY(i,  j  )
            if(Me%iFlowY(i+1,j) < 0.0) Me%Headwalls(n)%FlowEnteringCell = Me%Headwalls(n)%FlowEnteringCell - Me%iFlowY(i+1,j  )

            STAT_CALL = SewerGEMSEngine_setNodeSurchargeDepth(Me%Headwalls(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'setHeadWalls_SewerGems - ModuleRunOff - ERR10'

            !Set 2D water depth for headwall node to SewerGEMS SWMM 
            STAT_CALL = SewerGEMSEngine_setHeadwallWaterDepth(Me%Headwalls(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'setHeadWalls_SewerGems - ModuleRunOff - ERR20'
        enddo
    end subroutine setHeadWalls_SewerGems
    
    !--------------------------------------------------------------------------
    
    subroutine FlowFromManholes
        !--------------------------------------------------------------------------
        integer                     :: STAT_CALL, n, i, j
        !Begin---------------------------------------------------------------------
        do n = 1, Me%NumberOfManholes

            i = Me%Manholes(n)%I
            j = Me%Manholes(n)%J

            !Get SewerGEMS SWMM flow out of manhole
            STAT_CALL = SewerGEMSEngine_getNodeOverflow(Me%Manholes(n)%SWMM_ID, Me%Manholes(n)%Outflow)
            if (STAT_CALL /= SUCCESS_) stop 'FlowFromManholes - ModuleRunOff - ERR10'

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Manholes(n)%Outflow * Me%ExtVar%DT
            Me%TotalManholesVolume   = Me%TotalManholesVolume   + Me%Manholes(n)%Outflow * Me%ExtVar%DT
            
            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Me%Manholes(n)%Outflow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            STAT_CALL = SewerGEMSEngine_setNodeSurfaceDepth(Me%Manholes(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'FlowFromManholes - ModuleRunOff - ERR20'

        enddo
    end subroutine FlowFromManholes
    
    !--------------------------------------------------------------------------

    subroutine FlowFromToInlets
        !--------------------------------------------------------------------------
        integer                     :: STAT_CALL, n, i, j
        real                        :: Overflow
        !Begin---------------------------------------------------------------------
        do n = 1, Me%NumberOfInlets

            i = Me%Inlets(n)%I
            j = Me%Inlets(n)%J

            !Get SewerGEMS SWMM flow in or out of inlet
            STAT_CALL = SewerGEMSEngine_getNodeOverflow(Me%Inlets(n)%SWMM_ID, Overflow)
            if (STAT_CALL /= SUCCESS_) stop 'FlowFromToInlets - ModuleRunOff - ERR10'

            if(Overflow > 0.0)then
                if(Me%Inlets(n)%PotentialFlow > 0)then
                    Me%Inlets(n)%EffectiveFlow = Overflow - Me%Inlets(n)%PotentialFlow 
                else
                    Me%Inlets(n)%EffectiveFlow = Overflow
                endif
            else
                Me%Inlets(n)%EffectiveFlow = Me%Inlets(n)%PotentialFlow * -1.0
            endif

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT
            Me%TotalInletsVolume     = Me%TotalInletsVolume     + Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT
                
            if(Me%Inlets(n)%OutputResults)then
                if(Me%ExtVar%Now >= Me%Inlets(n)%NextOutputTime)then 

                    write(Me%Inlets(n)%OutputUnit,500)Me%Inlets(n)%OutputTime,          &
                                                      Me%Inlets(n)%FlowEnteringCell,    &
                                                      Me%Inlets(n)%PotentialFlow,       &
                                                      Me%Inlets(n)%EffectiveFlow*-1.0

                    Me%Inlets(n)%NextOutputTime = Me%Inlets(n)%NextOutputTime + Me%Inlets(n)%OutputTimeStep
                    Me%Inlets(n)%OutputTime     = Me%Inlets(n)%OutputTime     + Me%Inlets(n)%OutputTimeStep

                end if
            end if    

            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) + (Me%Inlets(n)%EffectiveFlow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            STAT_CALL = SewerGEMSEngine_setNodeSurfaceDepth(Me%Inlets(n)%SWMM_ID, Me%myWaterColumn (i, j))
            if (STAT_CALL /= SUCCESS_) stop 'FlowFromToInlets - ModuleRunOff - ERR20'
        enddo

        500 format(1x, f13.2, 1x, 3(1x, e20.12e3))
    end subroutine FlowFromToInlets
    
    !--------------------------------------------------------------------------
    
    subroutine FlowFromToOutfalls
        !--------------------------------------------------------------------------
        integer                     :: STAT_CALL, n, i, j
        !Begin---------------------------------------------------------------------
        do n = 1, Me%NumberOfOutfalls
            STAT_CALL = SewerGEMSEngine_getOutfallFlow(Me%Outfalls(n)%SWMM_ID, Me%Outfalls(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'FlowFromToOutfalls - ModuleRunOff - ERR10'
        
            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume + Me%Outfalls(n)%Flow * Me%ExtVar%DT
            Me%TotalOutfallsVolume   = Me%TotalOutfallsVolume   + Me%Outfalls(n)%Flow * Me%ExtVar%DT
        enddo
    end subroutine FlowFromToOutfalls
    
    !--------------------------------------------------------------------------
    
    subroutine FlowFromHeadWalls
        !--------------------------------------------------------------------------
        integer                     :: STAT_CALL, n, i, j
        !Begin---------------------------------------------------------------------
        do n = 1, Me%NumberOfHeadwalls

            i = Me%Headwalls(n)%I
            j = Me%Headwalls(n)%J

            !Get SewerGEMS Headwall downstream conduit flow
            STAT_CALL = SewerGEMSEngine_getHeadwallDownstreamFlow(Me%Headwalls(n)%SWMM_DownstreamLinkID, Me%Headwalls(n)%Flow)
            if (STAT_CALL /= SUCCESS_) stop 'FlowFromToOutfalls - ModuleRunOff - ERR10'

            Me%TotalStormWaterVolume = Me%TotalStormWaterVolume - Me%Headwalls(n)%Flow * Me%ExtVar%DT
            Me%TotalHeadwallsVolume  = Me%TotalHeadwallsVolume  - Me%Headwalls(n)%Flow * Me%ExtVar%DT
            
            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Me%Headwalls(n)%Flow * Me%ExtVar%DT)

            if(Me%myWaterVolume (i, j) < 0.0)then
                Me%myWaterVolume (i, j) = 0.0
                Me%MassError(i, j) = Me%MassError(i, j) + Me%myWaterVolume(i,j)
            endif

            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)

            if(Me%Headwalls(n)%OutputResults)then
                if(Me%ExtVar%Now >= Me%Headwalls(n)%NextOutputTime)then

                    write(Me%Headwalls(n)%OutputUnit,600)Me%Headwalls(n)%OutputTime,        &
                                                         Me%myWaterLevel  (i, j),           &
                                                         Me%myWaterColumn (i, j),           &
                                                         Me%Headwalls(n)%FlowEnteringCell

                    Me%Headwalls(n)%NextOutputTime = Me%Headwalls(n)%NextOutputTime + Me%Headwalls(n)%OutputTimeStep
                    Me%Headwalls(n)%OutputTime     = Me%Headwalls(n)%OutputTime     + Me%Headwalls(n)%OutputTimeStep

                end if
            end if

        enddo

        600 format(1x, f13.2, 1x, 3(1x, e20.12e3))
    end subroutine FlowFromHeadWalls
    
    !--------------------------------------------------------------------------

    subroutine FlowIntoChannels(LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: DifLevel
        real                                        :: Slope, AverageCellLength, dVol
        real                                        :: Area, HydraulicRadius, MaxFlow
        real                                        :: ChannelFreeVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength 
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsVolume
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)


        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR02'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR03'        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR04'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'   


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        !$OMP PARALLEL PRIVATE(I,J, DifLevel, Slope, AverageCellLength, dVol, Area, HydraulicRadius, MaxFlow, ChannelFreeVolume)


        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. ChannelsActiveState(i, j) == BasinPoint) then

                !Checks for Flow from Land -> Channel
                AverageCellLength  = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0

            
                if (ChannelsWaterLevel (i, j) < Me%myWaterLevel(i, j) .and. Me%myWaterColumn(i, j) > Me%MinimumWaterColumn) then

                    if (ChannelsWaterLevel (i, j) > Me%ExtVar%Topography(i, j)) then
                        DifLevel           = Me%myWaterLevel(i, j) - ChannelsWaterLevel (i, j)
                    else
                        DifLevel           = Me%myWaterColumn(i, j)
                    endif

                    !Volume which can enter the channel
                    ChannelFreeVolume = ChannelsMaxVolume(i, j) - ChannelsVolume (i, j)
                
                    !Channel almost empty... put all water into channel    
!                    if (ChannelFreeVolume / ChannelsMaxVolume(i, j) > 0.01) then

                        !Volume to channel: minimum between free volume and current volume in cell
!                        dVol = min(ChannelFreeVolume, Me%myWaterVolume (i, j))

                        !Flow to channel - positive if enters
!                        Me%lFlowToChannels(i, j) = dVol / LocalDT

!                    else
                
                        Slope                      = AdjustSlope(DifLevel / (AverageCellLength / 4.0))

                        Area                       = DifLevel * ChannelsNodeLength(i, j)
                    
                        HydraulicRadius            = Area / ChannelsNodeLength(i, j)
                
                        !Minium between friction (manning) and critical flow
                        Me%lFlowToChannels(i, j)   = min(Area * HydraulicRadius**(2./3.) * sqrt(Slope) /  &
                                                         Me%OverlandCoefficient(i,j), &
                                                         Area * sqrt(Gravity * DifLevel))
                        
                     
                        !MaxFlow = 0.5 * (DifLevel) * Me%ExtVar%GridCellArea(i, j) / LocalDT

                        MaxFlow = sqrt(Gravity * Me%myWaterColumn(i, j)) * Me%myWaterColumn(i, j) * ChannelsNodeLength(i, j)
                   
                        if (Me%lFlowToChannels(i, j) > MaxFlow) then
                            Me%lFlowToChannels(i, j) = MaxFlow
                        endif
                        
!                    endif
                    
                    
                              
                    !dVol
                    dVol                       = Me%lFlowToChannels(i, j) * LocalDT
                    
                    !Updates Water Volume
                    Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   - dVol 
                    
                    !Updates Water Column
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    

                else
                
                    Me%lFlowToChannels(i, j) = 0.0
                
                endif

            
            endif

        enddo
        enddo        
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR10'
        
        
    end subroutine FlowIntoChannels   
    
    !--------------------------------------------------------------------------

    subroutine FlowFromChannels
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: ChannelHeight
        real                                        :: WCR, dVol, VolExcess, NewLevel
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real   , dimension(:, :), pointer           :: ChannelsSurfaceWidth
        real   , dimension(:, :), pointer           :: ChannelsBankSlope
        real   , dimension(:, :), pointer           :: ChannelsBottomLevel
        real                                        :: a0, a1, a2
        real                                        :: x1, x2, MaxFlow, Flow
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)


        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR02'

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR03'

        call GetChannelsBankSlope (Me%ObjDrainageNetwork, ChannelsBankSlope, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR04'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR05'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'        


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        !$OMP PARALLEL PRIVATE(I,J, ChannelHeight, WCR, dVol, VolExcess, NewLevel, a0, a1, a2, x1, x2, MaxFlow)


        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. ChannelsActiveState(i, j) == BasinPoint) then

                if (ChannelsWaterLevel (i, j) > Me%myWaterLevel(i, j)) then
                
                    ChannelHeight = Me%ExtVar%Topography(i, j) - ChannelsBottomLevel(i, j)                                       
                    !ChannelSlope  = (ChannelsTopWidth(i, j) - ChannelsBottomWidth(i, j)) / ChannelHeight
                    !ChannelSurfaceWidth = ChannelsBottomWidth(i,j) + 2.* ChannelSlope * ChannelHeight
                    
                    !Water Column in River above Topo
                    WCR           = ChannelsWaterLevel (i, j) - Me%ExtVar%Topography(i, j)
                    
                    !Volume above Topography
                    VolExcess    = ChannelsBankSlope(i,j) * WCR * WCR * ChannelsNodeLength(i, j)       &
                                    + WCR * ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) +    &
                                    Me%myWaterVolume(i, j)

                    if (ChannelsBankSlope(i,j) <= AlmostZero) then
                        !Rectangular
                        a1 = ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) + Me%ExtVar%GridCellArea(i, j)
                        NewLevel = VolExcess / a1
                        NewLevel = NewLevel + Me%ExtVar%Topography(i, j)

                    else
                        !Trapezoidal - formula resolvente
                        a0 = ChannelsBankSlope(i,j) * ChannelsNodeLength(i, j)
                        a1 = ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) + Me%ExtVar%GridCellArea(i, j)
                        a2 = -1.0 * VolExcess
                                    
                        !Solves Polynominal
                        x1            = (-a1 + sqrt(a1**2. - 4.*a0*a2)) / (2.*a0)
                        x2            = (-a1 - sqrt(a1**2. - 4.*a0*a2)) / (2.*a0)                        

                        if (x1 > 0. .and. x1 < WCR) then
                            NewLevel  = x1 + Me%ExtVar%Topography(i, j)
                        else
                            NewLevel  = x2 + Me%ExtVar%Topography(i, j)
                        endif
                    endif

                    
                    dVol = (NewLevel - Me%myWaterLevel(i, j)) *  Me%ExtVar%GridCellArea(i, j)
                    
!                    Me%iFlowToChannels(i, j)    = -dVol / Me%ExtVar%DT 
                    !Revision David 10/4/10
                    !Usually for each cell flow has only one direction
                    !But may exist the special case where at the beggining channel level is lower than
                    !runoff level, but with the exchange, the channel level got bigger
                    !and a flow addition (subtraction) is needed    
                    !Me%iFlowToChannels(i, j)    = Me%iFlowToChannels(i, j) -dVol / Me%ExtVar%DT     
                    Flow = -dVol / Me%ExtVar%DT
                    
                    !Limits flow to critical one
                    MaxFlow = -1.0 * sqrt(Gravity * WCR) * WCR * ChannelsNodeLength(i, j)
                    
                    if (Flow > MaxFlow) then
                        Flow = MaxFlow
                    endif
                    
                    Me%iFlowToChannels(i, j)    = Me%iFlowToChannels(i, j) + Flow
                    
                    Me%myWaterVolume (i, j)     = Me%myWaterVolume (i, j) - (Flow *  Me%ExtVar%DT)
                    
                    Me%myWaterColumn  (i, j)    = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                    Me%myWaterLevel (i, j)      = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)

                
                endif

            
            endif

        enddo
        enddo        
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBankSlope, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        
   
    
    end subroutine FlowFromChannels
    
    !--------------------------------------------------------------------------
    
    !Same as 6 but with new mapping that is independent on 1D model used (e.g. Drainage Network or SWMM)
    subroutine OverLandChannelInteraction_6_NewMapping()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        !logical                                     :: FoundBankGridPoint
        !integer                                     :: STAT_CALL
        real                                        :: Flow
        !type(T_BankGridPoint), pointer              :: BankGridPoint
        type(T_MarginGridPoint), pointer            :: MarginGridPoint
        integer                                     :: i, j, itarget, jtarget
        real                                        :: RiverLevel
        integer                                     :: BoundaryFaces
        real                                        :: dh, CellWidth, FluxWidth, Area, sign        
        
        
        !output
        call SetMatrixValue(Me%MarginFlowToChannels, Me%Size, null_real)        
        
        
        !Go for MarginGridPoints and compute 1D-2D flow 
        !It will integrate on closest BankGridPoint's for DN or associated NodeGridPoint for OpenMI
        MarginGridPoint => Me%FirstMarginGridPoint
        

        do while (associated(MarginGridPoint))

            !Abreviation
            i = MarginGridPoint%GridI
            j = MarginGridPoint%GridJ
            RiverLevel = MarginGridPoint%RiverLevel
            itarget = MarginGridPoint%GridIIntegrateFlux
            jtarget = MarginGridPoint%GridJIntegrateFlux
            
            !To compute?
            if (RiverLevel > null_real / 2.0) then
                
                CellWidth = (Me%ExtVar%DUX(i, j) + Me%ExtVar%DVY(i, j) ) / 2.0

                BoundaryFaces = Me%ExtVar%BasinPoints(i,j-1) + Me%ExtVar%BasinPoints(i,j+1) + Me%ExtVar%BasinPoints(i-1,j) + Me%ExtVar%BasinPoints(i+1,j)
                FluxWidth = BoundaryFaces * CellWidth
                
                !2 cases:
                ! 1. Level inside channel below topography -> Weir equation
                ! 2. Level inside channel above topography -> Kinematic Wave
                if (RiverLevel < Me%ExtVar%Topography(i, j)) then
                            
                    !Only consider flow if level above minimum
                    if (Me%myWaterColumn (i, j) > Me%MinimumWaterColumn) then                            
                                
                        dh = Me%myWaterColumn (i, j)
                        !Weir equation with 0.4 as coeficient.
                        Flow  = 0.4 * CellWidth  * sqrt(2.0 * Gravity) * dh ** 1.5

                        !Maximum empty cell or in between levels (if river level is close to topography)
                        Flow     = min(Flow, (Me%myWaterColumn (i, j) - Me%MinimumWaterColumn) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT,     &
                                        (Me%myWaterLevel (i, j) - RiverLevel) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT)
                    else
                        Flow = 0.0
                    endif
                else
                            
                    if (ABS(RiverLevel - Me%myWaterLevel(i, j)) > Me%MinimumWaterColumn) then
                        
                        dh = Me%myWaterLevel(i, j) - RiverLevel
                        if (dh.LT.0.0) then
                            sign = -1.0
                        else
                            sign = 1.0
                        end if                        
                        area  = FluxWidth * (RiverLevel - Me%ExtVar%Topography(i, j)) + (Me%myWaterLevel (i, j) - Me%ExtVar%Topography(i, j)) / 2.0
                        Flow  = Area *  FluxWidth ** (2./3.) * sign * sqrt(ABS(dh)/cellwidth) / Me%OverlandCoefficient(i, j)
                
                        !Maximum equal levels
                        Flow = sign * min( ABS(Flow), ABS(Me%myWaterLevel (i, j) - RiverLevel) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT)
                    else
                        Flow = 0.0
                    endif
                endif       
                
            else
                Flow = 0.0
            endif     
            
            !!Important!! flow to channel may have other sources than this, so a sum is needed
            !Put the flow in integrated BankGriPoint - target i and j
            Me%iFlowToChannels(itarget, jtarget) = Me%iFlowToChannels(itarget, jtarget) + Flow

            !output
            Me%MarginFlowToChannels(i, j) = Flow
            
            !Updates Variables
            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)                
            
                    
            MarginGridPoint => MarginGridPoint%Next
                
        enddo         
      

    end subroutine OverLandChannelInteraction_6_NewMapping    
    
    
    
 
    
    
    !--------------------------------------------------------------------------
    
    
    !Same as 5 but simplified code, removed the weir when it goes from river to surface
    !and water moving to the lowest topography when column is low
    subroutine OverLandChannelInteraction_6

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: Flow
        real   , dimension(:, :), pointer           :: ChannelsVolume
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real                                        :: dh, cellwidth, width, area, sign
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        real  , dimension(:, :), pointer            :: ChannelsSurfaceWidth


        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR04'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'   

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR02'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'        

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR03'


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. &   !RiverPoint
                ChannelsActiveState  (i, j) == BasinPoint .and. &   !Active
                ChannelsMaxVolume    (i, j) > 0.0) then             !Not the outlet

                cellwidth = (Me%ExtVar%DUX(i, j) + Me%ExtVar%DVY(i, j) ) / 2.0
                width =     (ChannelsNodeLength(i,j) + cellwidth) / 2.0

                
                
                !2 cases:
                ! 1. Level inside channel below topography -> Weir equation
                ! 2. Level inside channel above topography -> Kinematic Wave
                if (ChannelsWaterLevel(i, j) < Me%ExtVar%Topography(i, j)) then
                            
                    !Only consider flow if level above minimum
                    if (Me%myWaterColumn (i, j) > Me%MinimumWaterColumn) then                            
                                
                        dh = Me%myWaterColumn (i, j)
                        !Weir equation with 0.4 as coeficient.
                        Flow  = 0.4 * cellwidth  * sqrt(2.0 * Gravity) * dh ** 1.5

                        !Maximum empty cell or in between levels (if river level is close to topography)
                        Flow     = min(Flow, (Me%myWaterColumn (i, j) - Me%MinimumWaterColumn) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT,     &
                                        (Me%myWaterLevel (i, j) - ChannelsWaterLevel(i, j)) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT)
                    else
                        Flow = 0.0
                    endif
                else
                            
                    if (ABS(ChannelsWaterLevel(i, j) - Me%myWaterLevel(i, j)) > Me%MinimumWaterColumn) then
                        
                        dh = Me%myWaterLevel(i, j) - ChannelsWaterLevel(i, j);
                        if (dh.LT.0.0) then
                            sign = -1.0
                        else
                            sign = 1.0
                        end if                        
                        area  = width * (ChannelsWaterLevel(i, j) - Me%ExtVar%Topography(i, j)) + (Me%myWaterLevel (i, j) - Me%ExtVar%Topography(i, j)) / 2.0
                        Flow  = area *  width ** (2./3.) * sign * sqrt(ABS(dh)/cellwidth) / Me%OverlandCoefficient(i, j)
                
                        !Maximum equal levels
                        Flow = sign * min(ABS(Flow), ABS(Me%myWaterLevel (i, j) - ChannelsWaterLevel(i, j)) / 2.0 * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT)
                    else
                        Flow = 0.0
                    endif
                endif                
                
                !!Important!! flow to channel may have other sources than this, so a sum is needed
                Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + Flow

                !Updates Variables
                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)                        
                       
                           
            endif

        enddo
        enddo        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR020'        

    end subroutine OverLandChannelInteraction_6    
    
    
    !--------------------------------------------------------------------------
    
    !Method to use celerity as the base for transport water in river runoff interaction
    subroutine OverLandChannelInteraction_2
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: Flow, MaxFlow
        real   , dimension(:, :), pointer           :: ChannelsVolume
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real                                        :: dh, dh_new, WaveHeight, Celerity
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        real  , dimension(:, :), pointer            :: ChannelsSurfaceWidth
        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR04'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'   

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR02'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'        

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR03'


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. &   !RiverPoint
                ChannelsActiveState  (i, j) == BasinPoint .and. &   !Active
                ChannelsMaxVolume    (i, j) > 0.0) then             !Not the outlet


                !dh > 0, flow to channels, dh < 0, flow from channels
                dh         =  Me%myWaterLevel(i, j) - ChannelsWaterLevel(i, j)
                WaveHeight =  max(Me%myWaterLevel(i, j), ChannelsWaterLevel(i, j)) - Me%ExtVar%Topography(i,j)

                Celerity = sqrt(Gravity * WaveHeight)
                
                if (dh > 0) then
                
                    if (Me%myWaterColumn (i, j) > Me%MinimumWaterColumn) then
                    
                        !flux is occuring between dh and with celerity 
                        !m3/s = m/s (celerity) * m2 (Area = (dh * L) * 2)
                        Flow    = Celerity * 2.0 * ChannelsNodeLength(i, j) * min(dh, WaveHeight)
                        
                        !MaxFlow = Me%myWaterVolume (i, j) / Me%ExtVar%DT
                        !if channel level lower than topography - limit is all volume (waveheight is water column)
                        !if channel level higher than topography limit is dh
                        Maxflow = min(dh, WaveHeight) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                    else
                    
                        Flow = 0.0
                        MaxFlow = 0.0
                    
                    endif
                else
                    !Implicit computation of new dh based on celerity dx transport
!                    dh_new = (ChannelsSurfaceWidth(i,j) * dh) /                      &
!                    (ChannelsSurfaceWidth(i,j) + 2 * min (Celerity * Me%ExtVar%DT, 0.5 * Me%ExtVar%DUX(i,j)))
                    
                    !Compute new water height above runoff column based on the distance that water 
                    !will be spread in one dt (surface width + 2 celerity paths - in both ways)
                    ![m] = [m] * [m] / [m] . this is the same as working with volumes where river lenght
                    !would be multiplied in both num and den. dh_new is estimated based on same volume spreading on
                    !wider area
                    dh_new = (ChannelsSurfaceWidth(i,j) * dh) /                      &
                    (ChannelsSurfaceWidth(i,j) + 2 * Celerity * Me%ExtVar%DT)
                    
                    !maximum spread where in one time step all the water above runoff column
                    !will spread along all the cell (DUX)
                    !in case that channel top width is == DUX no flow occurs so this was abandoned
                    !dh_min = (ChannelsSurfaceWidth(i,j) * dh) /                      &
                    !(Me%ExtVar%DUX(i,j))
                    
                    !m3/s = h * L * Length / s
                    Flow    = -1. * (dh_new - dh) * ChannelsSurfaceWidth(i,j) * ChannelsNodeLength(i,j) / Me%ExtVar%DT
                    
                    !MaxFlow = -1. * (dh_min - dh) * ChannelsSurfaceWidth(i,j) * ChannelsNodeLength(i,j) / Me%ExtVar%DT                    
                    !maximum is the channel water above runoff going all to runoff 
                    MaxFlow = dh * ChannelsSurfaceWidth(i,j) * ChannelsNodeLength(i,j) / Me%ExtVar%DT    
                endif

                if (abs(Flow) > abs(MaxFlow)) then
                    Flow = MaxFlow
                endif   
                    
                !!Important!! flow to channel may have other sources than this, so a sum is needed
                Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + Flow

                !Updates Volumes
                !Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - Me%iFlowToChannels    (i, j) * Me%ExtVar%DT
                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
                
                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)

                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)
                           
            endif

        enddo
        enddo        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR020'        

    
    end subroutine OverLandChannelInteraction_2
    
    !--------------------------------------------------------------------------    

    subroutine CheckStability (Restart)

        !Arguments-------------------------------------------------------------
        logical                                     :: Restart

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, n_restart
        real                                        :: variation        

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "CheckStability")

        n_restart = 0
        Restart = .false.
        
        !Verifies negative volumes
        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNKJ)
do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                if (Me%myWaterVolume(i, j) < -AllmostZero) then
!                    write(*,*) '-----'
!                    write(*,*) 'OldVolume ', Me%myWaterVolumeOld(i, j)
!                    write(*,*) 'Negative Volume - Me%myWaterVolume (', i, ', ', j, ') =', Me%myWaterVolume (i, j)
!                    write(*,*) '-----'
                    Restart = .true.                 
                        !exit do1  //Commented this exit because don't know how it begave with OpenMP
                else if (Me%myWaterVolume (i, j) < 0.0) then  
                    Me%myWaterVolume (i, j) = 0.0                 
                endif
            endif
        enddo
        enddo do1
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL


        if ((.not. Restart) .and. Me%CV%Stabilize) then

            
            !$OMP PARALLEL PRIVATE(I,J,variation)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNKJ) REDUCTION(+ : n_restart)
do2:        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
                if (Me%StabilityPoints(i, j) == BasinPoint) then
            
                    if ((.not. Me%CV%CheckDecreaseOnly) .or. (Me%myWaterVolumeOld(i, j) > Me%myWaterVolume(i, j))) then
                        if (Me%myWaterVolumeOld(i, j) / Me%ExtVar%GridCellArea(i, j) >= Me%CV%MinimumValueToStabilize) then
                            
                            variation = abs(Me%myWaterVolume(i, j) - Me%myWaterVolumeOld(i, j)) / Me%myWaterVolumeOld(i, j)
                            
                            if (variation > Me%CV%StabilizeFactor) then
                                !Debug routine - may be usefull for using in debug situation
                                !call DebugStability (i,j,variation)                                
                                
                                n_restart = n_restart + 1
                                
                            endif
                        endif
                    endif
                endif
            enddo
            enddo do2
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL


            if (n_restart > Me%CV%MinToRestart) then
                Restart = .true.
            endif                                 

        endif
        
        if (Restart) then        
            Me%CV%NextNiteration = max(int(Me%CV%NextNiteration * Me%CV%DTSplitFactor), Me%CV%NextNiteration + 1)
                 
            if (Me%CV%NextNiteration >= Me%CV%MaxIterations) then
                 write(*,*)'Number of iterations above maximum: ', Me%CV%NextNiteration
                 stop 'CheckStability - ModuleRunoff - ERR010'
            endif                          
        endif           
        
        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "CheckStability")

    end subroutine CheckStability
 
    !--------------------------------------------------------------------------

    subroutine DebugStability(i,j, variation)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: I, J
        real                                        :: variation
        !Local-----------------------------------------------------------------
        character (Len = 5)                         :: str_i, str_j
        character (Len = 15)                        :: str_1, str_2, str_3
        character (len = StringLength)              :: string_to_be_written 
        
        write(str_i, '(i3)') i 
        write(str_j, '(i3)') j
        write(str_1, '(ES10.3)') Me%myWaterVolumeOld(I,J)  
        write(str_2, '(ES10.3)') Me%myWaterVolume(I,J)   
        write(str_3, '(ES10.3)') variation                            
        
        string_to_be_written = ' '//str_i//','//str_j//' '//str_1//' '//str_2//' '//str_3
        
        call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)           
    
    
    end subroutine DebugStability
    
    !--------------------------------------------------------------------------
    !FUNCTION: This routine updates the water level, column and volume at each iteration
    !step if convergence is not met
    !INPUT: old water column
    !RESULT: updated 2D fields of water level, column and volume to initial values
    subroutine LocalWaterColumn (WaterColumn)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :), pointer              :: WaterColumn

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                Me%myWaterVolume(i, j) = WaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)

                Me%myWaterColumn(i, j) = WaterColumn(i, j)

                Me%myWaterLevel (i, j) = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
            endif
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL        
        
    end subroutine LocalWaterColumn            

    !--------------------------------------------------------------------------

    subroutine IntegrateFlow (LocalDT, SumDT)

        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT, SumDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: CHUNK
        real(8)                                     :: sumDischarge
        
        !----------------------------------------------------------------------

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        sumDischarge = Me%TotalDischargeFlowVolume

        !$OMP PARALLEL PRIVATE(I,J) 

        !Integrates along X Directions
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%iFlowX(i, j) = (Me%iFlowX(i, j) * SumDT + Me%lFlowX(i, j) * LocalDT) / &
                              (SumDT + LocalDT)
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Integrates along Y Directions
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%iFlowY(i, j) = (Me%iFlowY(i, j) * SumDT + Me%lFlowY(i, j) * LocalDT) / &
                              (SumDT + LocalDT)
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Integrates Flow to Channels
        if (Me%ObjDrainageNetwork /= 0) then
           !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%iFlowToChannels(i, j) = (Me%iFlowToChannels(i, j) * SumDT + Me%lFlowToChannels(i, j) * LocalDT) / &
                                           (SumDT + LocalDT)
            enddo
            enddo
            !$OMP END DO NOWAIT
        endif
        
        !Integrates Flow At boundary 
        !if (Me%ImposeBoundaryValue) then
        !   !$OMP DO SCHEDULE(DYNAMIC, CHUNK) REDUCTION(+:sumBoundary)
        !    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        !    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        !        Me%iFlowBoundary(i, j) = (Me%iFlowBoundary(i, j) * SumDT + Me%lFlowBoundary(i, j) * LocalDT) / &
        !                                 (SumDT + LocalDT)
        !
        !        sumBoundary = sumBoundary + (Me%iFlowBoundary(i, j) * LocalDT)
        !    enddo
        !    enddo
        !    !$OMP END DO NOWAIT
        !endif

        !Integrates Flow Discharges
        if (Me%Discharges) then
           !$OMP DO SCHEDULE(DYNAMIC, CHUNK) REDUCTION(+:sumDischarge)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB                
                Me%iFlowDischarge(i, j) = (Me%iFlowDischarge(i, j) * SumDT + Me%lFlowDischarge(i, j) * LocalDT) / &
                                          (SumDT + LocalDT)
                sumDischarge = sumDischarge + (Me%lFlowDischarge(i, j) * LocalDT)
            enddo
            enddo
            !$OMP END DO NOWAIT
        endif

        !$OMP END PARALLEL        

        Me%TotalDischargeFlowVolume = sumDischarge

    end subroutine IntegrateFlow

    !--------------------------------------------------------------------------
    
!    !New formulation: not instantaneous but flow computed based on celerity- updated to allow water in
!    subroutine ImposeBoundaryValue_Old ()
!    
!        !Arguments-------------------------------------------------------------
!        !real                                        :: LocalDT
!
!        !Local-----------------------------------------------------------------
!        integer                                     :: i, j, di, dj
!        integer                                     :: ILB, IUB, JLB, JUB
!        real                                        :: dh, dVOl
!        !logical                                     :: NearBoundary
!        real                                        :: AreaZX, AreaZY !, Width
!        real                                        :: WaveHeight, Celerity, MaxFlow
!
!        !Routes water outside the watershed if water is higher then a given treshold values
!        ILB = Me%WorkSize%ILB
!        IUB = Me%WorkSize%IUB
!        JLB = Me%WorkSize%JLB
!        JUB = Me%WorkSize%JUB
!        
!        !Default is zero
!        Me%iFlowBoundary = 0.0
!        
!        !Sets Boundary values
!        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!            if (Me%ExtVar%BasinPoints(i, j)  == BasinPoint .and.          &      !BasinPoint
!                Me%BoundaryCells     (i,j)   == BasinPoint .and.          &      !BoundaryPoints
!                Me%ExtVar%Topography (i, j)  < Me%MaxDtmForBoundary .and. &      !Low land point where to imposes BC
!                Me%myWaterLevel      (i, j)  > Me%BoundaryValue) then            !Level higher then imposed level
!
!!                !Check if near a boundary point
!!                NearBoundary = .false.
!!                do dj = -1, 1
!!                do di = -1, 1
!!                    if (Me%ExtVar%BasinPoints(i+di, j+dj) == 0) then
!!                        NearBoundary = .true.
!!                    endif
!!                enddo
!!                enddo
!!
!!                if (NearBoundary) then
!
!                    !Necessary Variation in height - always positive because only evaluates cell as so
!                    dh = Me%myWaterLevel (i, j) - Me%BoundaryValue
!                    
!                    if (dh > Me%MinimumWaterColumn) then
!                        
!!                        !Cell Width
!!                        Width                  = (Me%ExtVar%DYY(i, j) + Me%ExtVar%DXX(i, j)) / 2.0
!                        
!                        !celerity is limited by water column size and not dh
!                        WaveHeight = Me%myWaterColumn(i, j) 
!
!!                        !Area for Flow
!!                        Area                   = Width * min(dh, WaveHeight)
!                        
!                        Celerity = sqrt(Gravity * WaveHeight)
!
!                        !Flow to set cell equal to Boundary Value
!                        !m3/s                  = 
!!                        Me%lFlowBoundary(i, j) = Min(0.5 * dh * Me%ExtVar%GridCellArea(i, j) / LocalDT,         &
!!                                                     0.5 * Area * sqrt(Gravity * dh))
!
!                        !U direction - use middle area because in closed faces does not exist AreaU
!                        !flow negative (exiting runoff)
!                        AreaZX = Me%ExtVar%DVY(i,j) * Me%myWaterColumn(i,j)
!                        do dj = 0, 1
!                            if ((Me%ComputeFaceU(i,j+dj) == 0)) then
!                                Me%iFlowBoundary(i, j) = Me%iFlowBoundary(i, j) - AreaZX * Celerity
!                            endif
!                        enddo
!
!                        !V direction - use middle area because in closed faces does not exist AreaV
!                        AreaZY = Me%ExtVar%DUX(i,j) * Me%myWaterColumn(i,j)                      
!                        do di = 0, 1
!                            if ((Me%ComputeFaceV(i+di,j) == 0)) then
!                                Me%iFlowBoundary(i, j) = Me%iFlowBoundary(i, j) - AreaZY * Celerity
!                            endif
!                        enddo
!                        
!                        !cant remove more than up to boundary or water column if boundary lower than topography 
!                        !negative flow 
!                        !m3/s = m * m2 / s
!                        MaxFlow = - min(dh, WaveHeight) *  Me%ExtVar%GridCellArea(i, j) / Me%ExtVar%DT
!                        
!                        !m3/s = m2 * m/s 
!                        Me%iFlowBoundary(i, j)  = max (Me%iFlowBoundary(i, j), MaxFlow)
!                        
!                        !dVol
!                        dVol = Me%iFlowBoundary(i, j) * Me%ExtVar%DT
!                            
!                        !Updates Water Volume
!                        Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   + dVol 
!                            
!                        !Updates Water Column
!                        Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)  / Me%ExtVar%GridCellArea(i, j)
!
!                        !Updates Water Level
!                        Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)  + Me%ExtVar%Topography(i, j)
!                        
!                    else
!                    
!                        Me%iFlowBoundary(i, j) = 0.0
!                    
!                    endif
!
!!                endif
!           
!            endif
!        enddo
!        enddo
!
!    
!    end subroutine ImposeBoundaryValue_Old
    
    !--------------------------------------------------------------------------

    subroutine ImposeBoundaryValue ()
    
        !Arguments-------------------------------------------------------------
        !real                                        :: LocalDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, di, dj
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: dh, dVOl
        !logical                                     :: NearBoundary
        !real                                        :: AreaZX, AreaZY !, Width
        real                                        :: WaveHeight, Celerity, MaxFlow
        real                                        :: InflowVolume, OutflowVolume

        !Routes water outside the watershed if water is higher then a given treshold values
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        !Default is zero
        Me%iFlowBoundary = 0.0

        Me%BoundaryFlowVolume = 0.0
        InflowVolume          = 0.0 
        OutflowVolume         = 0.0

        
        !Sets Boundary values
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j)  == BasinPoint .and.          &      !BasinPoint
                Me%BoundaryCells     (i,j)   == BasinPoint .and.          &      !BoundaryPoints
                !Me%ExtVar%Topography (i, j)  < Me%MaxDtmForBoundary .and. &      !Low land point where to imposes BC
                (Me%AllowBoundaryInflow                    .or.           &
                 Me%myWaterLevel      (i, j)  > Me%WaterLevelBoundaryValue(i,j))) then          !Level higher then imposed level

                !Necessary Variation in height (negative if higher outside)
                dh = Me%myWaterLevel (i, j) - Me%WaterLevelBoundaryValue(i,j)
                    
                if (abs(dh) > Me%MinimumWaterColumn) then
                    
                    !celerity is limited by water column on the flow direction (higher level)
                    WaveHeight = max(Me%myWaterLevel(i, j), Me%WaterLevelBoundaryValue(i,j)) - Me%ExtVar%Topography (i, j)
                    
                    ![m/s] = [m/s2 * m]^1/2 = [m/s]
                    Celerity = sqrt(Gravity * WaveHeight)

                    !U direction - use middle area because in closed faces does not exist AreaU
                    !if dh negative, minimum is dh, if positive minimum is dh until boundary
                    !level is lower than terrain and wave height is used (water column)
                    !dh negative flow positive (entering runoff)
                    do dj = 0, 1
                        if ((Me%ComputeFaceU(i,j+dj) == 0)) then
                            ![m3/s] = [m3/s] - [m] * [m] * [m/s]
                            Me%iFlowBoundary(i, j) = Me%iFlowBoundary(i, j) - Me%ExtVar%DVY(i,j) *      &
                                                     min(dh, WaveHeight) * Celerity
                        endif
                    enddo

                    !V direction - use middle area because in closed faces does not exist AreaV
                    do di = 0, 1
                        if ((Me%ComputeFaceV(i+di,j) == 0)) then
                            ![m3/s] = [m3/s] - [m] * [m] * [m/s]
                            Me%iFlowBoundary(i, j) = Me%iFlowBoundary(i, j) - Me%ExtVar%DUX(i,j) *      &
                                                     min(dh, WaveHeight) * Celerity
                        endif
                    enddo
                    
                    !cant remove more than up to boundary or water column if boundary lower than topography 
                    !or add more up to boundary level if boundary level higher
                    !m3/s = m * m2 / s
                    MaxFlow = - min(dh, WaveHeight) *  Me%ExtVar%GridCellArea(i, j) / Me%ExtVar%DT
                    
                    !Me%iFlowBoundary(i, j)  = max (Me%iFlowBoundary(i, j), MaxFlow)
                    if (abs(Me%iFlowBoundary(i, j)) > abs(MaxFlow)) then
                        Me%iFlowBoundary(i, j) = MaxFlow
                    endif                    
                    
                    !dVol
                    dVol = Me%iFlowBoundary(i, j) * Me%ExtVar%DT
                        
                    !Updates Water Volume
                    Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   + dVol 

                    Me%BoundaryFlowVolume     = Me%BoundaryFlowVolume + dVol

                    if(dVol > 0.0)then
                        Me%TotalBoundaryInflowVolume    = Me%TotalBoundaryInflowVolume + dVol
                    else
                        Me%TotalBoundaryOutflowVolume   = Me%TotalBoundaryOutflowVolume + dVol
                    endif
                        
                    !Updates Water Column
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)  / Me%ExtVar%GridCellArea(i, j)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)  + Me%ExtVar%Topography(i, j)
                    
                else
                
                    Me%iFlowBoundary(i, j) = 0.0
                
                endif

            endif
        enddo
        enddo

        Me%TotalBoundaryFlowVolume  = Me%TotalBoundaryFlowVolume + Me%BoundaryFlowVolume

    end subroutine ImposeBoundaryValue
    
    !--------------------------------------------------------------------------
    
    !old formulation with instantaneous going to boundary level
    subroutine ImposeBoundaryValue_v2
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j !, di, dj
        integer                                     :: ILB, IUB, JLB, JUB
!        logical                                     :: NearBoundary
        real                                        :: OldVolume, dVol

        !Routes water outside the watershed if water is higher then a given treshold values
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        Me%BoundaryFlowVolume = 0.0
        
        !Sets Boundary values
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j)  == BasinPoint       .and.    &   !BasinPoint
                Me%BoundaryCells     (i,j)   == BasinPoint       .and.    &   !BoundaryPoints
                !Me%ExtVar%Topography (i, j)  < Me%MaxDtmForBoundary .and. &   !Low land point where to imposes BC
                (Me%AllowBoundaryInflow                          .or.     &
                 Me%myWaterLevel      (i, j)  > Me%WaterLevelBoundaryValue(i,j))) then        !Level higher then imposed level

!                !Check if near a boundary point
!                NearBoundary = .false.
!                do dj = -1, 1
!                do di = -1, 1
!                    if (Me%ExtVar%BasinPoints(i+di, j+dj) == 0) then
!                        NearBoundary = .true.
!                    endif
!                enddo
!                enddo
!
!                if (NearBoundary) then

                    !Necessary Variation in height
                    Me%myWaterLevel (i, j) = max(Me%WaterLevelBoundaryValue(i,j), Me%ExtVar%Topography (i, j))

                    !Updates Water Column
                    Me%myWaterColumn(i, j) = Me%myWaterLevel (i, j) - Me%ExtVar%Topography (i, j) 
                    
                    !Updates Volume and BoundaryFlowVolume
                    OldVolume              = Me%myWaterVolume(i, j)

                    dVol = Me%myWaterVolume(i, j) - OldVolume

                    if(dVol > 0.0)then
                        Me%TotalBoundaryInflowVolume    = Me%TotalBoundaryInflowVolume + dVol
                    else
                        Me%TotalBoundaryOutflowVolume   = Me%TotalBoundaryOutflowVolume + dVol
                    endif
                    
                    !m3 = m * m2
                    Me%myWaterVolume(i, j) = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
                    
                    !m3 = m3 + (m3 - m3)
                    Me%BoundaryFlowVolume  = Me%BoundaryFlowVolume + (OldVolume - Me%myWaterVolume(i, j))
                    
                    !m3/s = m3 / s - always negative exiting runoff
                    Me%iFlowBoundary(i, j) = (Me%myWaterVolume(i, j) - OldVolume) / Me%ExtVar%DT
                    
!                endif
           
            endif
        enddo
        enddo

        Me%TotalBoundaryFlowVolume  = Me%TotalBoundaryFlowVolume + Me%BoundaryFlowVolume

    
    end subroutine ImposeBoundaryValue_v2       
    !--------------------------------------------------------------------------
    
    subroutine ComputeCenterValues 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: CHUNK
        real                                        :: FlowX, FlowY

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ComputeCenterValues")

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
       
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
            
        if(.not. Me%ExtVar%Distortion) then

            if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ComputeCenterValues - CenterVelocity")

            !$OMP PARALLEL PRIVATE(I,J,FlowX,FlowY)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    FlowX = (Me%iFlowX(i, j) + Me%iFlowX(i, j+1)) / 2.0
                    FlowY = (Me%iFlowY(i, j) + Me%iFlowY(i+1, j)) / 2.0
                    
                    Me%CenterFlowX(i, j) = FlowX * Me%GridCosAngleX + FlowY * Me%GridCosAngleY
                    Me%CenterFlowY(i, j) = FlowX * Me%GridSinAngleX + FlowY * Me%GridSinAngleY
                
                    if (Me%myWaterColumn (i,j) > AllmostZero) then
                        Me%CenterVelocityX (i, j) = Me%CenterFlowX (i,j) / ( Me%ExtVar%DYY(i, j) * Me%myWaterColumn (i,j) )
                        Me%CenterVelocityY (i, j) = Me%CenterFlowY (i,j) / ( Me%ExtVar%DXX(i, j) * Me%myWaterColumn (i,j) )
                    else
                        Me%CenterVelocityX(i,j) = 0.0
                        Me%CenterVelocityY(i,j) = 0.0
                    end if
                else
                    Me%CenterFlowX(i,j)     = 0.0
                    Me%CenterFlowY(i,j)     = 0.0
                    Me%CenterVelocityX(i,j) = 0.0
                    Me%CenterVelocityY(i,j) = 0.0
                endif

            enddo
            enddo
            !$OMP END DO NOWAIT 
            !$OMP END PARALLEL

            if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ComputeCenterValues - CenterVelocity")

        else
            
            !$OMP PARALLEL PRIVATE(I,J,FlowX,FlowY)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    FlowX = (Me%iFlowX(i, j) + Me%iFlowX(i, j+1)) / 2.0
                    FlowY = (Me%iFlowY(i, j) + Me%iFlowY(i+1, j)) / 2.0
                    
                    Me%CenterFlowX(i, j) = FlowX * cos(Me%ExtVar%RotationX(i, j)) + FlowY * cos(Me%ExtVar%RotationY(i, j))
                    Me%CenterFlowY(i, j) = FlowX * sin(Me%ExtVar%RotationX(i, j)) + FlowY * sin(Me%ExtVar%RotationY(i, j))
                
                    if (Me%myWaterColumn (i,j) > AllmostZero) then
                        Me%CenterVelocityX (i, j) = Me%CenterFlowX (i,j) / ( Me%ExtVar%DYY(i, j) * Me%myWaterColumn (i,j))
                        Me%CenterVelocityY (i, j) = Me%CenterFlowY (i,j) / ( Me%ExtVar%DXX(i, j) * Me%myWaterColumn (i,j))
                    else
                        Me%CenterVelocityX(i,j) = 0.0
                        Me%CenterVelocityY(i,j) = 0.0
                    end if

                else

                    Me%CenterFlowX(i,j)     = 0.0
                    Me%CenterFlowY(i,j)     = 0.0
                    Me%CenterVelocityX(i,j) = 0.0
                    Me%CenterVelocityY(i,j) = 0.0

                endif

            enddo
            enddo
            !$OMP END DO NOWAIT 
            !$OMP END PARALLEL
            
        endif

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ComputeCenterValues - Modulus")

        
        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                Me%FlowModulus(i, j) = sqrt (Me%CenterFlowX(i, j)**2. + Me%CenterFlowY(i, j)**2.)
                
                if (Me%myWaterColumn (i,j) > AllmostZero) then
                    Me%VelocityModulus (i, j) = sqrt (Me%CenterVelocityX(i, j)**2.0 + Me%CenterVelocityY(i, j)**2.0)
                else
                    Me%VelocityModulus(i,j) = 0.0
                end if

                if(Me%Output%WriteMaxFlowModulus) then
                    if (Me%FlowModulus(i, j) > Me%Output%MaxFlowModulus(i, j)) then
                        Me%Output%MaxFlowModulus(i, j) = Me%FlowModulus(i, j)
                    end if
                end if
            else
                Me%FlowModulus(i,j)     = 0.0
                Me%VelocityModulus(i,j) = 0.0
            endif

        enddo
        enddo
        !$OMP END DO NOWAIT 
        !$OMP END PARALLEL
        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ComputeCenterValues - Modulus")

        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ComputeCenterValues")
        
    end subroutine ComputeCenterValues 

    !--------------------------------------------------------------------------
    
    subroutine ComputeNextDT (Niter)

        !Arguments-------------------------------------------------------------
        integer                                     :: Niter        
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL, CHUNK
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: nextDTCourant, aux
        real                                        :: nextDTVariation, MaxDT
        logical                                     :: VariableDT
        real                                        :: CurrentDT

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ComputeNextDT")


        call GetVariableDT(Me%ObjTime, VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNextDT - ModuleRunOff -  ERR010'

        call GetMaxComputeTimeStep(Me%ObjTime, MaxDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNextDT - ModuleRunOff -  ERR020'

        nextDTCourant   = -null_real
        nextDTVariation = -null_real
        
        if (VariableDT) then

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB
            
            if (Me%CV%LimitDTCourant) then
                        
                !$OMP PARALLEL PRIVATE(I,J,aux)
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK) REDUCTION(MIN:nextDTCourant)
                do j = JLB, JUB
                do i = ILB, IUB
                        
                    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint .and. Me%myWaterColumn (i,j) > Me%MinimumWaterColumn) then

                        !vel = sqrt(Gravity * Me%myWaterColumn (i,j))
                        
                        !if (vel .gt. 0.0) then
                        
                            !spatial step, in case of dx = dy, dist = sqrt(2) * dx
                            !dist = sqrt ((Me%ExtVar%DZX(i, j)**2.0) + (Me%ExtVar%DZY(i, j)**2.0))
                            aux = sqrt ((Me%ExtVar%DZX(i, j)**2.0) + (Me%ExtVar%DZY(i, j)**2.0)) * &
                                   Me%CV%MaxCourant / sqrt(Gravity * Me%myWaterColumn (i,j)) 
                        
                            nextDTCourant = min(nextDTCourant, aux)
                            
                        !endif
                            
                    endif

                enddo
                enddo
                !$OMP END DO NOWAIT 
                !$OMP END PARALLEL

            endif
            
            if (Niter == 1) then
            
                nextDTVariation = Me%ExtVar%DT * Me%CV%DTFactorUp
                Me%CV%NextNiteration = Niter
                
            elseif (Niter <= Me%CV%MinIterations) then                            
            
                if (Niter > Me%CV%LastGoodNiteration) then

                    nextDTVariation = Me%ExtVar%DT
                    Me%CV%NextNiteration = Niter

                else
                
                    nextDTVariation = Me%ExtVar%DT * Me%CV%DTFactorUp
                    Me%CV%NextNiteration = Niter

                endif
                
            else
            
                if (Niter >= Me%CV%StabilizeHardCutLimit) then
                
                    nextDTVariation = (Me%ExtVar%DT / Niter) * Me%CV%MinIterations
                    Me%CV%NextNiteration = Me%CV%MinIterations
                    
                elseif (Niter > Me%CV%LastGoodNiteration) then
                
                    nextDTVariation = Me%ExtVar%DT / Me%CV%DTFactorDown
                    Me%CV%NextNiteration = max(int(nextDTVariation / Me%CV%CurrentDT), 1)
                    
                else
                
                    nextDTVariation = Me%ExtVar%DT
                    Me%CV%NextNiteration = max(min(int(Niter / Me%CV%DTSplitFactor), Niter - 1), 1)
                    
                endif 
                               
            endif
            
            CurrentDT = nextDTVariation / Me%CV%NextNiteration                                     
                      
            Me%CV%NextDT = min(min(nextDTVariation, nextDTCourant), MaxDT)
            
            if (Me%CV%NextDT < nextDTVariation) then                
                Me%CV%NextNiteration = max(int(Me%CV%NextDT/CurrentDT), 1)
            endif
                       
        else
        
            Me%CV%NextDT = Me%ExtVar%DT
            Me%CV%NextNiteration = Niter            
        
        endif
        
        Me%CV%LastGoodNiteration = Niter
        Me%CV%CurrentDT          = Me%CV%NextDT / Me%CV%NextNiteration

        if (Me%StormWaterModel) then
            if (Me%StormWaterModelDT < Me%CV%NextDT) then        
                Me%CV%NextDT = Me%StormWaterModelDT            
            endif
        end if

        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ComputeNextDT")

    end subroutine ComputeNextDT

    !--------------------------------------------------------------------------
    
    subroutine RunOffOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,  dimension(:), pointer                :: AuxFlow
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        real, dimension(6)  , target                :: AuxTime
        real, dimension(:)  , pointer               :: TimePointer
        integer                                     :: dis
        logical                                     :: dbg = .false.

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "RunOffOutput")

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
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                 "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D      = TimePointer,                    &
                                 OutputNumber = Me%OutPut%NextOutPut,           &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR20'

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR30'
            
            
                        
            !Writes mask with grid cells above minimum water column height
            call HDF5WriteData   (Me%ObjHDF5, "/Grid/OpenPoints",              &
                                  "OpenPoints", "-",                            &
                                  Array2D      = Me%OpenPoints,                 &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR031'
            
            !Writes Flow values
            !Writes the Water Column - should be on runoff
            call HDF5WriteData   (Me%ObjHDF5, "/Results/water column",          &
                                  "water column", "m",                          &
                                  Array2D      = Me%MyWaterColumn,              &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR040'
            
            
            

       
            !Writes the Water Level
            call HDF5WriteData   (Me%ObjHDF5, "/Results/water level",           &
                                  "water level", "m",                           &
                                  Array2D      = Me%MyWaterLevel,               &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR050'
            


            !Writes Flow X
            call HDF5WriteData   (Me%ObjHDF5,                                       &
                                  "/Results/flow X",                                &
                                  "flow X",                                         &   
                                  "m3/s",                                           &
                                  Array2D      = Me%CenterFlowX,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR60'

            
            !Writes Flow Y
            call HDF5WriteData   (Me%ObjHDF5,                                       &
                                  "/Results/flow Y",                                &
                                  "flow Y",                                         &   
                                  "m3/s",                                           &
                                  Array2D      = Me%CenterFlowY,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR70'

             !Writes Flow Modulus
            call HDF5WriteData   (Me%ObjHDF5,                                       &
                                  "/Results/"//trim(GetPropertyName (FlowModulus_)),&
                                  trim(GetPropertyName (FlowModulus_)),             &   
                                  "m3/s",                                           &
                                  Array2D      = Me%FlowModulus,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR80'

             !Writes Velocity X 
            call HDF5WriteData   (Me%ObjHDF5,                                          &
                                  "/Results/"//trim(GetPropertyName (VelocityU_)),     &
                                  trim(GetPropertyName (VelocityU_)),                  &
                                  "m/s",                                               &
                                  Array2D      = Me%CenterVelocityX,                   &
                                  OutputNumber = Me%OutPut%NextOutPut,                 &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR90'

             !Writes Velocity Y 
            call HDF5WriteData   (Me%ObjHDF5,                                          &
                                  "/Results/"//trim(GetPropertyName (VelocityV_)),     &
                                  trim(GetPropertyName (VelocityV_)),                  &
                                  "m/s",                                               &
                                  Array2D      = Me%CenterVelocityY,                   &
                                  OutputNumber = Me%OutPut%NextOutPut,                 &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR100'

            !Writes Velocity Modulus
            call HDF5WriteData   (Me%ObjHDF5,                                                &
                                  "/Results/"//trim(GetPropertyName (VelocityModulus_)),     &
                                  trim(GetPropertyName (VelocityModulus_)),                  &
                                  "m/s",                                                     &
                                  Array2D      = Me%VelocityModulus,                         &
                                  OutputNumber = Me%OutPut%NextOutPut,                       &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR110'

            
            if (Me%Use1D2DInteractionMapping) then
                
                if (dbg) then
                    
                !River level from 1D model in river nodes
                call HDF5WriteData   (Me%ObjHDF5, "//Results/node river level", &
                                        "node river level", "m",                  &
                                        Array2D      = Me%NodeRiverLevel,         &
                                        OutputNumber = Me%OutPut%NextOutPut,      &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR191'   
                    
                !River level from 1D model interpolated in margin 2D cells
                call HDF5WriteData   (Me%ObjHDF5, "//Results/margin river level", &
                                        "margin river level", "m",                  &
                                        Array2D      = Me%MarginRiverlevel,         &
                                        OutputNumber = Me%OutPut%NextOutPut,        &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR192'      
                    
                !Flow to river in margin 2D cells
                call HDF5WriteData   (Me%ObjHDF5, "//Results/margin flow to river", &
                                        "margin flow to river", "m3/s",               &
                                        Array2D      = Me%MarginFlowToChannels,       &
                                        OutputNumber = Me%OutPut%NextOutPut,          &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR193'        
                    
                !Flow to river integrated in river cells
                call HDF5WriteData   (Me%ObjHDF5, "//Results/node flow to river", &
                                        "node flow to river", "m3/s",               &
                                        Array2D      = Me%iFlowToChannels,          &
                                        OutputNumber = Me%OutPut%NextOutPut,        &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR194'
                
                end if
                    
            endif
           
            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR200'

            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

        endif
        

        if (Me%OutPut%TimeSerieDischON) then
        
            if (Me%ExtVar%Now >=  Me%OutPut%NextOutPutDisch) then
            
                do dis = 1, Me%OutPut%DischargesNumber
       
                    allocate(AuxFlow(Me%OutPut%TS_Numb_DischProp))
                    
                    AuxFlow(1:Me%OutPut%TS_Numb_DischProp) = Me%OutPut%TimeSerieDischProp(dis,1:Me%OutPut%TS_Numb_DischProp)
                    
                    call WriteTimeSerieLine(Me%OutPut%TimeSerieDischID(dis), AuxFlow, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR210'
                    
                    deallocate(AuxFlow)
                    
                enddo    

                Me%OutPut%NextOutPutDisch = Me%OutPut%NextOutPutDisch + Me%Output%OutPutDischDT
                
            endif                
        endif              

         if (MonitorPerformance) call StopWatch ("ModuleRunOff", "RunOffOutput")
        
    end subroutine RunOffOutput

    !--------------------------------------------------------------------------
    
    subroutine OutputTimeSeries

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !----------------------------------------------------------------------
       
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%MyWaterLevel,                                   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR01'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%MyWaterColumn,                                  &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR02'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%CenterFlowX,                                    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR03'

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%CenterFlowY,                                    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR04'

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%FlowModulus,                                    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR05'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%CenterVelocityX,                                &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR06'

        
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%CenterVelocityY,                                &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR07'

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%VelocityModulus,                                &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR08'
     
        if (Me%Use1D2DInteractionMapping) then
            
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D = Me%NodeRiverLevel,                             &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR20'
            
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D = Me%MarginRiverlevel,                           &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR30'
            
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D = Me%MarginFlowToChannels,                       &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR40'
            
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D = Me%iFlowToChannels,                            &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputTimeSeries - ModuleRunoff - ERR50'            
            
        endif
   
    end subroutine OutputTimeSeries
    
    !--------------------------------------------------------------------------

    subroutine ComputeBoxesWaterFluxes

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, CHUNK, i, j
        real, dimension(:,:), pointer           :: WaterVolume
        !----------------------------------------------------------------------
       
        if (MonitorPerformance) call StartWatch ("ModuleRunoff", "ComputeBoxesWaterFluxes")
        
        call BoxDif(Me%ObjBoxDif,                                                    &
                    Me%iFlowX,                                                       &
                    Me%iFlowY,                                                       &
                    'runoff_water',                                                  &
                    Me%ExtVar%BasinPoints,                                           &
                    STAT = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                 &
           stop 'Subroutine ComputeBoxesWaterFluxes - ModuleRunoff. ERR01'
        
        
        allocate(WaterVolume(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
        WaterVolume = null_real
        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExtVar%BasinPoints(i, j) == 1) then
                
                ![m3] = [m] * [m2]
                WaterVolume(i,j) = Me%myWaterColumn(i,j) * Me%ExtVar%GridCellArea(i,j)
                
            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        
        call BoxDif(Me%ObjBoxDif,                                                    &
                    WaterVolume,                                                     &
                    'runoff_water',                                                  &
                    Me%ExtVar%BasinPoints,                                           &
                    STAT = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                 &
           stop 'Subroutine ComputeBoxesWaterFluxes - ModuleRunoff. ERR02'

        deallocate (WaterVolume)
        
        
        if (MonitorPerformance) call StopWatch ("ModuleRunoff", "ComputeBoxesWaterFluxes")

    end subroutine ComputeBoxesWaterFluxes

    !--------------------------------------------------------------------------

    subroutine OutputFlooding

        !Locals----------------------------------------------------------------
        integer                                 :: ILB,IUB, JLB, JUB, i, j
        integer                                 :: STAT_CALL
        real, dimension(:,:), pointer           :: ChannelsWaterLevel, ChannelsVelocity
        real, dimension(:,:), pointer           :: ChannelsTopArea
        real                                    :: SumArea, WeightedVelocity
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------        

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Begin-----------------------------------------------------------------

   
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
   
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                !Water Column of overland flow
                if (Me%myWaterColumn(i, j) > Me%Output%MaxWaterColumn(i, j)) then
                    Me%Output%MaxWaterColumn(i, j) = Me%myWaterColumn(i, j)
                    
                    !Velocity at MaxWater column
                    Me%Output%VelocityAtMaxWaterColumn(i,j) =  Me%VelocityModulus (i, j)

                    Me%Output%TimeOfMaxWaterColumn(i,j) = Me%ExtVar%Now - Me%BeginTime
                   
                endif
                if ((Me%myWaterColumn(i, j) * (Me%VelocityModulus (i, j) + Me%Output%FloodRiskVelCoef))        &
                     > Me%Output%MaxFloodRisk(i,j)) then
                    Me%Output%MaxFloodRisk(i,j) = Me%myWaterColumn(i, j)                                       &
                                                  * (Me%VelocityModulus (i, j) + Me%Output%FloodRiskVelCoef)
                endif

            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        if (Me%ObjDrainageNetwork /= 0) then

            call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputFlooding - ModuleRunOff - ERR01' 
            
            call GetChannelsTopArea  (Me%ObjDrainageNetwork, ChannelsTopArea, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputFlooding - ModuleRunOff - ERR02'              

            call GetChannelsVelocity  (Me%ObjDrainageNetwork, ChannelsVelocity, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputFlooding - ModuleRunOff - ERR03'             
            
            do j = JLB, JUB
            do i = ILB, IUB
       
                !Water Column of River Network
                if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) then
                    if (ChannelsWaterLevel(i, j) - Me%ExtVar%Topography(i, j) > Me%Output%MaxWaterColumn(i, j)) then
                        Me%Output%MaxWaterColumn(i, j) = ChannelsWaterLevel(i, j) - Me%ExtVar%Topography(i, j)
                        
                        SumArea = Me%ExtVar%GridCellArea(i,j) + ChannelsTopArea(i,j)
                        
                        WeightedVelocity = (Me%VelocityModulus (i, j) * Me%ExtVar%GridCellArea(i,j) +   &
                                            ChannelsVelocity(i,j) * ChannelsTopArea(i,j) ) / SumArea
                        
                        !weighted velocity with river
                        Me%Output%VelocityAtMaxWaterColumn(i,j) = WeightedVelocity
                        
                        if ((Me%Output%MaxWaterColumn(i, j) *  (WeightedVelocity + Me%Output%FloodRiskVelCoef))     &
                              > Me%Output%MaxFloodRisk(i,j)) then
                            Me%Output%MaxFloodRisk(i,j) = Me%Output%MaxWaterColumn(i, j)                            &
                                                          * (WeightedVelocity + Me%Output%FloodRiskVelCoef)
                        endif
                    endif                                        
                    
                endif

            enddo
            enddo

            call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVelocity, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputFlooding - ModuleRunOff - ERR04'

            call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputFlooding - ModuleRunOff - ERR05'

            call UnGetDrainageNetwork  (Me%ObjDrainageNetwork, ChannelsTopArea, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputFlooding - ModuleRunOff - ERR06'                
            
        endif


    end subroutine OutputFlooding

    !---------------------------------------------------------------------------
    
    subroutine OutputFloodPeriod

        !Locals----------------------------------------------------------------
        integer                                 :: ILB,IUB, JLB, JUB, i, j, n
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------        

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Begin-----------------------------------------------------------------

   
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        if(Me%Output%nFloodPeriodLimits > 0)then
            
            !$OMP PARALLEL PRIVATE(I,J, n)
            do n = 1, Me%Output%nFloodPeriodLimits

                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do j = JLB, JUB
                do i = ILB, IUB
   
                    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                        if (Me%myWaterColumn(i, j) > Me%Output%FloodPeriodWaterColumnLimits(n)) then
                            Me%Output%FloodPeriods(i, j, n) = Me%Output%FloodPeriods(i, j, n) + Me%ExtVar%DT
                        endif

                    endif

                enddo
                enddo
                !$OMP END DO NOWAIT
            enddo
            !$OMP END PARALLEL

        else
           
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB
   
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                    !Flooded cell
                    if (Me%myWaterColumn(i, j) > Me%Output%FloodPeriodWaterColumnLimit) then
                        Me%Output%FloodPeriod(i, j) = Me%Output%FloodPeriod(i, j) + Me%ExtVar%DT
                    endif

                endif

            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

        endif

    end subroutine OutputFloodPeriod
   
    !-----------------------------------------------------------------------------    
    
    subroutine OutputFloodArrivalTime

        !Locals----------------------------------------------------------------
        integer                                 :: ILB,IUB, JLB, JUB, i, j
        integer                                 :: CHUNK
        real                                    :: Sum

        !Begin-----------------------------------------------------------------        

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Begin-----------------------------------------------------------------

   
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        Sum   = 0.0

        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK) REDUCTION(+:sum)
        do j = JLB, JUB
        do i = ILB, IUB
   
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                if(Me%myWaterColumn(i, j) > Me%Output%FloodArrivalWaterColumnLimit)then

                    Sum = Sum + Me%ExtVar%GridCellArea(i,j)

                    if(Me%Output%FloodArrivalTime(i, j) < 0.0)then
                        Me%Output%FloodArrivalTime(i, j) = Me%ExtVar%Now - Me%BeginTime
                    endif
                endif
            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        Me%Output%TotalFloodedArea = Sum

        if(Me%Output%TotalFloodedArea > Me%Output%MaxTotalFloodedArea)then
            Me%Output%MaxTotalFloodedArea       = Me%Output%TotalFloodedArea
            Me%Output%TimeOfMaxTotalFloodedArea = Me%ExtVar%Now - Me%BeginTime
        endif

    end subroutine OutputFloodArrivalTime

    !-----------------------------------------------------------------------------    

!    subroutine  WriteChannelsLevelData
!
!        !Local-------------------------------------------------------------------
!        integer                                                 :: ILB,IUB, JLB, JUB
!        integer                                                 :: STAT_CALL,i,j
!        integer, dimension (:,:), pointer                       :: ChannelsID
!        character(len=StringLength), dimension (:,:), pointer   :: ChannelsStationName
!
!        !------------------------------------------------------------------------
!
!        call GetChannelsID  (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR01'
!
!        call GetChannelsStationName  (Me%ObjDrainageNetwork, ChannelsStationName, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR02'
!
!        call GetRiverPoints (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR02a'
!
!
!        !GetNodeID
!        !GetNodeStationName
!
!        open(UNIT=UnitMax, FILE=Me%MaxWaterColumnFile, ACTION='WRITE', STATUS='REPLACE', IOSTAT=STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR03'
!
!
!
!        write(UnitMax,*) 'NodeID     MaxWaterColumn DateTime            StationName'
!
!        ILB = Me%WorkSize%ILB
!        IUB = Me%WorkSize%IUB
!        JLB = Me%WorkSize%JLB
!        JUB = Me%WorkSize%JUB
!        
!        do j = JLB, JUB
!        do i = ILB, IUB
!
!            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) &
!                write(UnitMax,100) ChannelsID(i,j), Me%MaxWaterColumn(i,j), Me%MaxWaterColumnTime(i,j), &
!                trim(adjustl(ChannelsStationName(i,j)))
!
!        enddo
!        enddo
!       
!        close(UnitMax)
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR04'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsStationName, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR05'
!
!        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR05a'
!
!
!        100 format(I10,1x, f16.3, 1x, A19, 1x, A)   
!
!    end subroutine  WriteChannelsLevelData


    !----------------------------------------------------------------------------
   
    real function AdjustSlope (Slope)
    
        !Arguments--------------------------------------------------------------
        real                                    :: Slope
        real                                    :: sign

        !Slope correction given by City of Albuquerque, 1997, p.22-26
        !http://www.hkh-friend.net.np/rhdc/training/lectures/HEGGEN/Tc_3.pdf


        if (Slope .LT. 0.0) then
            sign = -1.0
        else
            sign = 1.0
        end if

        Slope = abs (Slope)
        
        if (Slope.GE.0.04) then
            Slope = 0.05247 + 0.06363 * Slope - 0.182 * exp (-62.38 * Slope)
        end if
        
        AdjustSlope = sign * Slope
        

    end function AdjustSlope

    !----------------------------------------------------------------------------

    subroutine CalculateTotalStoredVolume

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, CHUNK
        real(8)                                     :: Sum

        !Begin-----------------------------------------------------------------

        CHUNK = ChunkJ

        Sum = 0.0

        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK) REDUCTION(+:sum)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
            if (Me%ExtVar%BasinPoints(i, j) == 1) then
                !m3 = m3  + m3
                Sum = Sum + Me%MyWaterVolume(i, j)

            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        Me%TotalStoredVolume        = Sum
        Me%VolumeStoredInSurface    = Sum

        if(Me%TotalStoredVolume > Me%MaxTotal2DVolume)then
            Me%MaxTotal2DVolume       = Me%TotalStoredVolume
            Me%TimeOfMaxTotal2DVolume = Me%ExtVar%Now - Me%BeginTime
        endif

        if(Me%StormWaterModel)then
            Me%Total1D2DVolume = Me%TotalStoredVolume + Me%Total1DVolume
        else
            Me%Total1D2DVolume = Me%TotalStoredVolume
        endif

        if(Me%Total1D2DVolume > Me%MaxTotal1D2DVolume)then
            Me%MaxTotal1D2DVolume       = Me%Total1D2DVolume
            Me%TimeOfMaxTotal1D2DVolume = Me%ExtVar%Now - Me%BeginTime
        endif

!999 format(a20,1x,3f20.6)
!        write(99,999) TimeToString(Me%ExtVar%Now), Me%Total1DVolume, Me%TotalStoredVolume, Me%Total1D2DVolume


    end subroutine CalculateTotalStoredVolume

    !--------------------------------------------------------------------------
    
    subroutine WriteFinalFile_Bin(IsFinalFile)
        
        !Arguments-------------------------------------------------------------
        logical                                     :: IsFinalFile
        !Local-----------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: FinalFile
        integer                                     :: STAT_CALL
        character(LEN = PathLength)                 :: FileName
        
        !----------------------------------------------------------------------

        !Gets Date
        call ExtractDate(Me%ExtVar%Now, Year_File, Month_File, Day_File,               &
                         Hour_File, Minute_File, Second_File)
        
        
        !if (Me%ExtVar%Now == Me%EndTime) then
        if (IsFinalFile .or. Me%Output%RestartOverwrite) then
            FileName = Me%Files%FinalFile
        else
            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%ExtVar%Now))//".fin")
        endif            
        
        call UnitsManager(FinalFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFileOld - ModuleRunoff - ERR01'

        open(Unit = FinalFile, File = FileName, Form = 'UNFORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFileOld - ModuleRunoff - ERR02'

        !Writes Date
        write(FinalFile) Year_File, Month_File, Day_File, Hour_File, Minute_File,       &
                         Second_File

        write(FinalFile)Me%myWaterColumn
        
        call UnitsManager(FinalFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFileOld - ModuleRunoff - ERR03'

    end subroutine WriteFinalFile_Bin

    !------------------------------------------------------------------------
    
    subroutine WriteFinalFile_Hdf(IsFinalFile)
        
        !Arguments-------------------------------------------------------------
        logical                                     :: IsFinalFile
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !integer                                     :: OutPutNumber
        integer                                     :: HDF5_CREATE
        character(LEN = PathLength)                 :: FileName
        integer                                     :: ObjHDF5
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        type (T_Time)                               :: Actual           
        !Begin----------------------------------------------------------------

        !Gets a pointer to Topography
        call GetGridData        (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR00'

        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR01'

          !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Checks if it's at the end of the run 
        !or !if it's supposed to overwrite the final HDF file
        !if ((Me%ExtVar%Now == Me%ExtVar%EndTime) .or. Me%Output%RestartOverwrite) then
        if (IsFinalFile .or. Me%Output%RestartOverwrite) then

            filename = trim(Me%Files%FinalFile)

        else

            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%ExtVar%Now))//".fin")

        endif


        ObjHDF5 = 0
        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5,                                                     &
                            trim(filename),                                              &
                            HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalFile - ModuleRunoff - ERR10'

        Actual   = Me%ExtVar%Now
         
        call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
        !Writes Time
        TimePtr => AuxTime
        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR11'

        call HDF5WriteData  (ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                             Array1D = TimePtr, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR12'


        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5,                                &
                              Me%WorkSize%ILB,                           &
                              Me%WorkSize%IUB,                           &
                              Me%WorkSize%JLB,                           &
                              Me%WorkSize%JUB,                           &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR02'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR25'

        !Writes the Grid
        call HDF5WriteData   (ObjHDF5, "/Grid", "Topography", "m",           &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR05'

        !WriteBasinPoints
        call HDF5WriteData   (ObjHDF5, "/Grid", "BasinPoints", "-",          &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR07'



        call HDF5SetLimits   (ObjHDF5,                                &
                                Me%WorkSize%ILB,                           &
                                Me%WorkSize%IUB,                           &
                                Me%WorkSize%JLB,                           &
                                Me%WorkSize%JUB,                           &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR10'

        call HDF5WriteData   (ObjHDF5,                                    &
                                "/Results/water column",                  &
                                "water column",                           &
                                "m",                                      &
                                Array2D = Me%myWaterColumn,               &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR14'
                

        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR030'

        !Unget
        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR90'  

        !UnGets Topography
        call UnGetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR100'
            
        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR0190'            

    end subroutine WriteFinalFile_Hdf

    !----------------------------------------------------------------------------    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine KillRunOff(RunOffID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: RunOffID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, STAT_CALL, dis, n, i, j
        integer                             :: RunOffLogFileID
        logical                             :: IsFinalFile
        integer, dimension(2)               :: MaxWC_IJ
        real                                :: X, Y

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunOffID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then


            nUsers = DeassociateInstance(mRUNOFF_,  Me%InstanceID)

            if (nUsers == 0) then

                !Writes file with final condition
                IsFinalFile = .true.
                if (Me%OutPut%RestartFormat == BIN_) then
                    call WriteFinalFile_Bin(IsFinalFile)
                else if (Me%OutPut%RestartFormat == HDF_) then
                    call WriteFinalFile_Hdf(IsFinalFile)
                endif

                call UnitsManager(RunOffLogFileID, OPEN_FILE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR00'

                open(UNIT = RunOffLogFileID, FILE = trim(Me%Files%RunOffLogFile), STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR01'

                write(RunOffLogFileID, *)"------------------------ MOHID LAND-----------------------"
                write(RunOffLogFileID, *)"---------------------- RUNOFF STATS ----------------------"
                write(RunOffLogFileID, *)

                write(RunOffLogFileID, *)"VERSION                        : 5"
                write(RunOffLogFileID, *)"INITIAL_TOTAL_VOLUME           : ", Me%InitialTotalVolume
                write(RunOffLogFileID, *)"FINAL_TOTAL_VOLUME             : ", Me%TotalStoredVolume
                write(RunOffLogFileID, *)"TOTAL_INFLOW_VOLUME            : ", Me%TotalDischargeFlowVolume
                write(RunOffLogFileID, *)"TOTAL_BOUNDARY_VOLUME          : ", Me%TotalBoundaryFlowVolume
                write(RunOffLogFileID, *)"TOTAL_BOUNDARY_INFLOW_VOLUME   : ", Me%TotalBoundaryInflowVolume
                write(RunOffLogFileID, *)"TOTAL_BOUNDARY_OUTFLOW_VOLUME  : ", Me%TotalBoundaryOutflowVolume
                write(RunOffLogFileID, *)"TOTAL_STORMWATER_VOLUME        : ", Me%TotalStormWaterVolume
                write(RunOffLogFileID, *)"TOTAL_RAINFALL_VOLUME          : ", Me%TotalRainfallVolume
                write(RunOffLogFileID, *)"TOTAL_INFILTRATION_VOLUME      : ", Me%TotalInfiltrationVolume
                write(RunOffLogFileID, *)"VOLUME_ERROR                   : ", Me%TotalStoredVolume -        &
                                                                                Me%InitialTotalVolume -       &
                                                                                Me%TotalDischargeFlowVolume - &
                                                                                Me%TotalBoundaryFlowVolume  - & 
                                                                                Me%TotalRainfallVolume      - & 
                                                                                Me%TotalStormWaterVolume    + & 
                                                                                Me%TotalInfiltrationVolume

                write(RunOffLogFileID, *)"AVERAGE_CUMULATIVE_INFILTRATION_IN_METERS            : ",         &
                                                                        Me%AvrgAccInfiltrationDepth
                write(RunOffLogFileID, *)"AVERAGE_CUMULATIVE_INFILTRATION_VOLUME_IN_M3         : ",         &
                                                                        Me%AvrgAccInfiltrationVolume

                write(RunOffLogFileID, *)"TOTAL_MANHOLES_VOLUME          : ", Me%TotalManholesVolume
                write(RunOffLogFileID, *)"TOTAL_INLETS_VOLUME            : ", Me%TotalInletsVolume
                write(RunOffLogFileID, *)"TOTAL_OUTFALLS_VOLUME          : ", Me%TotalOutfallsVolume
                write(RunOffLogFileID, *)"TOTAL_PONDS_VOLUME             : ", Me%TotalPondsVolume
                write(RunOffLogFileID, *)"TOTAL_OPENCHANNELS_VOLUME      : ", Me%TotalOpenChannelVolume
                write(RunOffLogFileID, *)"TOTAL_HEADWALLS_VOLUME         : ", Me%TotalHeadwallsVolume
                
                write(RunOffLogFileID, *)"MAX_1D_VOLUME                  : ", Me%MaxTotal1DVolume
                write(RunOffLogFileID, *)"MAX_2D_VOLUME                  : ", Me%MaxTotal2DVolume
                write(RunOffLogFileID, *)"MAX_1D2D_VOLUME                : ", Me%MaxTotal1D2DVolume
                
                write(RunOffLogFileID, *)"TIME_OF_MAX_1D_VOLUME          : ", Me%TimeOfMaxTotal1DVolume
                write(RunOffLogFileID, *)"TIME_OF_MAX_2D_VOLUME          : ", Me%TimeOfMaxTotal2DVolume
                write(RunOffLogFileID, *)"TIME_OF_MAX_1D2D_VOLUME        : ", Me%TimeOfMaxTotal1D2DVolume


                call WriteGridData  (Me%Files%MassErrorFile,                   &
                     COMENT1          = "MassErrorFile",                       &
                     COMENT2          = "MassErrorFile",                       &
                     HorizontalGridID = Me%ObjHorizontalGrid,                  &
                     FillValue        = -99.0,                                 &
                     OverWrite        = .true.,                                &
                     GridData2D_Real  = Me%MassError,                          &
                     STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR020'

                deallocate(Me%MassError)

        
                if(Me%Output%WriteMaxFlowModulus) then
                    call WriteGridData  (Me%Output%MaxFlowModulusFile,         &
                         COMENT1          = "MaxFlowModulusFile",              &
                         COMENT2          = "MaxFlowModulusFile",              &
                         HorizontalGridID = Me%ObjHorizontalGrid,              &
                         FillValue        = -99.0,                             &
                         OverWrite        = .true.,                            &
                         GridData2D_Real  = Me%Output%MaxFlowModulus,          &
                         STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR030'

                    deallocate(Me%Output%MaxFlowModulus)

                endif
                
                if (Me%Output%WriteMaxWaterColumn) then
                    
                    call WriteGridData  (Me%Output%MaxWaterColumnFile,         &
                         COMENT1          = "MaxWaterColumnFile",              &
                         COMENT2          = "MaxWaterColumnFile",              &
                         HorizontalGridID = Me%ObjHorizontalGrid,              &
                         FillValue        = -99.0,                             &
                         OverWrite        = .true.,                            &
                         GridData2D_Real  = Me%Output%MaxWaterColumn,          &
                         STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR040'

                    if(Me%Output%OutputFloodRisk)then
                        write(RunOffLogFileID, *)"MAX_WATER_COLUMN_IN_METERS     : ", maxval(Me%Output%MaxWaterColumn)

                        MaxWC_IJ = maxloc(Me%Output%MaxWaterColumn)

                        call GetCellZ_XY(Me%ObjHorizontalGrid, MaxWC_IJ(1)-1, MaxWC_IJ(2)-1, 0.5, 0.5, X, Y, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR049'

                        write(RunOffLogFileID, *)"MAX_WATER_COLUMN_I             : ", MaxWC_IJ(1)-1
                        write(RunOffLogFileID, *)"MAX_WATER_COLUMN_J             : ", MaxWC_IJ(2)-1
                        write(RunOffLogFileID, *)"MAX_WATER_COLUMN_X_IN_METERS   : ", X
                        write(RunOffLogFileID, *)"MAX_WATER_COLUMN_Y_IN_METERS   : ", Y

                    endif

                    call GetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR041'
                    
                    !Temporarilu use the max water column array to compute max water level
                    Me%Output%MaxWaterColumn = Me%Output%MaxWaterColumn + Me%ExtVar%Topography
                    
                    !Gets a pointer to Topography
                    call UnGetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR042'

                    call WriteGridData  (Me%Output%MaxWaterLevelFile,          &
                         COMENT1          = "MaxWaterLevelFile",               &
                         COMENT2          = "MaxWaterLevelFile",               &
                         HorizontalGridID = Me%ObjHorizontalGrid,              &
                         FillValue        = -99.0,                             &
                         OverWrite        = .true.,                            &
                         GridData2D_Real  = Me%Output%MaxWaterColumn,          &
                         STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR043'

                    if(Me%Output%OutputFloodRisk)then
                        write(RunOffLogFileID, *)"MAX_WATER_LEVEL_IN_METERS      : ", maxval(Me%Output%MaxWaterColumn)

                        MaxWC_IJ = maxloc(Me%Output%MaxWaterColumn)

                        call GetCellZ_XY(Me%ObjHorizontalGrid, MaxWC_IJ(1)-1, MaxWC_IJ(2)-1, 0.5, 0.5, X, Y, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR043a'

                        write(RunOffLogFileID, *)"MAX_WATER_LEVEL_I              : ", MaxWC_IJ(1)-1
                        write(RunOffLogFileID, *)"MAX_WATER_LEVEL_J              : ", MaxWC_IJ(2)-1
                        write(RunOffLogFileID, *)"MAX_WATER_LEVEL_X_IN_METERS    : ", X
                        write(RunOffLogFileID, *)"MAX_WATER_LEVEL_Y_IN_METERS    : ", Y

                    endif

                    deallocate(Me%Output%MaxWaterColumn)

                    call WriteGridData  (Me%Output%TimeOfMaxWaterColumnFile,   &
                         COMENT1          = "TimeOfMaxWaterColumnFile",        &
                         COMENT2          = "TimeOfMaxWaterColumnFile",        &
                         HorizontalGridID = Me%ObjHorizontalGrid,              &
                         FillValue        = -99.0,                             &
                         OverWrite        = .true.,                            &
                         GridData2D_Real  = Me%Output%TimeOfMaxWaterColumn,    &
                         STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR044'

                    if(Me%Output%OutputFloodRisk)then
                        write(RunOffLogFileID, *)"TIME_OF_MAX_WATER_COLUMN       : ", maxval(Me%Output%TimeOfMaxWaterColumn)
                    endif


                    deallocate(Me%Output%TimeOfMaxWaterColumn)

                    
                    if (Me%Output%WriteVelocityAtMaxWaterColumn) then
                        call WriteGridData  (Me%Output%VelocityAtMaxWaterColumnFile,&
                             COMENT1          = "VelocityAtMaxWaterColumnFile",     &
                             COMENT2          = "VelocityAtMaxWaterColumnFile",     &
                             HorizontalGridID = Me%ObjHorizontalGrid,               &
                             FillValue        = -99.0,                              &
                             OverWrite        = .true.,                             &
                             GridData2D_Real  = Me%Output%VelocityAtMaxWaterColumn, &
                             STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR050'    

                        if(Me%Output%OutputFloodRisk)then
                            write(RunOffLogFileID, *)"VEL_AT_MAX_WATER_COLUMN        : ", maxval(Me%Output%VelocityAtMaxWaterColumn)
                        endif

                        deallocate(Me%Output%VelocityAtMaxWaterColumn)
                    endif

                endif

                if (Me%Output%WriteMaxFloodRisk) then
                    call WriteGridData  (Me%Output%MaxFloodRiskFile,           &
                         COMENT1          = "MaxFloodRisk",                    &
                         COMENT2          = "MaxFloodRisk",                    &
                         HorizontalGridID = Me%ObjHorizontalGrid,              &
                         FillValue        = -99.0,                             &
                         OverWrite        = .true.,                            &
                         GridData2D_Real  = Me%Output%MaxFloodRisk,            &
                         STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR060'

                    deallocate(Me%Output%MaxFloodRisk)

                endif
                

                if (Me%Output%WriteFloodPeriod) then

                    if(Me%Output%nFloodPeriodLimits > 0)then

                        call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR061'

                        do n = 1, Me%Output%nFloodPeriodLimits

                            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                                if (Me%ExtVar%BasinPoints(i, j) == 1) then
                                    Me%Output%FloodPeriod(i,j) = Me%Output%FloodPeriods(i,j,n)
                                endif
                            enddo
                            enddo

                            call WriteGridData  (Me%Output%FloodPeriodFiles(n),     &
                                 COMENT1          = "FloodPeriod",                  &
                                 COMENT2          = "FloodPeriod",                  &
                                 HorizontalGridID = Me%ObjHorizontalGrid,           &
                                 FillValue        = -99.0,                          &
                                 OverWrite        = .true.,                         &
                                 GridData2D_Real  = Me%Output%FloodPeriod,          &
                                 STAT             = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR070'


                        enddo

                        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR071'

                        deallocate(Me%Output%FloodPeriods)

                    else

                        call WriteGridData  (Me%Output%FloodPeriodFile,            &
                             COMENT1          = "FloodPeriod",                     &
                             COMENT2          = "FloodPeriod",                     &
                             HorizontalGridID = Me%ObjHorizontalGrid,              &
                             FillValue        = -99.0,                             &
                             OverWrite        = .true.,                            &
                             GridData2D_Real  = Me%Output%FloodPeriod,             &
                             STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR073'

                    endif

                    deallocate(Me%Output%FloodPeriod)

                endif

                if (Me%Output%WriteFloodArrivalTime) then
                    call WriteGridData  (Me%Output%FloodArrivalTimeFile,       &
                         COMENT1          = "FloodArrivalTime",                &
                         COMENT2          = "FloodArrivalTime",                &
                         HorizontalGridID = Me%ObjHorizontalGrid,              &
                         FillValue        = -99.0,                             &
                         OverWrite        = .true.,                            &
                         GridData2D_Real  = Me%Output%FloodArrivalTime,        &
                         STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR075'

                    deallocate(Me%Output%FloodArrivalTime)

                    write(RunOffLogFileID, *)"TOTAL_FLOODED_AREA_AT_END        : ", Me%Output%TotalFloodedArea
                    write(RunOffLogFileID, *)"MAXIMUM_TOTAL_FLOODED_AREA       : ", Me%Output%MaxTotalFloodedArea
                    write(RunOffLogFileID, *)"MAXIMUM_TOTAL_FLOODED_AREA_TIME  : ", Me%Output%TimeOfMaxTotalFloodedArea

                endif

                write(RunOffLogFileID, *)
                write(RunOffLogFileID, *)"---------------------- RUNOFF STATS ----------------------"
                call UnitsManager(RunOffLogFileID, CLOSE_FILE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR76'


                if (Me%ObjDrainageNetwork /= 0) then
 
!                    if(Me%WriteMaxWaterColumn) call WriteChannelsLevelData

                    nUsers = DeassociateInstance (mDRAINAGENETWORK_, Me%ObjDrainageNetwork)
                    if (nUsers == 0) stop 'KillRunOff - RunOff - ERR080'
                endif
                
                if(Me%StormWaterModel)then
#ifdef _SEWERGEMSENGINECOUPLER_

                    STAT_CALL = SewerGEMSEngine_end()
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR081'
                    
                    STAT_CALL = SewerGEMSEngine_close()
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR082'

#ifdef _CONVERT_SWMM_TO_HDF5_
                    !Convert SewerGEMS SWMM outputs to HDF5 and time series format
                    call ConvertSewerGEMSEngine()
#endif _CONVERT_SWMM_TO_HDF5_


                    if(Me%NumberOfHeadwalls > 0)then
                        do n = 1, Me%NumberOfHeadwalls
                            if(Me%Headwalls(n)%OutputResults)then
                                call WriteDataLine(Me%Headwalls(n)%OutputUnit, '<EndTimeSerie>')
                                
                                call UnitsManager(Me%Headwalls(n)%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
                                if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR0831'
                            endif
                        enddo
                    endif

                    if(Me%NumberOfInlets > 0)then
                        do n = 1, Me%NumberOfInlets
                            if(Me%Inlets(n)%OutputResults)then
                                call WriteDataLine(Me%Inlets(n)%OutputUnit, '<EndTimeSerie>')
                                
                                call UnitsManager(Me%Inlets(n)%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
                                if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR083'
                            endif
                        enddo
                    endif

                    if(Me%NumberOfCrossSections > 0)deallocate(Me%CrossSections)
                    if(Me%NumberOfInlets        > 0)deallocate(Me%Inlets       )
                    if(Me%NumberOfManholes      > 0)deallocate(Me%Manholes     )
                    if(Me%NumberOfOutfalls      > 0)deallocate(Me%Outfalls     )
                    if(Me%NumberOfPonds         > 0)deallocate(Me%Ponds        )
                    if(Me%NumberOfIgnoredNodes  > 0)deallocate(Me%IgnoredNodes )
                    if(Me%NumberOfHeadwalls     > 0)deallocate(Me%Headwalls    )


#endif _SEWERGEMSENGINECOUPLER_
                endif
                

                if (Me%OutPut%Yes) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunOff - ERR090'
                endif
                
                if(Me%OutPut%TimeSeries) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunOff - ERR091'
                endif
                
                if (Me%Discharges) then
                    call Kill_Discharges(Me%ObjDischarges, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunOff - ERR100'
                    
                    if (Me%OutPut%TimeSerieDischON) then
                        do dis = 1, Me%OutPut%DischargesNumber
                            
                            call KillTimeSerie(TimeSerieID         = Me%OutPut%TimeSerieDischID(dis), &
                                                 STAT              = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunOff - ERR105'
                            
                        enddo                    
                        
                        deallocate(Me%OutPut%TimeSerieDischProp)
                        deallocate(Me%OutPut%TimeSerieDischID)                    
                        
                    endif                     
                endif

                if (Me%ImposeBoundaryValue) then

                    deallocate(Me%WaterLevelBoundaryValue)

                    if(Me%BoundaryImposedLevelInTime)then

                        if(Me%HasBoundaryLines)then

                            do n = 1, Me%NumberOfBoundaryLines
                                deallocate(Me%BoundaryLines(n)%I)
                                deallocate(Me%BoundaryLines(n)%J)

                                if(Me%BoundaryLines(n)%Variable)then
                                    call KillTimeSerie(Me%BoundaryLines(n)%TimeSerie%ObjTimeSerie, STAT = STAT_CALL)
                                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunoff - 109' 
                                endif
                            enddo

                        else
                            if (Me%ImposedLevelTS%TimeSerie%ObjTimeSerie /= 0) then
                                call KillTimeSerie(Me%ImposedLevelTS%TimeSerie%ObjTimeSerie, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunoff - ERR110' 
                            endif
                        endif
                    end if

                endif

                if (Me%Output%BoxFluxes) then
                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillRunOff - RunOff - ERR120'
                endif
                
                if (Me%ExtVar%Distortion) then

                    call UnGetHorizontalGrid(Me%ObjHorizontalGrid,           &
                                             Me%ExtVar%RotationX,            &
                                             STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'KillRunOff - RunOff - ERR121'


                    call UnGetHorizontalGrid(Me%ObjHorizontalGrid,           &
                                             Me%ExtVar%RotationY,            &
                                             STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'KillRunOff - RunOff - ERR122'

                endif
                
                !Deassociates External Instances
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR130'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR140'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjGridData)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR150'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR160'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR170'
                
                deallocate(Me%myWaterColumnOld)
                
                deallocate (Me%iFlowX)
                deallocate (Me%iFlowY)
                deallocate (Me%lFlowX)
                deallocate (Me%lFlowY)
                deallocate (Me%iFlowToChannels)
                deallocate (Me%lFlowToChannels)
                deallocate (Me%lFlowBoundary)
                deallocate (Me%iFlowBoundary)
                deallocate (Me%iFlowRouteDFour)

                deallocate (Me%BoundaryCells)

                nullify    (Me%iFlowX)
                nullify    (Me%iFlowY)
                nullify    (Me%lFlowX)
                nullify    (Me%lFlowY)
                nullify    (Me%iFlowToChannels)
                nullify    (Me%lFlowToChannels)
                nullify    (Me%lFlowBoundary)
                nullify    (Me%iFlowBoundary)
                nullify    (Me%iFlowRouteDFour)


                !Deallocates Instance
                call DeallocateInstance ()

                RunOffID   = 0
                STAT_      = SUCCESS_

            end if


        end if cd1


        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillRunOff

    !------------------------------------------------------------------------
#ifdef _CONVERT_SWMM_TO_HDF5_
    subroutine ConvertSewerGEMSEngine

        character(len = :, kind = c_char), allocatable  :: swmmInputFile, swmmBinaryFile, swmmHDF5File
        character(len = :, kind = c_char), allocatable  :: swmmTimeSeriesLocation, timeSeriesDir
    

        swmmInputFile           = trim(ADJUSTL(Me%Files%SWMMinp))//C_NULL_CHAR
        swmmBinaryFile          = trim(ADJUSTL(Me%Files%SWMMout))//C_NULL_CHAR
        swmmHDF5File            = trim(ADJUSTL(Me%Files%SWMMHDF))//'5'//C_NULL_CHAR
        timeSeriesDir           = trim(ADJUSTL(Me%Files%SWMMTimeSeriesDir))//C_NULL_CHAR
        swmmTimeSeriesLocation  = trim(ADJUSTL(Me%Files%SWMMTimeSeries))//C_NULL_CHAR
        call ConvertSewerGEMSEngineToDrainageNetwork(swmmInputFile, swmmBinaryFile, swmmHDF5File, &
                                                     timeSeriesDir, swmmTimeSeriesLocation)
        
     
    end subroutine ConvertSewerGEMSEngine
#endif _CONVERT_SWMM_TO_HDF5_
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_RunOff), pointer                    :: AuxObjRunOff
        type (T_RunOff), pointer                    :: PreviousObjRunOff

        !Updates pointers
        if (Me%InstanceID == FirstObjRunOff%InstanceID) then
            FirstObjRunOff => FirstObjRunOff%Next
        else
            PreviousObjRunOff => FirstObjRunOff
            AuxObjRunOff      => FirstObjRunOff%Next
            do while (AuxObjRunOff%InstanceID /= Me%InstanceID)
                PreviousObjRunOff => AuxObjRunOff
                AuxObjRunOff      => AuxObjRunOff%Next
            enddo

            !Now update linked list
            PreviousObjRunOff%Next => AuxObjRunOff%Next

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

    subroutine Ready (RunOffID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: RunOffID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (RunOffID > 0) then
            call LocateObjRunOff (RunOffID)
            ready_ = VerifyReadLock (mRUNOFF_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjRunOff (ObjRunOffID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjRunOffID

        !Local-----------------------------------------------------------------

        Me => FirstObjRunOff
        do while (associated (Me))
            if (Me%InstanceID == ObjRunOffID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleRunOff - LocateObjRunOff - ERR01'

    end subroutine LocateObjRunOff

    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar (StaticOnly)
        
        !Arguments-------------------------------------------------------------
        logical                                     :: StaticOnly

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Time Stuff
        call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR01'

        call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR02'

        !Gets Basin Points
        call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR03'
        
        !Gets cell slope
        call GetCellSlope   (Me%ObjBasinGeometry, Me%ExtVar%CellSlope, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR04'

        !Gets River Points
        call GetRiverPoints (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR05'

        !Gets Horizontal Grid
        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                     &
                               DUX    = Me%ExtVar%DUX,    DVY    = Me%ExtVar%DVY,        &
                               DUY    = Me%ExtVar%DUY,    DVX    = Me%ExtVar%DVX,        &
                               DXX    = Me%ExtVar%DXX,    DYY    = Me%ExtVar%DYY,        &
                               DZX    = Me%ExtVar%DZX,    DZY    = Me%ExtVar%DZY,        &
                               XX2D_Z = Me%ExtVar%XX2D_Z, YY2D_Z = Me%ExtVar%YY2D_Z,     &
                               STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR06'

        call GetGridCellArea  (Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea,             &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR06a'

        !Gets a pointer to Topography
        call GetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR07'

        if (.not. StaticOnly) then

            !Gets Boundary Points
            call GetBoundaries    (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR10'
        
        endif

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar(StaticOnly)
        
        !Arguments-------------------------------------------------------------
        logical                                     :: StaticOnly
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Unget Basin Points
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR01'

        !Unget River Points
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR02'

        !Unget Cell Slope
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%CellSlope, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR02a'

        !Unget Horizontal Grid
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR03'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR04'
        
        !Unget Horizontal Grid
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR03a'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR04a'
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DXX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR05'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DYY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR06'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR05'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR06'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%XX2D_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR07'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%YY2D_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR08'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR09'

        !Ungets the Topography
        call UngetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR10'

        if (.not. StaticOnly) then

            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR11'

        endif 
        
    end subroutine ReadUnLockExternalVar
    
end module ModuleRunOff
