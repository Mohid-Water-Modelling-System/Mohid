    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Land
    ! MODULE        : SWMM Coupler
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : October 2019
    ! REVISION      : Ricardo Canelas
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION
    !> Module that provides a specific interface to couple send and request data to
    !> a SWMM model.
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

    module SWMMCoupler

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleHorizontalGrid
    use iso_c_binding

    implicit none
    private

    interface

    subroutine swmm_open(inFile, rptFile, outFile) bind(C, name='swmm_open')
    use iso_c_binding
    character(kind = c_char) :: inFile(*)
    character(kind = c_char) :: rptFile(*)
    character(kind = c_char) :: outFile(*)
    end subroutine swmm_open

    subroutine swmm_start(saveResults) bind(C, name='swmm_start')
    use iso_c_binding
    integer(c_int) :: saveResults
    end subroutine swmm_start

    subroutine swmm_end() bind(C, name='swmm_end')
    end subroutine swmm_end

    subroutine swmm_close() bind(C, name='swmm_close')
    end subroutine swmm_close

    end interface

    interface
    subroutine swmm_getNumberOfNodes(nNodes) bind(C, name='swmm_getNumberOfNodes')
    use iso_c_binding
    integer(c_int) :: nNodes
    end subroutine swmm_getNumberOfNodes
    end interface

    interface
    subroutine swmm_getNodeType(id, nType) bind(C, name='swmm_getNodeTypeByID')
    use iso_c_binding
    integer(c_int) :: id
    integer(c_int) :: nType
    end subroutine swmm_getNodeType
    end interface

    interface
    subroutine swmm_getNodeXY(id, xx, yy) bind(C, name='swmm_getNodeXYByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: xx
    real(c_double) :: yy
    end subroutine swmm_getNodeXY
    end interface

    interface
    subroutine swmm_getIsNodeOpenChannel(id, isOpen) bind(C, name='swmm_getIsNodeOpenChannelByID')
    use iso_c_binding
    integer(c_int) :: id
    integer(c_int) :: isOpen
    end subroutine swmm_getIsNodeOpenChannel
    end interface

    interface
    subroutine swmm_getNodeHasLateralInflow(id, hasLatFlow) bind(C, name='swmm_getNodeHasLateralInflowByID')
    use iso_c_binding
    integer(c_int) :: id
    integer(c_int) :: hasLatFlow
    end subroutine swmm_getNodeHasLateralInflow
    end interface

    interface
    subroutine swmm_getInflowByNode(id, inflow) bind(C, name='swmm_getInflowByNodeByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: inflow
    end subroutine swmm_getInflowByNode
    end interface

    interface
    subroutine swmm_getOutflowByNode(id, outflow) bind(C, name='swmm_getOutflowByNodeByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: outflow
    end subroutine swmm_getOutflowByNode
    end interface

    interface
    subroutine swmm_getLevelByNode(id, level) bind(C, name='swmm_getLevelByNodeByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: level
    end subroutine swmm_getLevelByNode
    end interface

    interface
    subroutine swmm_setDownstreamWaterLevel(id, level) bind(C, name='swmm_setDownstreamWaterLevelByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: level
    end subroutine swmm_setDownstreamWaterLevel
    end interface

    interface
    subroutine swmm_setLateralInflow(id, inflow) bind(C, name='swmm_setLateralInflowByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: inflow
    end subroutine swmm_setLateralInflow
    end interface

    interface
    subroutine swmm_setPondedWaterColumn(id, level) bind(C, name='swmm_setPondedWaterColumnByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: level
    end subroutine swmm_setPondedWaterColumn
    end interface

    interface
    subroutine swmm_setStormWaterPotentialInflow(id, inflow) bind(C, name='swmm_setStormWaterPotentialInflowByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: inflow
    end subroutine swmm_setStormWaterPotentialInflow
    end interface

    interface
    subroutine swmm_setOpenXSectionInflow(id, inflow) bind(C, name='swmm_setOpenXSectionInflowByID')
    use iso_c_binding
    integer(c_int) :: id
    real(c_double) :: inflow
    end subroutine swmm_setOpenXSectionInflow
    end interface

    interface
    subroutine swmm_step(elapsedTime) bind(C, name='swmm_step')
    use iso_c_binding
    real(c_double) :: elapsedTime
    end subroutine swmm_step
    end interface

    interface
    subroutine swmm_step_imposed_dt(elapsedTime, imposedDt) bind(C, name='swmm_step_imposed_dt')
    use iso_c_binding
    real(c_double) :: elapsedTime
    real(c_double) :: imposedDt
    end subroutine swmm_step_imposed_dt
    end interface

    interface
    subroutine swmm_getdt(Dt) bind(C, name='swmm_getdt')
    use iso_c_binding
    real(c_double) :: Dt
    end subroutine swmm_getdt
    end interface

    interface
    subroutine ConvertSwmmToDrainageNetwork(f1, f2, f3, f4, f5) bind(C, name='ConvertSwmmToDrainageNetworkCaller')
    use iso_c_binding
    character(kind = c_char) :: f1(*)
    character(kind = c_char) :: f2(*)
    character(kind = c_char) :: f3(*)
    character(kind = c_char) :: f4(*)
    character(kind = c_char) :: f5(*)
    end subroutine ConvertSwmmToDrainageNetwork
    end interface


    type :: NodeTypes_enum          !< enums for node types
        integer :: JUNCTION = 0
        integer :: OUTFALL = 1
    end type NodeTypes_enum

    integer, parameter :: nodeId = 1
    integer, parameter :: cellID = 2
    integer, parameter :: cellI = 3
    integer, parameter :: cellJ = 4
    integer, parameter :: isJunction = 5
    integer, parameter :: isOutfall = 6

    !main public class
    type :: swmm_coupler_class                      !< SWMM Coupler class
        logical :: initialized = .false.                        !< initialized flag
        type(NodeTypes_enum) :: NodeTypes                       !< node type flags
        integer :: NumberOfNodes                                !< number of SWMM nodes
        integer :: NumberOfInDomainNodes                        !< number of SWMM  nodes in the domain
        real, allocatable, dimension(:,:)  :: nodeXY            !< position of the SWMM nodes
        integer, allocatable, dimension(:,:)  :: nodeIJ         !< position of the SWMM nodes in mesh coordinates
        integer, allocatable, dimension(:,:)  :: n2cMap         !< node to cell mappings
        integer, allocatable, dimension(:) :: junctionIDX       !< ids of junction SWMM nodes
        integer, allocatable, dimension(:) :: outfallIDX        !< ids of outfall SWMM nodes
        integer, allocatable, dimension(:) :: inflowIDX         !< ids of inflow SWMM nodes
        integer, allocatable, dimension(:) :: xsectionLevelsIDX !< ids of open cross section SWMM nodes
        integer, allocatable, dimension(:) :: lateralFlowIDX    !< ids of lateral flow accepting SWMM nodes
        logical, allocatable, dimension(:) :: xSectionOpen      !warning: 1-based index array (c equivalent is 0-based)
        character(len = StringLength)      :: SWMM_dat
        character(len = StringLength)      :: SWMM_rpt
        character(len = StringLength)      :: SWMM_out
        character(len = StringLength)      :: STORMWATER_HDF
        character(len = StringLength)      :: SWMM_timeSeries_location
        character(len = StringLength)      :: SWMM_timeSeries_dir
    contains
    procedure :: initialize => initSWMMCoupler
    procedure :: mapElements
    procedure :: finalize => finalizeSWMMCoupler
    procedure :: print => printSWMMCoupler
    !mapping procedures
    procedure, private :: inDomainNode
    !control procedures
    procedure, private :: initializeSWMM
    procedure, private :: getSWMMFilesPaths
    procedure :: defaultPerformTimeStep
    procedure :: runStep => PerformTimeStep
    procedure, private :: finalizeSWMM
    procedure, private :: convertSWMM
    !import data procedures
    procedure :: GetDt
    procedure :: GetInflow
    procedure :: GetOutflow
    procedure :: GetLevel
    procedure :: GetCellList
    procedure, private :: GetNumberOfNodes
    procedure, private :: GetNodeXY
    procedure, private :: GetNodeXYByID
    procedure, private :: GetNodeTypeByID
    procedure, private :: GetIsNodeOpenChannel
    procedure, private :: GetNodeHasLateralInflow
    procedure, private :: GetInflowByID
    procedure, private :: GetOutflowByID
    procedure, private :: GetLevelByID
    !export data procedures
    procedure :: SetOutletLevel
    procedure :: SetLateralInflow
    procedure :: SetWaterColumn
    procedure :: SetInletInflow
    procedure :: SetXSectionInflow
    procedure, private :: SetOutletLevelByID
    procedure, private :: SetLateralInflowByID
    procedure, private :: SetWaterColumnByID
    procedure, private :: SetInletInflowByID
    procedure, private :: SetXSectionInflowByID
    end type swmm_coupler_class


    !Public access vars
    public :: swmm_coupler_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Initializes the SWMM Coupler object
    !---------------------------------------------------------------------------
    subroutine initSWMMCoupler(self, mapArrayXY, mapArrayIJ, mapArrayID)
    class(swmm_coupler_class), intent(inout) :: self
    real, allocatable, dimension(:,:), intent(inout) :: mapArrayXY
    integer, allocatable, dimension(:,:), intent(inout) :: mapArrayIJ
    integer, allocatable, dimension(:), intent(inout) :: mapArrayID

    print*, 'Initializing SWMM coupler, please wait...'

    self%initialized = .true.
    call self%getSWMMFilesPaths()
    call self%initializeSWMM()
    call self%GetNumberOfNodes()
    call self%GetNodeXY()

    print*, 'SWMM number of nodes is ', self%NumberOfNodes

    !allocating map arrays to send to Basin Module to get filled/used
    allocate(mapArrayXY(self%NumberOfNodes,2))
    allocate(mapArrayIJ(self%NumberOfNodes,2))
    allocate(mapArrayID(self%NumberOfNodes))
    mapArrayXY = self%nodeXY

    end subroutine initSWMMCoupler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Maps the SWMM Coupler object elements
    !---------------------------------------------------------------------------
    subroutine mapElements(self, mapArrayIJ, mapArrayID)
    class(swmm_coupler_class), intent(inout) :: self
    integer, dimension(:,:), intent(inout) :: mapArrayIJ
    integer, dimension(:), intent(inout) :: mapArrayID
    integer :: nJunction = 0
    integer :: nOutfall = 0
    integer :: nInflow = 0
    integer :: nXSection = 0
    integer :: nLatFlow = 0
    integer :: i, idx
    integer :: idxj, idxo, idxi, idxx, idxl

    print*, 'Mapping coupling points, please wait...'

    !do i=1, self%NumberOfNodes
    !    print*, 'id=',i, 'x=',self%nodeXY(i,1), 'y=',self%nodeXY(i,2)
    !    print*, 'cell id=', mapArrayID(i), 'i=',mapArrayIJ(i,1), 'j=',mapArrayIJ(i,2)
    !end do

    allocate(self%n2cMap(self%NumberOfNodes,6))
    self%n2cMap = 0

    do i=1, self%NumberOfNodes
        self%n2cMap(i,nodeId) = i
        self%n2cMap(i,cellID) = mapArrayID(i)
        self%n2cMap(i,cellI) = mapArrayIJ(i,1)
        self%n2cMap(i,cellJ) = mapArrayIJ(i,2)
        if (mapArrayIJ(i,1) == null_int .or. mapArrayIJ(i,2) == null_int) self%n2cMap(i,2) = null_int !nodes outside of the domain
    end do

    self%NumberOfInDomainNodes = count(self%n2cMap(:,cellID) /= null_int)
    print*, 'SWMM number of nodes in domain is ', self%NumberOfInDomainNodes

    !do i=1, self%NumberOfNodes
    !    if (self%inDomainNode(i)) print*, self%n2cMap(i,:)
    !end do

    !building open section list
    allocate(self%xSectionOpen(self%NumberOfNodes))
    self%xSectionOpen = .false.
    do i=1, self%NumberOfNodes
        if (self%NodeTypes%junction == self%GetNodeTypeByID(i-1)) then
            if (self%GetIsNodeOpenChannel(i-1)) then
                if (self%inDomainNode(i)) self%xSectionOpen(i) = .true.
            end if
        end if
    end do

    !building id lists for O(1) access
    do i=1, self%NumberOfNodes
        idx = i-1
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            if (self%inDomainNode(i)) nJunction = nJunction + 1
        end if
        if (self%NodeTypes%outfall  == self%GetNodeTypeByID(idx)) then
            if (self%inDomainNode(i)) nOutfall  = nOutfall  + 1
        end if
        if (self%NodeTypes%junction  == self%GetNodeTypeByID(idx)) then
            if (self%inDomainNode(i)) then
                if (.not.self%xSectionOpen(i)) nInflow = nInflow + 1     !only closed nodes
            end if
        end if
        if (self%NodeTypes%junction  == self%GetNodeTypeByID(idx)) then
            if (self%inDomainNode(i)) then
                if (self%xSectionOpen(i)) nXSection = nXSection + 1      !only open nodes
            end if
        end if
        if (self%NodeTypes%junction  == self%GetNodeTypeByID(idx)) then
            if (self%GetNodeHasLateralInflow(idx)) then
                if (self%inDomainNode(i)) nLatFlow = nLatFlow + 1
            end if
        end if
    end do

    !print*, ""
    !print*, "nJunction", nJunction
    !print*, "nOutfall", nOutfall
    !print*, "nInflow", nInflow
    !print*, "nXSection", nXSection
    !print*, "nLatFlow", nLatFlow
    !print*, ""

    allocate(self%junctionIDX(nJunction))
    allocate(self%outfallIDX(nOutfall))
    allocate(self%inflowIDX(nInflow))
    allocate(self%xsectionLevelsIDX(nXSection))
    allocate(self%lateralFlowIDX(nLatFlow))
    idxj =1
    idxo =1
    idxi =1
    idxx =1
    idxl =1
    do i=1, self%NumberOfNodes
        idx = i-1
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            if (self%inDomainNode(i)) then
                self%junctionIDX(idxj) = idx
                idxj = idxj + 1
                self%n2cMap(i, isJunction) = 1
            end if
        end if
        if (self%NodeTypes%outfall == self%GetNodeTypeByID(idx)) then
            if (self%inDomainNode(i)) then
                self%outfallIDX(idxo) = idx
                idxo = idxo + 1
                self%n2cMap(i, isOutfall) = 1
            end if
        end if
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            if (.not.self%xSectionOpen(i)) then !only closed nodes
                if (self%inDomainNode(i)) then
                    self%inflowIDX(idxi) = idx
                    idxi = idxi + 1
                end if
            end if
        end if
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            if (self%xSectionOpen(i)) then      !only open nodes
                if (self%inDomainNode(i)) then
                    self%xsectionLevelsIDX(idxx) = idx
                    idxx = idxx + 1
                end if
            end if
        end if
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            if (self%GetNodeHasLateralInflow(i)) then
                if (self%inDomainNode(i)) then
                    self%lateralFlowIDX(idxl) = idx
                    idxl = idxl + 1
                end if
            end if
        end if
    end do

    end subroutine mapElements

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Returns true if SWMM node is in the coupled domain
    !---------------------------------------------------------------------------
    logical function inDomainNode(self, idx)
    class(swmm_coupler_class), intent(in) :: self
    integer, intent(in) :: idx
    inDomainNode = self%n2cMap(idx, cellID) /= null_int
    end function inDomainNode

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Finalizes the SWMM Coupler object and the SWMM model
    !---------------------------------------------------------------------------
    subroutine finalizeSWMMCoupler(self)
    class(swmm_coupler_class), intent(inout) :: self

    call self%finalizeSWMM()

    end subroutine finalizeSWMMCoupler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the files required for SWMM run
    !---------------------------------------------------------------------------
    subroutine getSWMMFilesPaths(self)
    class(swmm_coupler_class), intent(inout) :: self
    integer :: STAT_CALL, fileID, iflag, dpos
    character(len=StringLength) :: dummy

    call ReadFileName('STORMWATER_DAT', self%SWMM_dat,                         &
        Message = "SWMM input file", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getSWMMFilesPaths - STORMWATER_DAT keyword not found on main file list'

    call ReadFileName('STORMWATER_HDF', self%STORMWATER_HDF,                         &
        Message = "SWMM hdf5 output file", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getSWMMFilesPaths - STORMWATER_HDF keyword not found on main file list'
    
    dpos = scan(trim(self%STORMWATER_HDF),".", BACK= .true.)
    
    !name of output files is the same as the output hdf5
    if (dpos > 0) then
        self%SWMM_rpt = self%STORMWATER_HDF(1:dpos)//'rpt'
        self%SWMM_out = self%STORMWATER_HDF(1:dpos)//'swm'
    else
        stop 'SWMMCoupler::getSWMMFilesPaths - STORMWATER_HDF file extension not recognized'
    end if
    
    call ReadFileName('ROOT_SRT', self%SWMM_timeSeries_dir,                         &
        Message = "SWMM Time Series output directory", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getSWMMFilesPaths - ROOT_SRT keyword not found on main file list'
    
    call ConstructEnterData (fileID, self%SWMM_dat, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getSWMMFilesPaths - SWMM_dat file not read'

    call GetData(dummy,                                            &
        fileID, iflag,                                             &
        SearchType   = FromFile,                                   &
        keyword      = 'TIME_SERIE_LOCATION',                      &        
        ClientModule = 'SWMMCoupler',                              &
        STAT         = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getSWMMFilesPaths - SWMM_dat file not readable'    
    if (iflag == 0) then
        self%SWMM_timeSeries_location = dummy
    else
        self%SWMM_timeSeries_location = ""
        self%SWMM_timeSeries_dir = ""
    endif
    
    call GetData(dummy,                                            &
        fileID, iflag,                                             &
        SearchType   = FromFile,                                   &
        keyword      = 'MODEL_CONFIGURATION',                      &        
        ClientModule = 'SWMMCoupler',                              &
        STAT         = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getSWMMFilesPaths - SWMM_dat file not readable'
    if (iflag /= 0) then
        self%SWMM_dat = dummy
    endif    
    
    end subroutine getSWMMFilesPaths

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Initializes the SWMM model through a DLL call
    !---------------------------------------------------------------------------
    subroutine initializeSWMM(self)
    class(swmm_coupler_class), intent(in) :: self
    character(len = :, kind = c_char), allocatable :: inFile, rptFile, outFile
    integer(c_int) :: saveResults

    print*, 'Initializing SWMM, please wait...'

    inFile = trim(ADJUSTL(self%SWMM_dat))//C_NULL_CHAR
    rptFile = trim(ADJUSTL(self%SWMM_rpt))//C_NULL_CHAR
    outFile = trim(ADJUSTL(self%SWMM_out))//C_NULL_CHAR
    saveResults = 1

#ifdef _SWMMCoupler_
    call swmm_open(inFile, rptFile, outFile)
    call swmm_start(saveResults)
#endif

    print*,''

    end subroutine initializeSWMM

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> calls the default time step integrator from SWMM model through a DLL call
    !---------------------------------------------------------------------------
    subroutine defaultPerformTimeStep(self)
    class(swmm_coupler_class), intent(in) :: self
    real(c_double) :: elapsedTime

#ifdef _SWMMCoupler_
    call swmm_step(elapsedTime)
#endif

    end subroutine defaultPerformTimeStep

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> calls the time step integrator from SWMM model with an imposed Dt
    !> through a DLL call
    !---------------------------------------------------------------------------
    subroutine PerformTimeStep(self, dt)
    class(swmm_coupler_class), intent(in) :: self
    real(c_double), intent(in) :: dt
    real(c_double) :: elapsedTime

#ifdef _SWMMCoupler_
    call swmm_step_imposed_dt(elapsedTime, dt)
#endif

    end subroutine PerformTimeStep

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the SWMM computed dt through a DLL call
    !---------------------------------------------------------------------------
    real function GetDt(self)
    class(swmm_coupler_class), intent(in) :: self
    real(c_double) :: dt

#ifdef _SWMMCoupler_
    call swmm_getdt(dt)
#endif
    GetDt = dt

    end function GetDt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Finalizes the SWMM model through DLL calls
    !---------------------------------------------------------------------------
    subroutine finalizeSWMM(self)
    class(swmm_coupler_class), intent(in) :: self
    
    print*, 'Finalizing SWMM, please wait...'
#ifdef _SWMMCoupler_
    call swmm_end()
    call swmm_close()
    call self%convertSWMM()
#endif

    end subroutine finalizeSWMM
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Converts the SWMM output through DLL calls
    !---------------------------------------------------------------------------
    subroutine convertSWMM(self)
    class(swmm_coupler_class), intent(in) :: self
    character(len = :, kind = c_char), allocatable :: swmmInputFile, swmmBinaryFile, drainageNetworkFile, timeSeriesDir, swmmTimeSeriesLocation

#ifdef _SWMMCoupler_    
    swmmInputFile = trim(ADJUSTL(self%SWMM_dat))//C_NULL_CHAR
    swmmBinaryFile = trim(ADJUSTL(self%SWMM_out))//C_NULL_CHAR
    drainageNetworkFile = trim(ADJUSTL(self%STORMWATER_HDF))//'5'//C_NULL_CHAR
    timeSeriesDir = trim(ADJUSTL(self%SWMM_timeSeries_dir))//C_NULL_CHAR
    swmmTimeSeriesLocation = trim(ADJUSTL(self%SWMM_timeSeries_location))//C_NULL_CHAR
    call ConvertSwmmToDrainageNetwork(swmmInputFile, swmmBinaryFile, drainageNetworkFile, timeSeriesDir, swmmTimeSeriesLocation)
#endif

    end subroutine convertSWMM

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the number of SWMM nodes through a DLL call
    !---------------------------------------------------------------------------
    subroutine GetNumberOfNodes(self)
    class(swmm_coupler_class), intent(inout) :: self
    integer(c_int) :: nNodes

#ifdef _SWMMCoupler_
    call swmm_getNumberOfNodes(nNodes)
#endif
    self%NumberOfNodes = nNodes

    end subroutine GetNumberOfNodes

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the xy coordinates of all SWMM nodes
    !---------------------------------------------------------------------------
    subroutine GetNodeXY(self)
    class(swmm_coupler_class), intent(inout) :: self
    real, dimension(2) :: xy
    integer :: i

    allocate(self%NodeXY(self%NumberOfNodes,2))
    do i = 1, self%NumberOfNodes
        xy = self%GetNodeXYByID(i)
        self%NodeXY(i,1) = xy(1)
        self%NodeXY(i,2) = xy(2)
    end do

    end subroutine GetNodeXY

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the xy coordinates of a SWMM node through a DLL call
    !---------------------------------------------------------------------------
    function GetNodeXYByID(self, id) result(xy)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double) :: x, y
    real, dimension(2) :: xy

#ifdef _SWMMCoupler_
    call swmm_getNodeXY(id, x, y)
#endif
    xy(1) = x
    xy(2) = y

    end function GetNodeXYByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the SWMM node type through a DLL call
    !> @param[in] self, id
    !---------------------------------------------------------------------------
    integer function GetNodeTypeByID(self, id)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    integer(c_int) :: nType

#ifdef _SWMMCoupler_
    call swmm_getNodeType(id, nType)
#endif
    GetNodeTypeByID = nType

    end function GetNodeTypeByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the open channel status of a SWMM node through a DLL call
    !> @param[in] self, id
    !---------------------------------------------------------------------------
    logical function GetIsNodeOpenChannel(self, id)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    integer(c_int) :: isOpen

    GetIsNodeOpenChannel = .false.
#ifdef _SWMMCoupler_
    call swmm_getIsNodeOpenChannel(id, isOpen)
#endif
    if (isOpen == 1) GetIsNodeOpenChannel = .true.

    end function GetIsNodeOpenChannel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the has lateral flow status of a SWMM node through a DLL call
    !> @param[in] self, id
    !---------------------------------------------------------------------------
    logical function GetNodeHasLateralInflow(self, id)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    integer(c_int) :: isOpen

    GetNodeHasLateralInflow = .false.
#ifdef _SWMMCoupler_
    call swmm_getNodeHasLateralInflow(id, isOpen)
#endif
    if (isOpen == 1) GetNodeHasLateralInflow = .true.

    end function GetNodeHasLateralInflow

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets ponded inflow by SWMM node ID through a DLL call
    !> @param[in] self, id
    !---------------------------------------------------------------------------
    real function GetInflowByID(self, id)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double) :: inflow

#ifdef _SWMMCoupler_
    call swmm_getInflowByNode(id, inflow)
#endif
    GetInflowByID = inflow

    end function GetInflowByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets ponded inflow at all required SWMM nodes
    !---------------------------------------------------------------------------
    function GetInflow(self, HorizontalGridID) result(inflow)
    class(swmm_coupler_class), intent(in) :: self
    integer, intent(in) :: HorizontalGridID
    real, dimension(size(self%inflowIDX), 3) :: inflow
    integer :: i, ii, jj, cid, nnodes

    inflow = 0.0
    if (size(self%inflowIDX)>0) then
        do i=1, size(self%inflowIDX)
            cid = self%n2cMap(self%inflowIDX(i), cellID) !cell id from this node
            nnodes = count(self%n2cMap(:,cellID) == cid)         !number of nodes sharing the same cell
            call GetCellIJfromID(HorizontalGridID, ii, jj, cid)
            inflow(i,1) = ii
            inflow(i,2) = jj
            inflow(i,3) = self%GetInflowByID(self%inflowIDX(i))
        end do
    end if

    end function GetInflow

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets outflow by SWMM node ID through a DLL call
    !> @param[in] self, id
    !---------------------------------------------------------------------------
    real function GetOutflowByID(self, id)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double) :: outflow

#ifdef _SWMMCoupler_
    call swmm_getOutflowByNode(id, outflow)
#endif
    GetOutflowByID = outflow

    end function GetOutflowByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets outflow at all required SWMM nodes
    !---------------------------------------------------------------------------
    function GetOutflow(self, HorizontalGridID) result(outflow)
    class(swmm_coupler_class), intent(in) :: self
    integer, intent(in) :: HorizontalGridID
    real, dimension(size(self%outfallIDX), 3) :: outflow
    integer :: i, ii, jj, cid, nnodes

    outflow = 0.0
    if (size(self%outfallIDX)>0) then
        do i=1, size(self%outfallIDX)
            cid = self%n2cMap(self%outfallIDX(i), cellID) !cell id from this node
            nnodes = count(self%n2cMap(:,cellID) == cid)         !number of nodes sharing the same cell
            call GetCellIJfromID(HorizontalGridID, ii, jj, cid)
            outflow(i,1) = ii
            outflow(i,2) = jj
            outflow(i,3) = self%GetOutflowByID(self%outfallIDX(i))
        end do
    end if

    end function GetOutflow

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets level by SWMM node ID through a DLL call
    !> @param[in] self, id
    !---------------------------------------------------------------------------
    real function GetLevelByID(self, id)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double) :: level

#ifdef _SWMMCoupler_
    call swmm_getLevelByNode(id, level)
#endif
    GetLevelByID = level

    end function GetLevelByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets level at all required SWMM nodes
    !---------------------------------------------------------------------------
    function GetLevel(self, HorizontalGridID) result(level)
    class(swmm_coupler_class), intent(in) :: self
    integer, intent(in) :: HorizontalGridID
    real, dimension(size(self%xsectionLevelsIDX), 3) :: level
    integer :: i, ii, jj, cid, nnodes

    level = 0.0
    if (size(self%xsectionLevelsIDX)>0) then
        do i=1, size(self%xsectionLevelsIDX)
            cid = self%n2cMap(self%xsectionLevelsIDX(i), cellID) !cell id from this node
            nnodes = count(self%n2cMap(:,cellID) == cid)         !number of nodes sharing the same cell
            call GetCellIJfromID(HorizontalGridID, ii, jj, cid)
            level(i,1) = ii
            level(i,2) = jj
            level(i,3) = self%GetLevelByID(self%xsectionLevelsIDX(i))!/nnodes !contribution of this node to the average of this cell - WIP
        end do
    end if

    end function GetLevel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets outlet level at a SWMM node ID through a DLL call
    !> @param[in] self, id, outletLevel
    !---------------------------------------------------------------------------
    subroutine SetOutletLevelByID(self, id, outletLevel)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double), intent(in) :: outletLevel

#ifdef _SWMMCoupler_
    call swmm_setDownstreamWaterLevel(id, outletLevel)
#endif

    end subroutine SetOutletLevelByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets outlet level at all suited SWMM nodes
    !> @param[in] self, outletLevel
    !---------------------------------------------------------------------------
    subroutine SetOutletLevel(self, outletLevel, cellIDs)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: outletLevel
    integer, dimension(:), intent(in) :: cellIDs
    integer :: i, cid, cidx

    if (size(self%outfallIDX)>0) then
        do i=1, size(self%outfallIDX)
            cid = self%n2cMap(self%outfallIDX(i), cellID)         !cell id from this node
            cidx = findloc(cellIDs, value = cid, dim = 1)         !index of the cell in the data array
            if (cidx > 0) call self%SetOutletLevelByID(self%outfallIDX(i), outletLevel(cidx))
        end do
    end if
    end subroutine SetOutletLevel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets outlet level at a SWMM node ID through a DLL call
    !> @param[in] self, id, inflow
    !---------------------------------------------------------------------------
    subroutine SetLateralInflowByID(self, id, inflow)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double), intent(in) :: inflow

#ifdef _SWMMCoupler_
    call swmm_setLateralInflow(id, inflow)
#endif

    end subroutine SetLateralInflowByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets outlet level at all suited SWMM nodes
    !> @param[in] self, inflow
    !---------------------------------------------------------------------------
    subroutine SetLateralInflow(self, inflow, cellIDs)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: inflow
    integer, dimension(:), intent(in) :: cellIDs
    integer :: i, cid, nnodes, cidx

    if (size(self%junctionIDX)>0) then
        do i=1, size(self%junctionIDX)
            cid = self%n2cMap(self%junctionIDX(i), cellID)       !cell id from this node
            nnodes = count(self%n2cMap(:,cellID) == cid)         !number of nodes sharing the same cell
            cidx = findloc(cellIDs, value = cid, dim = 1)        !index of the cell in the data array
            if (cidx > 0) call self%SetLateralInflowByID(self%junctionIDX(i), inflow(cidx)/nnodes)
        end do
    end if

    end subroutine SetLateralInflow

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets water level at a SWMM node ID through a DLL call
    !> @param[in] self, id, level
    !---------------------------------------------------------------------------
    subroutine SetWaterColumnByID(self, id, level)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double), intent(in) :: level

#ifdef _SWMMCoupler_
    call swmm_setPondedWaterColumn(id, level)
#endif

    end subroutine SetWaterColumnByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets outlet level at all suited SWMM nodes
    !> @param[in] self, level
    !---------------------------------------------------------------------------
    subroutine SetWaterColumn(self, level, cellIDs)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: level
    integer, dimension(:), intent(in) :: cellIDs
    integer :: i, cid, cidx

    if (size(self%inflowIDX)>0) then
        do i=1, size(self%inflowIDX)
            cid = self%n2cMap(self%inflowIDX(i), cellID)         !cell id from this node
            cidx = findloc(cellIDs, value = cid, dim = 1)        !index of the cell in the data array
            if (cidx > 0) call self%SetWaterColumnByID(self%inflowIDX(i), level(cidx))
        end do
    end if
    end subroutine SetWaterColumn

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets stormwater inflow at a SWMM node ID through a DLL call
    !> @param[in] self, id, inflow
    !---------------------------------------------------------------------------
    subroutine SetInletInflowByID(self, id, inflow)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double), intent(in) :: inflow

#ifdef _SWMMCoupler_
    call swmm_setStormWaterPotentialInflow(id, inflow)
#endif

    end subroutine SetInletInflowByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets inflow at all suited SWMM nodes
    !> @param[in] self, inflow
    !---------------------------------------------------------------------------
    subroutine SetInletInflow(self, inflow, cellIDs)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: inflow
    integer, dimension(:), intent(in) :: cellIDs
    integer :: i, cid, nnodes, cidx

    if (size(self%inflowIDX)>0) then
        do i=1, size(self%inflowIDX)
            cid = self%n2cMap(self%inflowIDX(i), cellID)         !cell id from this node
            nnodes = count(self%n2cMap(:,cellID) == cid)         !number of nodes sharing the same cell
            cidx = findloc(cellIDs, value = cid, dim = 1)        !index of the cell in the data array
            if (cidx > 0) call self%SetInletInflowByID(self%inflowIDX(i), inflow(cidx)/nnodes)
        end do
    end if
    end subroutine SetInletInflow

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets open cross section inflow at a SWMM node ID through a DLL call
    !> @param[in] self, id, inflow
    !---------------------------------------------------------------------------
    subroutine SetXSectionInflowByID(self, id, inflow)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    real(c_double), intent(in) :: inflow

#ifdef _SWMMCoupler_
    call swmm_setOpenXSectionInflow(id, inflow)
#endif

    end subroutine SetXSectionInflowByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets open cross section inflow at all suited SWMM nodes
    !> @param[in] self, inflow
    !---------------------------------------------------------------------------
    subroutine SetXSectionInflow(self, inflow, cellIDs)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: inflow
    integer, dimension(:), intent(in) :: cellIDs
    integer :: i, cid, nnodes, cidx

    if (size(self%xsectionLevelsIDX)>0) then
        do i=1, size(self%xsectionLevelsIDX)
            cid = self%n2cMap(self%xsectionLevelsIDX(i), cellID) !cell id from this node
            nnodes = count(self%n2cMap(:,cellID) == cid)         !number of nodes sharing the same cell
            cidx = findloc(cellIDs, value = cid, dim = 1)        !index of the cell in the data array
            if (cidx > 0) call self%SetXSectionInflowByID(self%xsectionLevelsIDX(i), inflow(cidx)/nnodes)
        end do
    end if
    end subroutine SetXSectionInflow
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets mohid cell list for a number of SWMM nodes requiring cell data
    !> @param[in] self, dataName, cellIDs
    !---------------------------------------------------------------------------
    subroutine GetCellList(self, dataName, cellIDs)
    class(swmm_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: dataName
    integer, dimension(:), allocatable, intent(inout) :: cellIDs

    if (dataName == 'Outfalls') then
        allocate(cellIDs(count(self%n2cMap(:,isOutfall) == 1)))
        where(self%n2cMap(:,isOutfall) == 1) cellIDs = self%n2cMap(:,cellID)
    else
        stop 'swmm_coupler_class::GetCellList - requested data not mapped for coupling'
    end if

    end subroutine GetCellList


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Prints the SWMM Coupler object
    !---------------------------------------------------------------------------
    subroutine printSWMMCoupler(self)
    class(swmm_coupler_class), intent(in) :: self

    print*, 'SWMM Coupler'
    print*, 'Initialized - ', self%initialized

    end subroutine printSWMMCoupler


    end module SWMMCoupler