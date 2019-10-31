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

    subroutine swmm_getNumberOfNodes(nNodes) bind(C, name='swmm_getNumberOfNodes')
    use iso_c_binding
    integer(c_int) :: nNodes
    end subroutine swmm_getNumberOfNodes

    subroutine swmm_getNodeType(id, nType) bind(C, name='swmm_getNodeTypeByID')
    use iso_c_binding
    integer(c_int) :: id
    integer(c_int) :: nType
    end subroutine swmm_getNodeType

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
    

    type :: NodeTypes_enum          !< enums for node types
        integer :: JUNCTION = 0
        integer :: OUTFALL = 1
    end type NodeTypes_enum


    !main public class
    type :: swmm_coupler_class                      !< SWMM Coupler class
        logical :: initialized = .false.                        !< initialized flag
        type(NodeTypes_enum) :: NodeTypes                       !< node type flags
        integer :: NumberOfNodes                                !< number of SWMM nodes
        integer, allocatable, dimension(:) :: junctionIDX       !< ids of junction SWMM nodes
        integer, allocatable, dimension(:) :: outfallIDX        !< ids of outfall SWMM nodes
        integer, allocatable, dimension(:) :: inflowIDX         !< ids of inflow SWMM nodes
        integer, allocatable, dimension(:) :: xsectionLevelsIDX !< ids of open cross section SWMM nodes
        integer, allocatable, dimension(:) :: lateralFlowIDX    !< ids of lateral flow accepting SWMM nodes
        logical, allocatable, dimension(:) :: xSectionOpen      !warning: 1-based index array (c equivalent is 0-based)
        character(len = StringLength)      :: SWMM_dat
        character(len = StringLength)      :: SWMM_rpt
        character(len = StringLength)      :: SWMM_out
    contains
    procedure :: initialize => initSWMMCoupler
    procedure :: initializeSWMM
    procedure, private :: getFilesPaths
    procedure :: print => printSWMMCoupler
    !import data procedures
    procedure :: GetInflow
    procedure :: GetOutflow
    procedure :: GetLevel
    procedure, private :: GetNumberOfNodes
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
    !control procedures

    end type swmm_coupler_class


    !Public access vars
    public :: swmm_coupler_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Initializes the SWMM Coupler object
    !---------------------------------------------------------------------------
    subroutine initSWMMCoupler(self)
    class(swmm_coupler_class), intent(inout) :: self
    integer :: nJunction = 0
    integer :: nOutfall = 0
    integer :: nInflow = 0
    integer :: nXSection = 0
    integer :: nLatFlow = 0
    integer :: i, idx
    integer :: idxj, idxo, idxi, idxx, idxl = 1
    
    !real, allocatable, dimension(:) :: temp

    self%initialized = .true.
    !initialize SWMM
    call self%getFilesPaths()
    call self%initializeSWMM()

    call self%GetNumberOfNodes()

    print*, 'SWMM number of nodes is ', self%NumberOfNodes

    !building open section list
    allocate(self%xSectionOpen(self%NumberOfNodes))
    self%xSectionOpen = .false.
    do i=1, self%NumberOfNodes
        if (self%NodeTypes%junction == self%GetNodeTypeByID(i-1)) then
            if (self%GetIsNodeOpenChannel(i-1)) self%xSectionOpen(i) = .true.
        end if
    end do

    !building id lists for O(1) access
    do i=1, self%NumberOfNodes
        idx = i-1
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) nJunction = nJunction + 1
        if (self%NodeTypes%outfall  == self%GetNodeTypeByID(idx)) nOutfall  = nOutfall  + 1
        if (self%NodeTypes%junction  == self%GetNodeTypeByID(idx)) then
            if (.not.self%xSectionOpen(i)) nInflow = nInflow + 1     !only closed nodes
        end if
        if (self%NodeTypes%junction  == self%GetNodeTypeByID(idx)) then
            if (self%xSectionOpen(i)) nXSection = nXSection + 1      !only open nodes
        end if
        if (self%NodeTypes%junction  == self%GetNodeTypeByID(idx)) then
            if (self%GetNodeHasLateralInflow(idx)) nLatFlow = nLatFlow + 1
        end if
    end do

    !print*, "nJunction", nJunction
    !print*, "nOutfall", nOutfall
    !print*, "nInflow", nInflow
    !print*, "nXSection", nXSection
    !print*, "nLatFlow", nLatFlow

    allocate(self%junctionIDX(nJunction))
    allocate(self%outfallIDX(nOutfall))
    allocate(self%inflowIDX(nInflow))
    allocate(self%xsectionLevelsIDX(nXSection))
    allocate(self%lateralFlowIDX(nLatFlow))
    do i=1, self%NumberOfNodes
        idx = i-1
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            self%junctionIDX(idxj) = idx
            idxj = idxj + 1
        end if
        if (self%NodeTypes%outfall == self%GetNodeTypeByID(idx)) then
            self%outfallIDX(idxo) = idx
            idxo = idxo + 1
        end if
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            if (.not.self%xSectionOpen(i)) then !only closed nodes
                self%inflowIDX(idxi) = idx
                idxi = idxi + 1
            end if
        end if
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            if (self%xSectionOpen(i)) then      !only open nodes
                self%xsectionLevelsIDX(idxx) = idx
                idxx = idxx + 1
            end if
        end if
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) then
            if (self%GetNodeHasLateralInflow(i)) then
                self%lateralFlowIDX(idxl) = idx
                idxl = idxl + 1
            end if
        end if
    end do

    !print*, "self%xSectionOpen", self%xSectionOpen
    !print*, "self%junctionIDX", self%junctionIDX
    !print*, "self%outfallIDX", self%outfallIDX
    !print*, "self%inflowIDX", self%inflowIDX
    !print*, "self%xsectionLevelsIDX", self%xsectionLevelsIDX    
    !print*, "self%lateralFlowIDX", self%lateralFlowIDX    
    !!print*, "self%GetLevel()", self%GetLevel()
    !
    !allocate(temp(size(self%xsectionLevelsIDX)))
    !temp = 10.0
    !call self%SetXSectionInflow(temp)
    !
    !print*, 'initSWMMCoupler done'
    !read(*,*)

    end subroutine initSWMMCoupler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the files required for SWMM run
    !---------------------------------------------------------------------------
    subroutine getFilesPaths(self)
    class(swmm_coupler_class), intent(inout) :: self
    integer :: STAT_CALL

    call ReadFileName('SWMM_DAT', self%SWMM_dat,                         &
        Message = "SWMM input file", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getFilesPaths - SWMM_DAT keyword not found on main file list'

    call ReadFileName('SWMM_RPT', self%SWMM_rpt,                         &
        Message = "SWMM report file", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getFilesPaths - SWMM_RPT keyword not found on main file list'

    call ReadFileName('SWMM_OUT', self%SWMM_out,                         &
        Message = "SWMM output file", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'SWMMCoupler::getFilesPaths - SWMM_OUT keyword not found on main file list'

    end subroutine getFilesPaths

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Initializes the SWMM model through a DLL call
    !---------------------------------------------------------------------------
    subroutine initializeSWMM(self)
    class(swmm_coupler_class), intent(inout) :: self
    character(len = :, kind = c_char), allocatable :: inFile, rptFile, outFile

    print*, 'Initializing SWMM, please wait...'

    inFile = trim(ADJUSTL(self%SWMM_dat))//C_NULL_CHAR
    rptFile = trim(ADJUSTL(self%SWMM_rpt))//C_NULL_CHAR
    outFile = trim(ADJUSTL(self%SWMM_out))//C_NULL_CHAR

    call swmm_open(inFile, rptFile, outFile)

    print*,''

    end subroutine initializeSWMM

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the number of SWMM nodes through a DLL call
    !---------------------------------------------------------------------------
    subroutine GetNumberOfNodes(self)
    class(swmm_coupler_class), intent(inout) :: self
    integer(c_int) :: nNodes          !fortran target

    call swmm_getNumberOfNodes(nNodes)
    self%NumberOfNodes = nNodes

    end subroutine GetNumberOfNodes

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

    call swmm_getNodeType(id, nType)
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
    call swmm_getIsNodeOpenChannel(id, isOpen)
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
    call swmm_getNodeHasLateralInflow(id, isOpen)
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

    call swmm_getInflowByNode(id, inflow)
    GetInflowByID = inflow

    end function GetInflowByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets ponded inflow at all required SWMM nodes
    !---------------------------------------------------------------------------
    function GetInflow(self) result(inflow)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(size(self%inflowIDX)) :: inflow
    integer :: i

    inflow = 0.0
    if (size(self%inflowIDX)>0) then
        do i=1, size(self%inflowIDX)
            inflow(i) = self%GetInflowByID(self%inflowIDX(i))
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

    call swmm_getOutflowByNode(id, outflow)
    GetOutflowByID = outflow

    end function GetOutflowByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets outflow at all required SWMM nodes
    !---------------------------------------------------------------------------
    function GetOutflow(self) result(outflow)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(size(self%outfallIDX)) :: outflow
    integer :: i

    outflow = 0.0
    if (size(self%outfallIDX)>0) then
        do i=1, size(self%outfallIDX)
            outflow(i) = self%GetOutflowByID(self%outfallIDX(i))
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

    call swmm_getLevelByNode(id, level)
    GetLevelByID = level

    end function GetLevelByID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets level at all required SWMM nodes
    !---------------------------------------------------------------------------
    function GetLevel(self) result(level)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(size(self%xsectionLevelsIDX)) :: level
    integer :: i

    level = 0.0
    if (size(self%xsectionLevelsIDX)>0) then
        do i=1, size(self%xsectionLevelsIDX)
            level(i) = self%GetLevelByID(self%xsectionLevelsIDX(i))
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

    call swmm_setDownstreamWaterLevel(id, outletLevel)

    end subroutine SetOutletLevelByID
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets outlet level at all suited SWMM nodes
    !> @param[in] self, outletLevel
    !---------------------------------------------------------------------------
    subroutine SetOutletLevel(self, outletLevel)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: outletLevel
    integer :: i

    if (size(self%outfallIDX)>0) then
        if (size(self%outfallIDX) == size(outletLevel)) then
            do i=1, size(self%outfallIDX)
                call self%SetOutletLevelByID(self%outfallIDX(i), outletLevel(i))
            end do
        else
            stop 'SWMMCoupler::SetOutletLevel - size mismatch between data array and id array'
        end if        
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

    call swmm_setLateralInflow(id, inflow)

    end subroutine SetLateralInflowByID
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets outlet level at all suited SWMM nodes
    !> @param[in] self, inflow
    !---------------------------------------------------------------------------
    subroutine SetLateralInflow(self, inflow)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: inflow
    integer :: i

    if (size(self%junctionIDX)>0) then
        if (size(self%junctionIDX) == size(inflow)) then
            do i=1, size(self%junctionIDX)
                call self%SetLateralInflowByID(self%junctionIDX(i), inflow(i))
            end do
        else
            stop 'SWMMCoupler::SetLateralInflow - size mismatch between data array and id array'
        end if        
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

    call swmm_setPondedWaterColumn(id, level)

    end subroutine SetWaterColumnByID
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets outlet level at all suited SWMM nodes - not checking if node 
    !> should have or not
    !> @param[in] self, level
    !---------------------------------------------------------------------------
    subroutine SetWaterColumn(self, level)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: level
    integer :: i

    if (size(self%inflowIDX)>0) then
        if (size(self%inflowIDX) == size(level)) then
            do i=1, size(self%inflowIDX)
                call self%SetWaterColumnByID(self%inflowIDX(i), level(i))
            end do
        else
            stop 'SWMMCoupler::SetWaterColumn - size mismatch between data array and id array'
        end if        
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

    call swmm_setStormWaterPotentialInflow(id, inflow)

    end subroutine SetInletInflowByID
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets inflow at all suited SWMM nodes - not checking if node 
    !> should have or not
    !> @param[in] self, inflow
    !---------------------------------------------------------------------------
    subroutine SetInletInflow(self, inflow)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: inflow
    integer :: i

    if (size(self%inflowIDX)>0) then
        if (size(self%inflowIDX) == size(inflow)) then
            do i=1, size(self%inflowIDX)
                call self%SetInletInflowByID(self%inflowIDX(i), inflow(i))
            end do
        else
            stop 'SWMMCoupler::SetInletInflow - size mismatch between data array and id array'
        end if        
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

    call swmm_setOpenXSectionInflow(id, inflow)

    end subroutine SetXSectionInflowByID
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Sets open cross section inflow at all suited SWMM nodes
    !> @param[in] self, inflow
    !---------------------------------------------------------------------------
    subroutine SetXSectionInflow(self, inflow)
    class(swmm_coupler_class), intent(in) :: self
    real, dimension(:), intent(in) :: inflow
    integer :: i

    if (size(self%xsectionLevelsIDX)>0) then
        if (size(self%xsectionLevelsIDX) == size(inflow)) then
            do i=1, size(self%xsectionLevelsIDX)
                call self%SetXSectionInflowByID(self%xsectionLevelsIDX(i), inflow(i))
            end do
        else
            stop 'SWMMCoupler::SetXSectionInflow - size mismatch between data array and id array'
        end if        
    end if

    end subroutine SetXSectionInflow


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