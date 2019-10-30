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


    interface dllFuncs
    
    subroutine swmm_open(inFile, rptFile, outFile) bind(C, name='swmm_open')
    use iso_c_binding
    type(c_ptr) :: inFile
    type(c_ptr) :: rptFile
    type(c_ptr) :: outFile
    end subroutine swmm_open

    subroutine swmm_getNumberOfNodes(nNodes) bind(C, name='swmm_getNumberOfNodes')
    use iso_c_binding
    type(c_ptr) :: nNodes
    end subroutine swmm_getNumberOfNodes

    subroutine swmm_getNodeType(id) bind(C, name='swmm_getNodeType')
    use iso_c_binding
    type(c_ptr) :: id
    end subroutine swmm_getNodeType

    subroutine swmm_getIsNodeOpenChannel(id1) bind(C, name='swmm_getIsNodeOpenChannel')
    use iso_c_binding
    type(c_ptr) :: id1
    end subroutine swmm_getIsNodeOpenChannel

    subroutine swmm_getInflowByNode(id2) bind(C, name='swmm_getInflowByNode')
    use iso_c_binding
    type(c_ptr) :: id2
    end subroutine swmm_getInflowByNode

    subroutine swmm_getOutflowByNode(id3) bind(C, name='swmm_getOutflowByNode')
    use iso_c_binding
    type(c_ptr) :: id3
    end subroutine swmm_getOutflowByNode

    subroutine swmm_getLevelByNode(id4) bind(C, name='swmm_getLevelByNode')
    use iso_c_binding
    type(c_ptr) :: id4
    end subroutine swmm_getLevelByNode

    end interface dllFuncs


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
        integer, allocatable, dimension(:) :: xsectionLevelsIDX !< ids of open cros section SWMM nodes
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
    procedure, private :: GetInflowByID
    procedure, private :: GetOutflowByID
    procedure, private :: GetLevelByID
    !export data procedures
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
    integer :: i, idx
    integer :: idxj, idxo, idxi, idxx = 1

    self%initialized = .true.
    !initialize SWMM
    call self%getFilesPaths()
    call self%initializeSWMM()
    
    call self%GetNumberOfNodes()

    !building open section list
    allocate(self%xSectionOpen(self%NumberOfNodes))
    self%xSectionOpen = .false.
    do i=1, self%NumberOfNodes
        idx = i-1
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
    end do
    allocate(self%junctionIDX(nJunction))
    allocate(self%outfallIDX(nOutfall))
    allocate(self%inflowIDX(nInflow))
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
    end do

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
    type(c_ptr) :: inFileC, rptFileC, outFileC      !c pointer
    character(len = StringLength), target :: inFileF, rptFileF, outFileF    !fortran target
    
    print*, 'Initializing SWMM, please wait...'
    
    inFileF = trim(ADJUSTL(self%SWMM_dat))//C_NULL_CHAR
    rptFileF = trim(ADJUSTL(self%SWMM_rpt))//C_NULL_CHAR
    outFileF = trim(ADJUSTL(self%SWMM_out))//C_NULL_CHAR
    
    inFileF = "hello world"//C_NULL_CHAR
    
    print*, inFileF
    
    inFileC = c_loc(inFileF)            !pointing c pointer to fortran target location
    rptFileC = c_loc(rptFileF)
    outFileC = c_loc(outFileF)
    
    call swmm_open(inFileC, rptFileC, outFileC)
    
    print*,''

    end subroutine initializeSWMM
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the number of SWMM nodes through a DLL call
    !---------------------------------------------------------------------------
    subroutine GetNumberOfNodes(self)
    class(swmm_coupler_class), intent(inout) :: self
    type(c_ptr) :: nNodesC              !c pointer
    integer, target :: nNodesF          !fortran target

    nNodesC = c_loc(nNodesF)            !pointing c pointer to fortran target location
    !call swmm_getNumberOfNodes(nNodesC)
    self%NumberOfNodes = nNodesF

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
    type(c_ptr) :: nTypeC              !c pointer
    integer, target :: nTypeF          !fortran target

    nTypeC = c_loc(nTypeF)             !pointing c pointer to fortran target location
    !call swmm_getNodeType(id, nTypeC)
    GetNodeTypeByID = nTypeF

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
    type(c_ptr) :: isOpenC              !c pointer
    integer, target :: isOpenF          !fortran target

    GetIsNodeOpenChannel = .false.
    isOpenC = c_loc(isOpenF)             !pointing c pointer to fortran target location
    !call swmm_getIsNodeOpenChannel(id, isOpenC)
    if (isOpenF == 1) GetIsNodeOpenChannel = .true.

    end function GetIsNodeOpenChannel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets ponded inflow by SWMM node ID through a DLL call
    !> @param[in] self, id
    !---------------------------------------------------------------------------
    real function GetInflowByID(self, id)
    class(swmm_coupler_class), intent(in) :: self
    integer(c_int), intent(in) :: id
    type(c_ptr) :: inflowC   !c pointer
    real, target :: inflowF  !fortran target

    inflowC = c_loc(inflowF)             !pointing c pointer to fortran target location
    !call swmm_getInflowByNode(id, inflowC)
    GetInflowByID = inflowF

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
    type(c_ptr) :: outflowC   !c pointer
    real, target :: outflowF  !fortran target

    outflowC = c_loc(outflowF)             !pointing c pointer to fortran target location
    !call swmm_getOutflowByNode(id, outflowC)
    GetOutflowByID = outflowF

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
    type(c_ptr) :: levelC   !c pointer
    real, target :: levelF  !fortran target

    levelC = c_loc(levelF)             !pointing c pointer to fortran target location
    !call swmm_getLevelByNode(id, levelC)
    GetLevelByID = levelF

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
    !> Prints the SWMM Coupler object
    !---------------------------------------------------------------------------
    subroutine printSWMMCoupler(self)
    class(swmm_coupler_class), intent(in) :: self

    print*, 'SWMM Coupler'
    print*, 'Initialized - ', self%initialized

    end subroutine printSWMMCoupler


    end module SWMMCoupler