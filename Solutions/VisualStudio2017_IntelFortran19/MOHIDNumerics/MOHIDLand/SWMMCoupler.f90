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
    use iso_c_binding

    implicit none
    private


    interface dllFuncs

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

    end interface dllFuncs


    type :: NodeTypes_enum          !< enums for node types
        integer :: JUNCTION = 0
        integer :: OUTFALL = 1
    end type NodeTypes_enum


    !main public class
    type :: swmm_coupler_class                      !< SWMM Coupler class
        logical :: initialized = .false.                    !< initialized flag
        type(NodeTypes_enum) :: NodeTypes                   !< node type flags
        integer :: NumberOfNodes                            !< number of SWMM nodes
        integer, allocatable, dimension(:) :: junctionIDX   !< ids of junction SWMM nodes
        integer, allocatable, dimension(:) :: outfallIDX    !< ids of outfall SWMM nodes
        integer, allocatable, dimension(:) :: inflowIDX     !< ids of inflow SWMM nodes
        logical, allocatable, dimension(:) :: xSectionOpen  !warning: 1-based index array (c equivalent is 0-based)
    contains
    procedure :: initialize => initSWMMCoupler
    procedure :: print => printSWMMCoupler
    !import data procedures
    procedure, private :: GetNumberOfNodes
    procedure, private :: GetNodeTypeByID
    procedure, private :: GetIsNodeOpenChannel
    procedure :: GetInflow
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
    integer :: i, idx
    integer :: idxj, idxo, idxi = 1

    self%initialized = .true.
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

    !building id lists
    do i=1, self%NumberOfNodes
        idx = i-1
        if (self%NodeTypes%junction == self%GetNodeTypeByID(idx)) nJunction = nJunction + 1
        if (self%NodeTypes%outfall  == self%GetNodeTypeByID(idx)) nOutfall  = nOutfall  + 1
        if (self%NodeTypes%junction  == self%GetNodeTypeByID(idx)) then
            if (self%xSectionOpen(i)) nInflow = nInflow + 1
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
            if (self%xSectionOpen(i)) then
                self%inflowIDX(idxi) = idx
                idxi = idxi + 1
            end if            
        end if
    end do

    end subroutine initSWMMCoupler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the number of SWMM nodes by a DLL call
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
    !> Gets the SWMM node type by a DLL call
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
    !> Gets the open channel status of a SWMM node by a DLL call
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
    !> Gets the open channel status of a SWMM node by a DLL call
    !---------------------------------------------------------------------------
    function GetInflow(self) result(inflow)
    class(swmm_coupler_class), intent(in) :: self
    real, allocatable, dimension(:) :: inflow
    type(c_ptr), allocatable, dimension(:) :: inflowC   !c pointer
    real, allocatable, dimension(:), target :: inflowF  !fortran target
    integer :: i
    
    !allocate()
    !do i=1, size(self%inflowIDX)
        
    
    end function GetInflow


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