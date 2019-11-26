    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Land
    ! MODULE        : External Coupler
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : October 2019
    ! REVISION      : Ricardo Canelas
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION
    !> Module that provides a generic interface to specfic coupler implementations
    !> it should instantiate and control any coupling classes and provide a single
    !> API for the basin module
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

    module ModuleExternalCoupler

    use ModuleGlobalData
    use ModuleHorizontalGrid
    use ModuleSewerGEMSEngineCoupler

    implicit none
    private

    type :: external_coupler_class       !< External Coupler class
        logical :: initialized = .false.         !< initialized flag
        logical :: SewerGEMSEngineCoupling = .false.        !< use SewerGEMSEngine coupler flag
        type(SewerGEMSEngine_coupler_class) :: SewerGEMSEngineCoupler
    contains
    procedure :: initialize => initExternalCoupler
    procedure :: initializeCouplerToModel
    procedure :: finalize => finalizeExternalCoupler
    procedure :: mapElements
    procedure :: isModelCoupled
    procedure :: print => printExternalCoupler
    !control procedures
    procedure :: runStep
    !import data procedures
    procedure :: setValues
    !export data procedures
    procedure :: getValues
    procedure :: GetCellList
    procedure :: getCoupledDt
    end type external_coupler_class

    !Public access vars
    public :: external_coupler_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Initializes the External Coupler object
    !---------------------------------------------------------------------------
    subroutine initExternalCoupler(self)
    class(external_coupler_class), intent(inout) :: self

    self%initialized = .true.

    end subroutine initExternalCoupler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Initializes an External Coupler object
    !---------------------------------------------------------------------------
    subroutine initializeCouplerToModel(self, modelName, HorizontalGridID)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName
    real, allocatable, dimension(:,:) :: mapArrayXY
    integer, allocatable, dimension(:,:) :: mapArrayIJ
    integer, allocatable, dimension(:) :: mapArrayID
    integer, intent(inout) :: HorizontalGridID

    if (self%initialized) then
        if (modelName == 'SewerGEMSEngine') then
            self%SewerGEMSEngineCoupling = .true.
            !initialize the model coupler
            call self%SewerGEMSEngineCoupler%initialize(mapArrayXY, mapArrayIJ, mapArrayID)
            call GetXYArrayIJ(HorizontalGridID, mapArrayXY, mapArrayIJ)
            call GetCellIDfromIJArray(HorizontalGridID, mapArrayIJ, mapArrayID)
            call self%mapElements(modelName, mapArrayIJ, mapArrayID)
        end if
        !add more models here
    end if

    end subroutine initializeCouplerToModel
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Finalizes an External Coupler object and all contained model objects
    !---------------------------------------------------------------------------
    subroutine finalizeExternalCoupler(self)
    class(external_coupler_class), intent(inout) :: self
    
    if (self%initialized) then
        if (self%SewerGEMSEngineCoupling) then
            self%SewerGEMSEngineCoupling = .false.
            call self%SewerGEMSEngineCoupler%finalize()            
        end if
        !add more models here
        
        self%initialized = .false.
    end if
    
    
    end subroutine finalizeExternalCoupler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> maps elements in an External Coupler object
    !---------------------------------------------------------------------------
    subroutine mapElements(self, modelName, mapArrayIJ, mapArrayID)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName
    integer, dimension(:,:), intent(inout) :: mapArrayIJ
    integer, dimension(:), intent(inout) :: mapArrayID

    if (self%initialized) then
        if (modelName == 'SewerGEMSEngine') then
            !initialize the coupler mapper
            call self%SewerGEMSEngineCoupler%mapElements(mapArrayIJ, mapArrayID)
        end if
        !add more models here
    end if

    end subroutine mapElements

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Returns coupling status of a model
    !---------------------------------------------------------------------------
    logical function isModelCoupled(self, modelName)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName

    isModelCoupled = .false.
    if (self%initialized) then
        if (modelName == 'SewerGEMSEngine') then
            if (self%SewerGEMSEngineCoupler%initialized) isModelCoupled = .true.
            return
        end if
        !add more models here
    end if

    end function isModelCoupled

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Runs a time step on all coupled models
    !---------------------------------------------------------------------------
    subroutine runStep(self, dt)
    class(external_coupler_class), intent(inout) :: self
    real, intent(in) :: dt

    if (self%SewerGEMSEngineCoupler%initialized) then
        call self%SewerGEMSEngineCoupler%runStep(dt)
    end if
    !add more models here

    end subroutine runStep

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> returns the smallest dt of all coupled models
    !---------------------------------------------------------------------------
    real function getCoupledDt(self)
    class(external_coupler_class), intent(inout) :: self
    real :: localDt

    getCoupledDt = -null_real

    if (self%SewerGEMSEngineCoupler%initialized) then
        localDt = self%SewerGEMSEngineCoupler%GetDt()
        getCoupledDt = min(localDt, getCoupledDt)
    end if
    !add more models here

    end function getCoupledDt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> returns queries values from externaly coupled object model
    !> @param[in] self, modelName, HorizontalGridID, dataName, dataArray
    !---------------------------------------------------------------------------
    subroutine getValues(self, modelName, HorizontalGridID, dataName, dataArray)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName
    integer, intent(in) :: HorizontalGridID
    character(len = StringLength), intent(in) :: dataName
    real, allocatable, dimension(:,:) :: dataArray

    if (modelName == 'SewerGEMSEngine') then
        if (dataName == 'Inflow') then
            allocate(dataArray(size(self%SewerGEMSEngineCoupler%inflowIDX),3))
            dataArray = self%SewerGEMSEngineCoupler%GetInflow(HorizontalGridID)
            return
        else if (dataName == 'Outflow') then
            allocate(dataArray(size(self%SewerGEMSEngineCoupler%outfallIDX),3))
            dataArray = self%SewerGEMSEngineCoupler%GetOutflow(HorizontalGridID)
            return
        else if (dataName == 'xLevel') then
            allocate(dataArray(size(self%SewerGEMSEngineCoupler%xsectionLevelsIDX),3))
            dataArray = self%SewerGEMSEngineCoupler%GetLevel(HorizontalGridID)
            return
        else
            stop 'external_coupler_class::getValues - requested quantity does not exist'
        end if
    end if
    !add more models here

    end subroutine getValues

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> exports variables to externaly coupled object model from MOHID
    !> @param[in] self, modelName, dataName, dataArray, cellIDs
    !---------------------------------------------------------------------------
    subroutine setValues(self, modelName, dataName, dataArray, cellIDs)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName
    character(len = StringLength), intent(in) :: dataName
    real, dimension(:) :: dataArray
    integer, dimension(:) :: cellIDs

    if (modelName == 'SewerGEMSEngine') then
        if (dataName == 'OutletLevel') then
            call self%SewerGEMSEngineCoupler%SetOutletLevel(dataArray, cellIDs)
            return
        else if (dataName == 'LateralInflow') then
            call self%SewerGEMSEngineCoupler%SetLateralInflow(dataArray, cellIDs)
            return
        else if (dataName == 'WaterColumn') then
            call self%SewerGEMSEngineCoupler%SetWaterColumn(dataArray, cellIDs)
            return
        else if (dataName == 'InletInflow') then
            call self%SewerGEMSEngineCoupler%SetInletInflow(dataArray, cellIDs)
            return
        else if (dataName == 'XSectionFlow') then
            call self%SewerGEMSEngineCoupler%SetXSectionInflow(dataArray, cellIDs)
            return
        else
            stop 'external_coupler_class::setValues - requested quantity does not exist'
        end if
    end if
    !add more models here

    end subroutine setValues
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets mohid cell list for a specific variable transfer
    !> @param[in] self, modelName, dataName, cellIDs
    !---------------------------------------------------------------------------
    subroutine GetCellList(self, modelName, dataName, cellIDs)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName
    character(len = StringLength), intent(in) :: dataName
    integer, dimension(:), allocatable, intent(inout) :: cellIDs

    if (modelName == 'SewerGEMSEngine') then
        call self%SewerGEMSEngineCoupler%GetCellList(dataName, cellIDs)
        return        
    end if
    !add more models here

    end subroutine GetCellList

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Prints the External Coupler object
    !---------------------------------------------------------------------------
    subroutine printExternalCoupler(self)
    class(external_coupler_class), intent(in) :: self

    print*, 'External Coupler Module'
    print*, 'Initialized - ', self%initialized
    print*, 'Couple to SewerGEMS Engine - ', self%SewerGEMSEngineCoupling

    end subroutine printExternalCoupler


    end module ModuleExternalCoupler