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
    use SWMMCoupler

    implicit none
    private

    type :: external_coupler_class       !< External Coupler class
        logical :: initialized = .false.         !< initialized flag
        logical :: swmmCoupling = .false.        !< use swmm coupler flag
        type(swmm_coupler_class) :: SWMMCoupler
    contains
    procedure :: initialize => initExternalCoupler
    procedure :: initializeCouplerToModel
    procedure :: isModelCoupled
    procedure :: print => printExternalCoupler
    !control procedures
    procedure :: runStep
    !import data procedures
    procedure :: setValues
    !export data procedures
    procedure :: getValues    
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
    !> Initializes the External Coupler object
    !---------------------------------------------------------------------------
    subroutine initializeCouplerToModel(self, modelName)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName

    if (self%initialized) then
        if (modelName == 'SWMM') then
            self%swmmCoupling = .true.
            !initialize the model coupler
            call self%SWMMCoupler%initialize()
        end if
        !add more models here
    end if

    end subroutine initializeCouplerToModel
    
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
        if (modelName == 'SWMM') then            
            if (self%SWMMCoupler%initialized) isModelCoupled = .true.
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

     if (self%SWMMCoupler%initialized) then
        call self%SWMMCoupler%runStep(dt)
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
    
    if (self%SWMMCoupler%initialized) then
        localDt = self%SWMMCoupler%GetDt()
        getCoupledDt = min(localDt, getCoupledDt)        
    end if
    !add more models here

    end function getCoupledDt
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> returns queries values from externaly coupled object
    !> @param[in] self, modelName, dataName, dataArray
    !---------------------------------------------------------------------------
    subroutine getValues(self, modelName, dataName, dataArray)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName
    character(len = StringLength), intent(in) :: dataName
    real, allocatable, dimension(:) :: dataArray

    if (modelName == 'SWMM') then
        if (dataName == 'Inflow') then
            allocate(dataArray(size(self%SWMMCoupler%inflowIDX)))
            dataArray = self%SWMMCoupler%GetInflow()
            return
        else if (dataName == 'Outflow') then
            allocate(dataArray(size(self%SWMMCoupler%outfallIDX)))
            dataArray = self%SWMMCoupler%GetOutflow()
            return
        else if (dataName == 'xLevel') then
            allocate(dataArray(size(self%SWMMCoupler%xsectionLevelsIDX)))
            dataArray = self%SWMMCoupler%GetLevel()
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
    !> exports variables to externaly coupled object
    !> @param[in] self, modelName, dataName, dataArray
    !---------------------------------------------------------------------------
    subroutine setValues(self, modelName, dataName, dataArray)
    class(external_coupler_class), intent(inout) :: self
    character(len = StringLength), intent(in) :: modelName
    character(len = StringLength), intent(in) :: dataName
    real, dimension(:) :: dataArray

    if (modelName == 'SWMM') then
        if (dataName == 'OutletLevel') then
            call self%SWMMCoupler%SetOutletLevel(dataArray)
            return
        else if (dataName == 'LateralInflow') then
            call self%SWMMCoupler%SetLateralInflow(dataArray)
            return
        else if (dataName == 'WaterColumn') then
            call self%SWMMCoupler%SetWaterColumn(dataArray)
            return
        else if (dataName == 'InletInflow') then
            call self%SWMMCoupler%SetInletInflow(dataArray)
            return
        else if (dataName == 'XSectionFlow') then
            call self%SWMMCoupler%SetXSectionInflow(dataArray)
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
    !> Prints the External Coupler object
    !---------------------------------------------------------------------------
    subroutine printExternalCoupler(self)
    class(external_coupler_class), intent(in) :: self

    print*, 'External Coupler Module'
    print*, 'Initialized - ', self%initialized
    print*, 'Couple to SWMM - ', self%swmmCoupling

    end subroutine printExternalCoupler


    end module ModuleExternalCoupler