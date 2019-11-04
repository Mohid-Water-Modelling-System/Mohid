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
    procedure :: print => printExternalCoupler
    !control procedures
    procedure :: runStep
    !import data procedures
    procedure :: getCoupledDt
    !export data procedures
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
    !> Prints the External Coupler object
    !---------------------------------------------------------------------------
    subroutine printExternalCoupler(self)
    class(external_coupler_class), intent(in) :: self

    print*, 'External Coupler Module'
    print*, 'Initialized - ', self%initialized
    print*, 'Couple to SWMM - ', self%swmmCoupling

    end subroutine printExternalCoupler


    end module ModuleExternalCoupler