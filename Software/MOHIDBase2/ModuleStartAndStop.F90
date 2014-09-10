!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : Snow
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jul 2014
! REVISION      : Eduardo Jauch - v4.0
! DESCRIPTION   : Module to deal with the continuation of runs
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
   
module ModuleStartAndStop

    use ModuleGlobalData
    use ModuleCollection
    use ModuleTime
    use ModuleHDF5
    use ModuleFunctions,        only :  TimeToString, ChangeSuffix   
    use ModuleHorizontalGrid,   only :  WriteHorizontalGrid
   
    implicit none
    private 
   
    !Public Elements-----------------------------------------------------------        
    public  ::  T_StartAndStop    
    
    !Types---------------------------------------------------------------------
    type T_Property        
        integer                                         ::  Id          =   -1
        character(len=StringLength)                     ::  Name        =   '*****'
        character(len=StringLength)                     ::  Group       =   '*****'
        character(len=StringLength)                     ::  Units       =   '*****'
        logical                                         ::  Load        =   .true.
        logical                                         ::  Save        =   .false.
        integer                                         ::  ShapeType   =   1   !1 -> 1D, 2-> 1D !3 -> 3D
    end type T_Property

    type, extends(T_Property) :: T_PropertyInteger0D
        integer, pointer                                ::  Values
    end type T_PropertyInteger0D

    type, extends(T_Property) :: T_PropertyInteger1D
        integer, dimension(:), pointer                  ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1
    end type T_PropertyInteger1D

    type, extends(T_Property) :: T_PropertyInteger2D
        integer, dimension(:,:), pointer                ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1,         &
                                                            JLB         =   -1,         &
                                                            JUB         =   -1
    end type T_PropertyInteger2D

    type, extends(T_Property) :: T_PropertyInteger3D
        integer, dimension(:,:,:), pointer              ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1,         &
                                                            JLB         =   -1,         &
                                                            JUB         =   -1,         &
                                                            KLB         =   -1,         &
                                                            KUB         =   -1
    end type T_PropertyInteger3D
    
    type, extends(T_Property) :: T_PropertyReal0D
        real(4), pointer                                ::  Values
    end type T_PropertyReal0D
    
    type, extends(T_Property) :: T_PropertyReal1D
        real(4), dimension(:), pointer                  ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1        
    end type T_PropertyReal1D

    type, extends(T_Property) :: T_PropertyReal2D
        real(4), dimension(:,:), pointer                ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1,         &
                                                            JLB         =   -1,         &
                                                            JUB         =   -1
    end type T_PropertyReal2D

    type, extends(T_Property) :: T_PropertyReal3D
        real(4), dimension(:,:,:), pointer              ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1,         &
                                                            JLB         =   -1,         &
                                                            JUB         =   -1,         &
                                                            KLB         =   -1,         &
                                                            KUB         =   -1
    end type T_PropertyReal3D    
    
    type, extends(T_Property) :: T_PropertyDouble0D
        real(8), pointer                                ::  Values
    end type T_PropertyDouble0D
    
    type, extends(T_Property) :: T_PropertyDouble1D
        real(8), dimension(:), pointer                  ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1
    end type T_PropertyDouble1D

    type, extends(T_Property) :: T_PropertyDouble2D
        real(8), dimension(:,:), pointer                ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1,         &
                                                            JLB         =   -1,         &
                                                            JUB         =   -1
    end type T_PropertyDouble2D

    type, extends(T_Property) :: T_PropertyDouble3D
        real(8), dimension(:,:,:), pointer              ::  Values
        integer                                         ::  ILB         =   -1,         &
                                                            IUB         =   -1,         &
                                                            JLB         =   -1,         &
                                                            JUB         =   -1,         &
                                                            KLB         =   -1,         &
                                                            KUB         =   -1
    end type T_PropertyDouble3D     
    
    type T_StartAndStop
        private
            type (C_Collection)                         ::  pProperties
            character (StringLength)                    ::  pMessage = ""
            
            logical, public                             ::  StopOnWrongDate = .true.            
            logical, public                             ::  WriteHorizGrid  = .false.            
            integer, public                             ::  ObjHorizGrid = 0            
            character (PathLength), public              ::  File = "*****"            
            
    contains
            procedure, private                          ::  AddI
            procedure, private                          ::  AddR4
            procedure, private                          ::  AddR8
            procedure, private                          ::  Add1DI
            procedure, private                          ::  Add1DR4
            procedure, private                          ::  Add1DR8
            procedure, private                          ::  Add2DI
            procedure, private                          ::  Add2DR4
            procedure, private                          ::  Add2DR8
            procedure, private                          ::  Add3DI
            procedure, private                          ::  Add3DR4
            procedure, private                          ::  Add3DR8            
            generic                                     ::  Add =>              &
                                                            AddI,               &
                                                            AddR4,              &
                                                            AddR8,              &
                                                            Add1DI,             &
                                                            Add1DR4,            &
                                                            Add1DR8,            &
                                                            Add2DI,             &
                                                            Add2DR4,            &
                                                            Add2DR8,            &
                                                            Add3DI,             &
                                                            Add3DR4,            &
                                                            Add3DR8
            
            procedure                                   ::  Load
            procedure                                   ::  Save

            procedure                                   ::  Remove            
            procedure                                   ::  Clear
    
            procedure                                   ::  Message
            
            procedure, private                          ::  InitializeObject
            final                                       ::  KillObject
    
    end type T_StartAndStop
        
    !Interfaces----------------------------------------------------------------
    
    interface T_StartAndStop
        module procedure Constructor
    end interface T_StartAndStop
    
    !Subroutines---------------------------------------------------------------
    contains
    
    function Constructor (file, stop_on_wrong_date, write_horiz_grid, obj_horiz_grid)
    
        !Arguments-------------------------------------------------------------
        type(T_StartAndStop)                            ::  Constructor
        character (*)                                   ::  file
        logical, optional, intent (in)                  ::  stop_on_wrong_date
        logical, optional, intent (in)                  ::  write_horiz_grid
        integer, optional, intent (in)                  ::  obj_horiz_grid
        
        !Local-----------------------------------------------------------------
        logical                                         ::  stop_on_wrong_date_
        logical                                         ::  write_horiz_grid_
        integer                                         ::  obj_horiz_grid_
        
        !----------------------------------------------------------------------
        if (present (stop_on_wrong_date)) then
            stop_on_wrong_date_ = stop_on_wrong_date
        else
            stop_on_wrong_date_ = .true.
        endif
        
        if (present (write_horiz_grid)) then
            write_horiz_grid_ = write_horiz_grid
            if (write_horiz_grid) then
                if (present (obj_horiz_grid)) then
                    obj_horiz_grid_ = obj_horiz_grid
                else
                    stop 'T_StartAndStop - Constructor - ERR010'
                endif
            endif
        else
            write_horiz_grid_ = .false.
            obj_horiz_grid_   = 0
        endif
        
        call constructor%InitializeObject (stop_on_wrong_date, write_horiz_grid_, obj_horiz_grid_, file)
        
    end function Constructor
    
    !--------------------------------------------------------------------------
    
   subroutine AddI (this, id, values, load, save, group, units, name)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        integer, pointer                                ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyInteger0D::prop)
        
        select type (prop)
        type is (T_PropertyInteger0D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  0
            prop%Values     => values
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine AddI
    
    !--------------------------------------------------------------------------
    
    subroutine AddR4 (this, id, values, load, save, group, units, name)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        real(4), pointer                                ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyReal0D::prop)
        
        select type (prop)
        type is (T_PropertyReal0D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  0
            prop%Values     => values
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine AddR4
    
    !--------------------------------------------------------------------------
    
    subroutine AddR8 (this, id, values, load, save, group, units, name)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        real(8), pointer                                ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyDouble0D::prop)
        
        select type (prop)
        type is (T_PropertyDouble0D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  0
            prop%Values     => values
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine AddR8
    
    !--------------------------------------------------------------------------

    subroutine Add1DI (this, id, values, load, save, group, units, name, ilb, iub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        integer, dimension(:), pointer                  ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyInteger1D::prop)
        
        select type (prop)
        type is (T_PropertyInteger1D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  1
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add1DI
    
    !--------------------------------------------------------------------------
    
    subroutine Add1DR4 (this, id, values, load, save, group, units, name, ilb, iub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        real(4), dimension(:), pointer                  ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyReal1D::prop)
        
        select type (prop)
        type is (T_PropertyReal1D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  1
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add1DR4
    
    !--------------------------------------------------------------------------
    
    subroutine Add1DR8 (this, id, values, load, save, group, units, name, ilb, iub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        real(8), dimension(:), pointer                  ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyDouble1D::prop)
        
        select type (prop)
        type is (T_PropertyDouble1D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  1
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add1DR8
    
    !--------------------------------------------------------------------------
    
    subroutine Add2DI (this, id, values, load, save, group, units, name, ilb, iub, jlb, jub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        integer, dimension(:,:), pointer                ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub,        &
                                                            jlb,        &
                                                            jub
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyInteger2D::prop)
        
        select type (prop)
        type is (T_PropertyInteger2D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  2
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            if (present(jlb)) then
                prop%JLB    =  jlb
            else
                prop%JLB    =  lbound(values, 2)
            endif
            if (present(jub)) then
                prop%JUB    = jub
            else
                prop%JUB    =  ubound(values, 2)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add2DI
    
    !--------------------------------------------------------------------------    

    subroutine Add2DR4 (this, id, values, load, save, group, units, name, ilb, iub, jlb, jub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        real(4), dimension(:,:), pointer                ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub,        &
                                                            jlb,        &
                                                            jub
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyReal2D::prop)
        
        select type (prop)
        type is (T_PropertyReal2D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  2
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            if (present(jlb)) then
                prop%JLB    =  jlb
            else
                prop%JLB    =  lbound(values, 2)
            endif
            if (present(jub)) then
                prop%JUB    = jub
            else
                prop%JUB    =  ubound(values, 2)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add2DR4
    
    !--------------------------------------------------------------------------

    subroutine Add2DR8 (this, id, values, load, save, group, units, name, ilb, iub, jlb, jub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        real(8), dimension(:,:), pointer                ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub,        &
                                                            jlb,        &
                                                            jub
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyDouble2D::prop)
        
        select type (prop)
        type is (T_PropertyDouble2D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  2
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            if (present(jlb)) then
                prop%JLB    =  jlb
            else
                prop%JLB    =  lbound(values, 2)
            endif
            if (present(jub)) then
                prop%JUB    = jub
            else
                prop%JUB    =  ubound(values, 2)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add2DR8
    
    !--------------------------------------------------------------------------

    subroutine Add3DI (this, id, values, load, save, group, units, name, ilb, iub, jlb, jub, klb, kub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        integer, dimension(:,:,:), pointer              ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub,        &
                                                            jlb,        &
                                                            jub,        &
                                                            klb,        &
                                                            kub         
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyInteger3D::prop)
        
        select type (prop)
        type is (T_PropertyInteger3D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  3
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            if (present(jlb)) then
                prop%JLB    =  jlb
            else
                prop%JLB    =  lbound(values, 2)
            endif
            if (present(jub)) then
                prop%JUB    = jub
            else
                prop%JUB    =  ubound(values, 2)
            endif
            if (present(klb)) then
                prop%KLB    =  klb
            else
                prop%KLB    =  lbound(values, 3)
            endif
            if (present(kub)) then
                prop%KUB    =  kub
            else
                prop%KUB    =  ubound(values, 3)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add3DI
    
    !--------------------------------------------------------------------------
    
    subroutine Add3DR4 (this, id, values, load, save, group, units, name, ilb, iub, jlb, jub, klb, kub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        real(4), dimension(:,:,:), pointer              ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub,        &
                                                            jlb,        &
                                                            jub,        &
                                                            klb,        &
                                                            kub
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyReal3D::prop)
        
        select type (prop)
        type is (T_PropertyReal3D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  3
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            if (present(jlb)) then
                prop%JLB    =  jlb
            else
                prop%JLB    =  lbound(values, 2)
            endif
            if (present(jub)) then
                prop%JUB    = jub
            else
                prop%JUB    =  ubound(values, 2)
            endif
            if (present(klb)) then
                prop%KLB    =  klb
            else
                prop%KLB    =  lbound(values, 3)
            endif
            if (present(kub)) then
                prop%KUB    =  kub
            else
                prop%KUB    =  ubound(values, 3)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add3DR4
    
    !--------------------------------------------------------------------------
    
    subroutine Add3DR8 (this, id, values, load, save, group, units, name, ilb, iub, jlb, jub, klb, kub)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        real(8), dimension(:,:,:), pointer              ::  values
        logical, optional                               ::  load,       &
                                                            save
        character(len=*), optional                      ::  group,      &
                                                            units,      &
                                                            name
        integer, optional                               ::  ilb,        &
                                                            iub,        &
                                                            jlb,        &
                                                            jub,        &
                                                            klb,        &
                                                            kub
        
        !Local-----------------------------------------------------------------
        class(*), pointer                               ::  prop
        
        !----------------------------------------------------------------------
        allocate (T_PropertyDouble3D::prop)
        
        select type (prop)
        type is (T_PropertyDouble3D)              
            prop%Id         =  id
            if (present(name)) then
                prop%Name   =  name
            else
                prop%Name   =  GetPropertyName (id)
            endif
            prop%ShapeType  =  3
            prop%Values     => values
            if (present(ilb)) then
                prop%ILB    =  ilb
            else
                prop%ILB    =  lbound(values, 1)
            endif
            if (present(iub)) then
                prop%IUB    =  iub
            else
                prop%IUB    =  ubound(values, 1)
            endif
            if (present(jlb)) then
                prop%JLB    =  jlb
            else
                prop%JLB    =  lbound(values, 2)
            endif
            if (present(jub)) then
                prop%JUB    = jub
            else
                prop%JUB    =  ubound(values, 2)
            endif
            if (present(klb)) then
                prop%KLB    =  klb
            else
                prop%KLB    =  lbound(values, 3)
            endif
            if (present(kub)) then
                prop%KUB    =  kub
            else
                prop%KUB    =  ubound(values, 3)
            endif
            
            if (present(load)) then
                prop%Load = load
            else
                prop%Load = .true.
            endif
            if (present(save)) then
                prop%Save = save
            else
                prop%Save = .true.
            endif
            if (present(group)) then
                prop%Group = group
            else
                prop%Group = "Results"
            endif
            if (present(units)) then
                prop%Units = units
            else
                prop%Units = "-"
            endif            
        end select
        
        call this%pProperties%Add (prop)
    
    end subroutine Add3DR8
    
    !--------------------------------------------------------------------------
        
    function Load (this, now)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        type(T_Time)                                    ::  now
        logical                                         ::  Load
        
        !Local-----------------------------------------------------------------
        logical                                         ::  exist
        integer                                         ::  obj_hdf5,   &
                                                            HDF5_READ,  &
                                                            stat
        class(*), pointer                               ::  item
        real, dimension(6), target                      ::  AuxTime
        real, dimension(:), pointer                     ::  PtrTime
        integer, dimension(1), target                   ::  AuxInteger
        integer, dimension(:), pointer                  ::  PtrInteger
        real(4), dimension(1), target                   ::  AuxReal
        real(4), dimension(:), pointer                  ::  PtrReal
        real(8), dimension(1), target                   ::  AuxDouble
        real(8), dimension(:), pointer                  ::  PtrDouble
        type(T_Time)                                    ::  TimeInFile
        
        !----------------------------------------------------------------------
        inquire (FILE=trim(this%File), EXIST=exist)
        
        if (exist .and. (.not. this%pProperties%IsEmpty())) then
                        
            PtrInteger  => AuxInteger
            PtrReal     => AuxReal
            PtrDouble   => AuxDouble
            
            obj_hdf5 = 0
            call GetHDF5FileAccess(HDF5_READ = HDF5_READ)            
            call ConstructHDF5(obj_hdf5, trim(this%File), HDF5_READ, STAT=stat)
            if (stat /= SUCCESS_) then
                Load = .false.
                return
            endif
            
            if (this%StopOnWrongDate) then
                call HDF5SetLimits (obj_hdf5, 1, 6, STAT=stat)
                if (stat /= SUCCESS_) then
                    Load = .false.
                    return
                endif
                            
                PtrTime => AuxTime
                call HDF5ReadData (obj_hdf5, "/Time",                                       &
                                   "Time", Array1D=PtrTime,                                 &
                                   STAT=stat)
                if (stat /= SUCCESS_) then
                    Load = .false.
                    return
                endif
                
                call SetDate (TimeInFile, AuxTime(1), AuxTime(2), AuxTime(3),               &
                                          AuxTime(4), AuxTime(5), AuxTime(6))
                if (abs(TimeInFile-now) > 0.01) then
                    Load = .false.
                    this%pMessage = "Time in initialization file is different from run start time"
                    return
                endif
            endif
            
            item => this%pProperties%First ()
do1:        do while (associated (item))
                select type (item)
                class is (T_Property)                    
                    if (item%Load) then
                        select type (item)
                        type is (T_PropertyInteger0D)
                            
                            call HDF5SetLimits (obj_hdf5, 1, 1, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                            
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array1D=PtrInteger,               &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif                            
                            item%Values = PtrInteger(1)
                            
                        type is (T_PropertyInteger1D)
                    
                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array1D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                        type is (T_PropertyInteger2D)
                    
                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB,                         &
                                                item%JLB, item%JUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                exit do1
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array2D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                        type is (T_PropertyInteger3D)
                    
                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB,                         &
                                                item%JLB, item%JUB,                         &
                                                item%KLB, item%KUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array3D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif

                        type is (T_PropertyReal0D)
                            
                            call HDF5SetLimits (obj_hdf5, 1, 1, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                            
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array1D=PtrReal,                  &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif                            
                            item%Values = PtrReal(1)
                            
                        type is (T_PropertyReal1D)
                    
                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array1D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                        type is (T_PropertyReal2D)
                    
                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB,                         &
                                                item%JLB, item%JUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array2D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                        type is (T_PropertyReal3D)
                    
                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB,                         &
                                                item%JLB, item%JUB,                         &
                                                item%KLB, item%KUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array3D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif

                        type is (T_PropertyDouble0D)
                            
                            call HDF5SetLimits (obj_hdf5, 1, 1, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                            
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array1D=PtrDouble,                &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif                            
                            item%Values = PtrDouble(1)
                            
                        type is (T_PropertyDouble1D)
                    
                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array1D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                        type is (T_PropertyDouble2D)
                    
                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB,                         &
                                                item%JLB, item%JUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array2D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                        type is (T_PropertyDouble3D)

                            call HDF5SetLimits (obj_hdf5,                                   &
                                                item%ILB, item%IUB,                         &
                                                item%JLB, item%JUB,                         &
                                                item%KLB, item%KUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                            call HDF5ReadData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), Array3D=item%Values,              &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Load = .false.
                                return
                            endif
                    
                        end select
                    endif
                end select
                item => this%pProperties%next ()  
            enddo do1
            
            call KillHDF5(obj_hdf5, STAT=stat)
            if (stat /= SUCCESS_) then
                Load = .false.
                return
            endif            
            
            Load = .true.
        else
            Load = .false.
        endif                
    
    end function Load
    
    !--------------------------------------------------------------------------
    
    function Save (this, now, restart_file)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        type(T_Time)                                    ::  now
        logical, optional                               ::  restart_file
        logical                                         ::  Save
        
        !Local-----------------------------------------------------------------        
        logical                                         ::  restart_file_
        integer                                         ::  obj_hdf5,   &
                                                            HDF5_CREATE,&
                                                            stat
        class(*), pointer                               ::  item
        character(len=PathLength)                       ::  file
        real, dimension(6), target                      ::  AuxTime
        real, dimension(:), pointer                     ::  AuxTimePointer
        integer, dimension(1), target                   ::  AuxInteger
        integer, dimension(:), pointer                  ::  PtrInteger
        real(4), dimension(1), target                   ::  AuxReal
        real(4), dimension(:), pointer                  ::  PtrReal
        real(8), dimension(1), target                   ::  AuxDouble
        real(8), dimension(:), pointer                  ::  PtrDouble        
        
        !----------------------------------------------------------------------
        Save = .true.
        
        if (present(restart_file)) then
            restart_file_ = restart_file
        else
            restart_file = .false.
        endif
        
        if (.not. this%pProperties%IsEmpty()) then
            
            PtrInteger  => AuxInteger
            PtrReal     => AuxReal
            PtrDouble   => AuxDouble
            
            if (restart_file_) then
                file = ChangeSuffix(trim(this%File), "_"//trim(TimeToString(now))//".fin")
            else
                file = trim(this%File)
            endif
            
            obj_hdf5 = 0
            call GetHDF5FileAccess(HDF5_CREATE = HDF5_CREATE)
            call ConstructHDF5(obj_hdf5, trim(file), HDF5_CREATE, STAT=stat)
            if (stat /= SUCCESS_) then
                Save = .false.
                return
            endif
            
            call ExtractDate(now, AuxTime(1), AuxTime(2), AuxTime(3),                   &
                                  AuxTime(4), AuxTime(5), auxTime(6))
            AuxTimePointer => AuxTime
            
            if (this%WriteHorizGrid) then
                call WriteHorizontalGrid(this%ObjHorizGrid, obj_hdf5, STAT=stat)
                if (stat /= SUCCESS_) then
                    Save = .false.
                    this%pMessage = 'Error during Horizontal Grid recording to HDF file'
                    return
                endif
            endif
            
            call HDF5SetLimits (obj_hdf5, 1, 6, STAT=stat)
            if (stat /= SUCCESS_) then
                Save = .false.
                return
            endif
            
            call HDF5WriteData(obj_hdf5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                               Array1D=AuxTimePointer, STAT=stat)
            if (stat /= SUCCESS_) then
                Save = .false.
                return
            endif
            
            item => this%pProperties%First ()
do1:        do while (associated (item))
                select type (item)
                class is (T_Property)
                    if (item%Save) then
                        select type (item)
                        type is (T_PropertyInteger0D)
                    
                            call HDF5SetLimits (obj_hdf5, 1, 1, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            PtrInteger(1) = item%Values
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array1D=PtrInteger,       &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                            
                        type is (T_PropertyInteger1D)
                    
                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array1D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        type is (T_PropertyInteger2D)
                    
                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB,                             &
                                                item%JLB, item%JUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array2D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        type is (T_PropertyInteger3D)
                    
                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB,                             &
                                                item%JLB, item%JUB,                             &
                                                item%KLB, item%KUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array3D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        type is (T_PropertyReal0D)
                    
                            call HDF5SetLimits (obj_hdf5, 1, 1, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            PtrReal(1) = item%Values
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array1D=PtrReal,          &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                            
                        type is (T_PropertyReal1D)
                    
                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array1D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        type is (T_PropertyReal2D)
                    
                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB,                             &
                                                item%JLB, item%JUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array2D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        type is (T_PropertyReal3D)
                    
                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB,                             &
                                                item%JLB, item%JUB,                             &
                                                item%KLB, item%KUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array3D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        type is (T_PropertyDouble0D)
                    
                            call HDF5SetLimits (obj_hdf5, 1, 1, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            PtrDouble(1) = item%Values
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array1D=PtrDouble,        &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                            
                        type is (T_PropertyDouble1D)
                    
                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array1D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        type is (T_PropertyDouble2D)
                    
                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB,                             &
                                                item%JLB, item%JUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array2D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        type is (T_PropertyDouble3D)

                            call HDF5SetLimits (obj_hdf5,                                       &
                                                item%ILB, item%IUB,                             &
                                                item%JLB, item%JUB,                             &
                                                item%KLB, item%KUB, STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                            call HDF5WriteData (obj_hdf5, "/"//trim(item%Group),      &
                                               trim(item%Name), item%Units, Array3D=item%Values,      &
                                               STAT=stat)
                            if (stat /= SUCCESS_) then
                                Save = .false.
                                return
                            endif
                    
                        end select
                    endif
                end select
                item => this%pProperties%next ()  
            enddo do1            
            
            call HDF5FlushMemory (obj_hdf5, STAT=stat)
            if (stat /= SUCCESS_) then
                Save = .false.
                return
            endif
            
            call KillHDF5(obj_hdf5, STAT=stat)
            if (stat /= SUCCESS_) then
                Save = .false.
                return
            endif
        else
            Save = .false.
        endif
        
    end function Save
    
    !--------------------------------------------------------------------------
    
    subroutine Remove (this, id)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        integer                                         ::  id
        
        !Local-----------------------------------------------------------------
        class (T_Property), pointer                     ::  item
        
        !----------------------------------------------------------------------
        
        if (.not. this%pProperties%IsEmpty()) then
            
            item => this%pProperties%First ()
do1:        do while (associated (item))
                if (item%Id == id) exit do1
                item => this%pProperties%Next ()
            enddo do1
            
            if (associated (item)) &
                call this%pProperties%Remove ()
        endif        
    
    end subroutine Remove
    
    !--------------------------------------------------------------------------
    
    subroutine Clear (this)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        
        !----------------------------------------------------------------------
        call this%pProperties%Clear ()
    
    end subroutine Clear
    
    !--------------------------------------------------------------------------
    
    function Message (this)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        character (StringLength)                        ::  Message
        
        !----------------------------------------------------------------------    
        Message = this%pMessage
    
    end function Message
    
    !--------------------------------------------------------------------------
    
    subroutine InitializeObject (this, stop_on_wrong_date, write_horiz_grid, obj_horiz_grid, file)
    
        !Arguments-------------------------------------------------------------
        class(T_StartAndStop)                           ::  this
        logical                                         ::  stop_on_wrong_date
        logical                                         ::  write_horiz_grid
        integer                                         ::  obj_horiz_grid
        character (*)                                   ::  file
        
        !----------------------------------------------------------------------
        this%StopOnWrongDate = stop_on_wrong_date
        this%File             = file
        this%pProperties      = C_Collection (own_items = .true.)
        this%WriteHorizGrid   = write_horiz_grid
        this%ObjHorizGrid     = obj_horiz_grid
        
    end subroutine InitializeObject
    
    !--------------------------------------------------------------------------
    
    subroutine KillObject (this)
    
        !Arguments-------------------------------------------------------------
        type(T_StartAndStop)                            ::  this
        
        !----------------------------------------------------------------------
        call this%Clear ()
        !deallocate (this%pProperties)
        
    end subroutine KillObject
    
    !--------------------------------------------------------------------------
    
end module ModuleStartAndStop   