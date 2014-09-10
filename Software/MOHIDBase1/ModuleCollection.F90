!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Snow
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jul 2014
! REVISION      : Eduardo Jauch - v4.0
! DESCRIPTION   : Module the provides a linked list class
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
   
module ModuleCollection

    use ModuleGlobalData    
   
    implicit none
    private 
   
    !Public Elements-----------------------------------------------------------        
    public  ::  C_Collection
    
    !Types---------------------------------------------------------------------
    type T_CollectionItem
        type (T_CollectionItem), pointer                ::  Next => null(),         &
                                                            Prior => null()
        class (*), pointer                              ::  Data => null()
    end type T_CollectionItem
    
    type C_Collection
        private            
            type (T_CollectionItem), pointer            ::  pFirstItem => null(),   &
                                                            pLastItem => null(),    &
                                                            pCurrentItem => null()
            
            integer                                     ::  pCount = 0,             &
                                                            pCurrentIndex = -1
            
            logical, public                             ::  OwnItems = .false.
            
        contains
            procedure                                   ::  Add
            procedure, private                          ::  RemoveCurrent
            procedure, private                          ::  RemoveByIndex
            generic                                     ::  Remove =>               &
                                                            RemoveCurrent,          &
                                                            RemoveByIndex
            procedure                                   ::  Clear
            procedure                                   ::  Insert
            
            procedure                                   ::  GoTo
            procedure                                   ::  GoToNext
            procedure                                   ::  GoToPrior
            procedure                                   ::  GoToFirst
            procedure                                   ::  GoToLast
            procedure                                   ::  Next
            procedure                                   ::  Prior
            procedure                                   ::  First
            procedure                                   ::  Last
            procedure, private                          ::  CheckIndex
            
            procedure                                   ::  HasNext
            procedure                                   ::  HasPrior
            procedure                                   ::  Count
            procedure                                   ::  IsEmpty                        
            
            procedure                                   ::  Item
            procedure                                   ::  ItemIndex
            
            procedure, private                          ::  InitializeObject            
            final                                       ::  KillObject
    
    end type C_Collection
        
    !Interfaces----------------------------------------------------------------
    
    interface C_Collection
        module procedure Constructor
    end interface C_Collection
    
    !Subroutines--------------------------------------------------------------- 
    contains
    
    function Constructor (own_items)
    
        !Arguments-------------------------------------------------------------
        type(C_Collection)                              :: Constructor
        logical, optional, intent (in)                  :: own_items
        
        !----------------------------------------------------------------------
        if (present (own_items)) then
            call constructor%InitializeObject (own_items)
        else
            call constructor%InitializeObject (.false.)
        endif        
        
    end function
    
    !--------------------------------------------------------------------------
    
    subroutine Add (this, item)
    
        !Arguments-------------------------------------------------------------
        class (C_Collection)                            ::  this
        class (*), pointer, intent (in)                 ::  item        
        
        !Local-----------------------------------------------------------------
        type (T_CollectionItem), pointer                ::  new_collection_item
        
        !----------------------------------------------------------------------        

        allocate (new_collection_item)        
        new_collection_item%Data => item
        
        this%pCount = this%pCount + 1

        if (.not. associated(this%pLastItem)) then
           
            new_collection_item%Next  => null()
            new_collection_item%Prior => null()
            
            this%pFirstItem => new_collection_item
            this%pLastItem  => new_collection_item
            
        else
            
            new_collection_item%Next  => null()
            new_collection_item%Prior => this%pLastItem
            
            this%pLastItem%Next => new_collection_item
            this%pLastitem      => new_collection_item
            
        endif
        
    end subroutine Add
    
    !--------------------------------------------------------------------------
    
    subroutine RemoveCurrent (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        
        !Local-----------------------------------------------------------------
        type (T_CollectionItem), pointer                ::  temp
        
        !----------------------------------------------------------------------                
        if (associated (this%pCurrentItem)) then
            
            temp => this%pCurrentItem
            nullify (this%pCurrentItem)
            
            if (associated (temp%Next)) then                
                temp%Next%Prior   => temp%Prior
                this%pCurrentItem => temp%Next                
            endif
            
            if (associated (temp%Prior)) then                
                temp%Prior%Next => temp%Next
                
                if (.not. associated (this%pCurrentItem)) then                    
                    this%pCurrentItem => temp%Prior                    
                endif                
            endif
            
            if (this%OwnItems) then
                deallocate (temp%Data)
            endif
            
            deallocate (temp)
            this%pCount = this%pCount - 1
            
            if (this%pCurrentIndex > this%pCount) then                
                this%pCurrentIndex = this%pCount
                this%pLastItem => this%pCurrentItem
                
                if (this%pCount == 1) &
                    this%pFirstItem => this%pCurrentItem
            endif
            
            if (this%pCurrentIndex <= 0) then
                nullify (this%pFirstItem)
                nullify (this%pLastItem)
            endif
            
        endif
        
    end subroutine RemoveCurrent

    !--------------------------------------------------------------------------
    
    subroutine RemoveByIndex (this, index)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        integer                                         ::  index
        
        !Local-----------------------------------------------------------------
        integer                                         ::  old_index
        
        !----------------------------------------------------------------------        
        if (this%CheckIndex (index)) then
            
            old_index = this%pCurrentIndex
            call this%GoTo (index)
            call this%Remove
            
            if (old_index > this%pCount) then
                call this%GoToLast            
            elseif (old_index > this%pCurrentIndex) then
                call this%GoTo (old_index - 1)
            else
                call this%GoTo (old_index)
            endif
        
        endif
        
    end subroutine RemoveByIndex
    
    !--------------------------------------------------------------------------
    
    subroutine Clear (this)

        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        
        !Local-----------------------------------------------------------------
        call this%GoToLast
        do while (.not. this%IsEmpty())
            call this%Remove
        enddo
        
        !----------------------------------------------------------------------
        
    end subroutine Clear
    
    !--------------------------------------------------------------------------
    
    subroutine Insert (this, item, index)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        class(*), pointer                               ::  item
        integer                                         ::  index
        
        !Local-----------------------------------------------------------------
        type (T_CollectionItem), pointer                ::  new_collection_item
        integer                                         ::  old_index
        
        !----------------------------------------------------------------------
        if (this%CheckIndex (index)) then
            
            if (index <= this%pCurrentIndex) then
                old_index = this%pCurrentIndex + 1
            else
                old_index = this%pCurrentIndex
            endif                        
            
            allocate (new_collection_item)
            new_collection_item%Data => item
        
            call this%GoTo (index)
            
            if (this%HasPrior()) then
                this%pCurrentItem%Prior%Next => new_collection_item
                new_collection_item%Prior => this%pCurrentItem%Prior
            endif
            
            if (this%HasNext()) then
                this%pCurrentItem%Next%Prior => new_collection_item
                new_collection_item%Next => this%pCurrentItem%Next
            endif                        
            
            this%pCount = this%pCount + 1                                    
            
            call this%GoTo (old_index)                        
            
        endif
        
    end subroutine Insert
    
    !--------------------------------------------------------------------------
    
    subroutine GoTo (this, index)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        integer                                         ::  index
        
        !----------------------------------------------------------------------
        if (this%CheckIndex (index)) then
            call this%GoToFirst
            do while (this%ItemIndex() /= index)
                call this%GoToNext
            enddo
        endif
        
    end subroutine GoTo
    
    !--------------------------------------------------------------------------
    
    subroutine GoToNext (this)
        
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
                
        !----------------------------------------------------------------------        
        if (this%HasNext()) then
            
            this%pCurrentItem => this%pCurrentItem%Next
            this%pCurrentIndex = this%pCurrentIndex + 1
            
        endif
        
    end subroutine GoToNext
    
    !--------------------------------------------------------------------------
    
    subroutine GoToPrior (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        
        !----------------------------------------------------------------------
        if (this%HasPrior()) then
            
            this%pCurrentItem => this%pCurrentItem%Prior
            this%pCurrentIndex = this%pCurrentIndex - 1
            
        endif
        
    end subroutine GoToPrior

    !--------------------------------------------------------------------------
    
    subroutine GoToFirst (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        
        !----------------------------------------------------------------------
        if (this%pCount > 0) then
            this%pCurrentItem => this%pFirstItem
            this%pCurrentIndex = 1
        endif
    
    end subroutine GoToFirst

    !--------------------------------------------------------------------------
    
    subroutine GoToLast (this)
        
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        
        !----------------------------------------------------------------------
        if (this%pCount > 0) then
            this%pCurrentItem => this%pLastItem
            this%pCurrentIndex = this%pCount
        endif
        
    end subroutine GoToLast

    !--------------------------------------------------------------------------
    
    function Next (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        class(*), pointer                               ::  Next
        
        !Local-----------------------------------------------------------------
        integer                                         ::  old_index
        
        !----------------------------------------------------------------------
        
        old_index = this%pCurrentIndex        
        call this%GoToNext
        
        if (old_index < this%pCurrentIndex) then
            Next => this%pCurrentItem%Data
        else
            Next => null()
        endif
        
    end function Next
    
    !--------------------------------------------------------------------------
    
    function Prior (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        class(*), pointer                               ::  Prior
        
        !Local-----------------------------------------------------------------
        integer                                         ::  old_index
        
        !----------------------------------------------------------------------
        
        old_index = this%pCurrentIndex        
        call this%GoToPrior
        
        if (old_index > this%pCurrentIndex) then
            Prior => this%pCurrentItem%Data
        else
            Prior => null()
        endif
        
    end function Prior
    
    !--------------------------------------------------------------------------
    
    function First (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        class(*), pointer                               ::  First
        
        !Local-----------------------------------------------------------------        
        
        !----------------------------------------------------------------------
               
        call this%GoToFirst
        
        if (this%pCurrentIndex == 1) then
            First => this%pCurrentItem%Data
        else
            First => null()
        endif
        
    end function First
    
    !--------------------------------------------------------------------------
    
    function Last (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        class(*), pointer                               ::  Last
        
        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------
               
        call this%GoToLast
        
        if (this%pCurrentIndex >= 1 .and. this%pCurrentIndex == this%pCount) then
            Last => this%pCurrentItem%Data
        else
            Last => null()
        endif
        
    end function Last
    
    !--------------------------------------------------------------------------
    
    function CheckIndex (this, index)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        integer                                         ::  index
        logical                                         ::  CheckIndex
        
        !----------------------------------------------------------------------
        if (index > 0 .and. index <= this%pCount) then
            CheckIndex = .true.
        else
            CheckIndex = .false.
        endif
        
    end function CheckIndex
    
    !--------------------------------------------------------------------------
    
    function HasNext (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        logical                                         ::  HasNext
        
        !----------------------------------------------------------------------
                
        HasNext = .false.
        
        if (associated (this%pCurrentItem)) then
            if (associated (this%pCurrentItem%Next)) then
                HasNext = .true.
            endif
        endif
        
    end function HasNext

    !--------------------------------------------------------------------------
    
    function HasPrior (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        logical                                         ::  HasPrior
        
        !----------------------------------------------------------------------
                
        HasPrior = .false.
        
        if (associated (this%pCurrentItem)) then
            if (associated (this%pCurrentItem%Prior)) then
                HasPrior = .true.
            endif
        endif
        
    end function HasPrior
    
    !--------------------------------------------------------------------------
    
    function Count (this)
        
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        integer                                         ::  Count
        
        !----------------------------------------------------------------------
        Count = this%pCount
        
    end function Count
    
    !--------------------------------------------------------------------------
    
    function IsEmpty (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        logical                                         ::  IsEmpty
        
        !----------------------------------------------------------------------
        if (this%pCount > 0) then
            IsEmpty = .false.
        else
            IsEmpty = .true.
        endif
        
    end function IsEmpty

    !--------------------------------------------------------------------------
    
    function Item (this)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        class(*), pointer                               ::  Item
        
        !----------------------------------------------------------------------
        if (associated(this%pCurrentItem)) then
            Item => this%pCurrentItem%Data
        else
            Item => null()
        endif
        
    end function Item
    
    !--------------------------------------------------------------------------
    
    function ItemIndex (this)
        
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this        
        integer                                         ::  ItemIndex
        
        !----------------------------------------------------------------------
        ItemIndex = this%pCurrentIndex

    end function ItemIndex
    
    !--------------------------------------------------------------------------
    
    subroutine InitializeObject (this, own_items)
    
        !Arguments-------------------------------------------------------------
        class(C_Collection)                             ::  this
        logical, intent(in)                             ::  own_items
        
        !----------------------------------------------------------------------
        this%OwnItems = own_items
        
    end subroutine InitializeObject
    
    !--------------------------------------------------------------------------
    
    subroutine KillObject (this)
    
        !Arguments-------------------------------------------------------------
        type(C_Collection)                              ::  this
        
        !----------------------------------------------------------------------
        call this%Clear
        
    end subroutine KillObject
    
    !--------------------------------------------------------------------------
    
end module ModuleCollection