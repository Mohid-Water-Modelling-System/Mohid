!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model - Results consolidation
! PROJECT       : DDC - Domain Decomposition Consolidation
! PROGRAM       : MainDDC
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          :  2013
! REVISION      : Ricardo Miranda
! BASED ON WORK : Paul Hiemstra
! DESCRIPTION   : Program to consolidate results from MPI run with domain decomposition
!
!
!------------------------------------------------------------------------------


MODULE ModuleHashTable

    use ModuleGlobalData

    implicit none
    private

    public  :: hash_init

    public  :: hash_count
    private ::    hash_count2

    public  :: hash_get
    private ::      hash_get2
    public  :: hash_get_last
    private ::      hash_get_last2
    public  :: hash_get_first
    public  :: hash_get_first_exists
    public  :: hash_get_first_key
    public  :: hash_get_next
    public  :: hash_get_next_exists
    public  :: hash_get_next_key
    private ::      hash_get_next_key2

    public  :: hash_set
    private ::      hash_set2
    private ::      hash_push

    public  :: hash_pop
    private ::      hash_pop2

    public  :: KillHash_map

    public  :: T_HashTable
    type       T_HashTable
        private
        type(T_HashList), pointer           :: HashList
    end type  T_HashTable

    private :: T_HashList
    type       T_HashList
        type(T_HashList), pointer           :: Next

        CHARACTER(PathLength)               :: key
        INTEGER                             :: value_

        integer, dimension(4)               :: DomainSize       = NULL_INT
        integer, dimension(4)               :: WindowPosition   = NULL_INT
        integer, dimension(4)               :: WindowFrame      = NULL_INT
        integer                             :: ObjID            = NULL_INT
    end type  T_HashList

    contains


    function hash_init()

        !Function----------------------------------------------------------------
        type (T_HashTable), pointer                     :: hash_init

        !Local-------------------------------------------------------------------
        type (T_HashTable), pointer                     :: NewObjHashTable

        !------------------------------------------------------------------------

        allocate(NewObjHashTable)
        nullify (NewObjHashTable%HashList)

        hash_init => NewObjHashTable

        !------------------------------------------------------------------------

    END function hash_init

    !--------------------------------------------------------------------------

    integer function hash_count(Me)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me

        !------------------------------------------------------------------------

if1 :   if (associated(Me)) then
if2 :       if (associated(Me%HashList)) then
                hash_count = hash_count2(Me%HashList, count_ = 1)
            else if2
                hash_count = 0
            endif if2
        else if1
            hash_count = 0
        endif if1

        !------------------------------------------------------------------------

    end function hash_count

    !--------------------------------------------------------------------------

    recursive integer function hash_count2(HashList, count_)

        !External----------------------------------------------------------------
        type (T_HashList), pointer                      :: HashList
        integer                                         :: count_

        !------------------------------------------------------------------------

if2 :   if (associated(HashList%Next)) then
            hash_count2 = hash_count2(HashList%Next, count_ = count_+1)
        else if2
            hash_count2 = count_
        endif if2

        !------------------------------------------------------------------------

    end function hash_count2

    !--------------------------------------------------------------------------

    function hash_push(key, value_)

        !Function----------------------------------------------------------------
        type(T_HashList), pointer          :: hash_push

        !External----------------------------------------------------------------
        CHARACTER(*), INTENT(IN)           :: key
        INTEGER     , INTENT(IN)           :: value_

        !Local-------------------------------------------------------------------
        type(T_HashList), pointer          :: NewHashList

        !------------------------------------------------------------------------

        allocate(NewHashList)
        nullify (NewHashList%Next)

        NewHashList%key     = adjustl(trim(key))
        NewHashList%value_  = value_

        hash_push => NewHashList

        !------------------------------------------------------------------------

    end function hash_push

    !--------------------------------------------------------------------------

    subroutine hash_set(Me, key, value_)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer :: Me
        character(*), intent(IN)    :: key
        integer, optional           :: value_
        integer                     :: value_l

        !------------------------------------------------------------------------

if1 :   if (.NOT. present(value_)) then
            value_l = hash_count(Me) + 1
        else if1
            value_l = value_
        end if if1

        call hash_set2(Me%HashList, key    = key,                               &
                                    value_ = value_l)

        !------------------------------------------------------------------------

    end subroutine hash_set

    !--------------------------------------------------------------------------

    recursive SUBROUTINE hash_set2(HashList, key, value_)

        !External----------------------------------------------------------------
        type (T_HashList), pointer                      :: HashList
        CHARACTER(*), INTENT(IN)                        :: key
        INTEGER     , INTENT(IN)                        :: value_

        !------------------------------------------------------------------------

if2 :   if (.NOT. associated(HashList)) then
            HashList => hash_push(key    = key,                                 &
                                  value_ = value_)
        else if2
if1 :       if (TRIM(HashList%key) == TRIM(key)) then
                HashList%value_ = value_
            else if1
                call hash_set2(HashList%Next, key    = key,                     &
                                              value_ = value_)
            endif if1
        endif if2

        !------------------------------------------------------------------------

    END SUBROUTINE hash_set2

    !--------------------------------------------------------------------------

    integer function hash_get_last(Me)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList)) then
            hash_get_last = hash_get_last2(Me%HashList)
        else if2
            hash_get_last =-1 * NOT_FOUND_ERR_
        endif if2

        !------------------------------------------------------------------------

    end function hash_get_last

    !--------------------------------------------------------------------------

    integer recursive function hash_get_last2(HashList)

        !External----------------------------------------------------------------
        type (T_HashList), pointer                      :: HashList

        !------------------------------------------------------------------------

if2 :   if (associated(HashList%Next)) then
            hash_get_last2 = hash_get_last2(HashList%Next)
        else if2
            hash_get_last2 = HashList%value_
        endif if2

        !------------------------------------------------------------------------

    end function hash_get_last2

    !--------------------------------------------------------------------------

    integer function hash_get(Me, key)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList)) then
            hash_get = hash_get2(Me%HashList, key = key)
        else if2
            hash_get =-1 * NOT_FOUND_ERR_
        endif if2

        !------------------------------------------------------------------------

    END function hash_get

    !--------------------------------------------------------------------------

    recursive integer function hash_get2(HashList, key)

        !External----------------------------------------------------------------
        type (T_HashList), pointer                      :: HashList
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if1 :   if (adjustl(trim(HashList%key)) .EQ. adjustl(trim(key))) then
            hash_get2 = HashList%value_
        else if1
if2 :       if (associated(HashList%Next)) then
                hash_get2 = hash_get2(HashList%Next, key = key)
            else if2
                hash_get2 =-1 * NOT_FOUND_ERR_
            endif if2
        endif if1

        !------------------------------------------------------------------------

    END function hash_get2

    !--------------------------------------------------------------------------

    logical function hash_get_first_exists(Me)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList)) then
            hash_get_first_exists = .TRUE.
        else if2
            hash_get_first_exists =.FALSE.
        endif if2

        !------------------------------------------------------------------------

    END function hash_get_first_exists
    !--------------------------------------------------------------------------

    integer function hash_get_first(Me)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList)) then
            hash_get_first = Me%HashList%value_
        else if2
            hash_get_first =-1 * NOT_FOUND_ERR_
        endif if2

        !------------------------------------------------------------------------

    END function hash_get_first

    !--------------------------------------------------------------------------

    function hash_get_first_key(Me)

        !Function----------------------------------------------------------------
        CHARACTER(PathLength)                           :: hash_get_first_key

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList)) then
            hash_get_first_key = adjustl(trim(Me%HashList%key))
        else if2
            hash_get_first_key = NULL_STR
        endif if2

        !------------------------------------------------------------------------

    END function hash_get_first_key

    !--------------------------------------------------------------------------

    integer function hash_get_next(Me, key)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList%Next)) then
            hash_get_next = hash_get2(Me%HashList%Next, key = key)
        else if2
            hash_get_next =-1 * NOT_FOUND_ERR_
        endif if2

        !------------------------------------------------------------------------

    END function hash_get_next

    !--------------------------------------------------------------------------

    logical function hash_get_next_exists(Me, key)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList)) then
            hash_get_next_exists = hash_get_next_exists2(Me%HashList, key)
        else if2
            hash_get_next_exists =.FALSE.
        endif if2

        !------------------------------------------------------------------------

    END function hash_get_next_exists

    !--------------------------------------------------------------------------

    logical recursive function hash_get_next_exists2(HashList, key)

        !External----------------------------------------------------------------
        type (T_HashList), pointer                      :: HashList
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if1 :   if (adjustl(trim(HashList%key)) .EQ. adjustl(trim(key))) then
if3 :       if (associated(HashList%Next)) then
                hash_get_next_exists2 = .TRUE.
            else if3
                hash_get_next_exists2 = .FALSE.
            endif if3
        else if1
if4 :       if (associated(HashList%Next)) then
                hash_get_next_exists2 = hash_get_next_exists2(HashList%Next, key)
            else if4
                hash_get_next_exists2 = .FALSE.
            endif if4
        endif if1

        !------------------------------------------------------------------------

    END function hash_get_next_exists2

    !--------------------------------------------------------------------------

    function hash_get_next_key(Me, key)

        !Function----------------------------------------------------------------
        CHARACTER(PathLength)                           :: hash_get_next_key

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList)) then
            hash_get_next_key = hash_get_next_key2(Me%HashList, key)
        else if2
            hash_get_next_key = NULL_STR
        endif if2

        !------------------------------------------------------------------------

    END function hash_get_next_key

    !--------------------------------------------------------------------------

    recursive function hash_get_next_key2(HashList, key)

        !Function----------------------------------------------------------------
        CHARACTER(PathLength)                           :: hash_get_next_key2

        !External----------------------------------------------------------------
        type (T_HashList), pointer                      :: HashList
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if1 :   if (adjustl(trim(HashList%key)) .EQ. adjustl(trim(key))) then
            hash_get_next_key2 = HashList%Next%Key
        else if1
if2 :       if (associated(HashList%Next)) then
                hash_get_next_key2 = hash_get_next_key2(HashList%Next, key = key)
            else if2
                hash_get_next_key2 = NULL_STR
            endif if2
        endif if1

        !------------------------------------------------------------------------

    END function hash_get_next_key2

    !--------------------------------------------------------------------------

    integer function hash_pop(Me, key)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if2 :   if (associated(Me%HashList)) then
            hash_pop = hash_pop2(Me%HashList, key = key)
        else if2
            hash_pop =-1 * NOT_FOUND_ERR_
        endif if2

        !------------------------------------------------------------------------

    END function hash_pop

    !--------------------------------------------------------------------------

    recursive integer function hash_pop2(HashList, key)

        !External----------------------------------------------------------------
        type (T_HashList), pointer                      :: HashList
        CHARACTER(*), INTENT(IN)                        :: key

        !------------------------------------------------------------------------

if1 :   if (adjustl(trim(HashList%key)) .EQ. adjustl(trim(key))) then
            HashList => HashList%Next
            hash_pop2 = SUCCESS_
        else if1
if2 :       if (associated(HashList%Next)) then
                hash_pop2 = hash_pop2(HashList%Next, key = key)
            else if2
                hash_pop2 =-1 * NOT_FOUND_ERR_
            endif if2
        endif if1

        !------------------------------------------------------------------------

    END function hash_pop2

    !--------------------------------------------------------------------------

    SUBROUTINE KillHash_map(Me)

        !External----------------------------------------------------------------
        type (T_HashTable), pointer                     :: Me

        !------------------------------------------------------------------------

        if (associated(Me%HashList)) call deallocateHashList(Me%HashList)

        deallocate (Me)
        nullify    (Me)

        !------------------------------------------------------------------------

    END SUBROUTINE KillHash_map

    !--------------------------------------------------------------------------

    recursive SUBROUTINE deallocateHashList(HashList)

        !External----------------------------------------------------------------
        type (T_HashList), pointer                      :: HashList

        !------------------------------------------------------------------------

        if (associated(HashList%Next)) call deallocateHashList(HashList%Next)

        DEALLOCATE(HashList)
        nullify   (HashList)

        !------------------------------------------------------------------------

    END SUBROUTINE deallocateHashList

    !--------------------------------------------------------------------------

END MODULE ModuleHashTable
