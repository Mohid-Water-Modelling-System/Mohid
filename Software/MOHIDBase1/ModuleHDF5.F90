!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : HDF 5
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to Write/Read HDF 5 Files
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

Module ModuleHDF5

    use ModuleGlobalData
    use ModuleFunctions     , only :  minival, maxival, SetMatrixValue
    use hdf5

#ifdef _GUI_    
    use dfwin
    use comctl32
#endif

    implicit none 

    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  ::  ConstructHDF5
    private ::      AllocateInstance
    
    !Modifier
    public  ::  HDF5CreateGroup
    private ::      CreateMinMaxAttribute
    private ::      CheckGroupExistence
    public  ::  HDF5SetLimits
    public  ::  HDF5WriteData
    private ::      PrepareWrite
    private ::      FinishWrite
    private ::      UpdateMinMaxAttribute
    private ::      ConstructDSName
    public  ::  HDF5WriteGlobalAttribute
    public  ::  HDF5ReadGlobalAttribute
    public  ::  HDF5ReadData
    public  ::  HDF5ReadHyperSlab
    public  ::  HDF5ReadWindow 
    public  ::  HDF5WriteWindow   
    public  ::  HDF5FlushMemory

    !Destructor
    public  ::  KillHDF5
    private ::      DeallocateInstance

    !Selector
    public  ::  GetHDF5FileID
    public  ::  GetHDF5FileAccess
    public  ::  GetHDF5GroupID
    public  ::  GetHDF5GroupNumberOfItems
    public  ::  GetHDF5GroupExist
    public  ::  GetHDF5DataSetExist
    public  ::  GetHDF5ObjectInfo
    public  ::  GetHDF5ArrayDimensions
    public  ::  GetHDF5ArrayDim
    public  ::  GetHDF5FileOkToRead
    public  ::  GetHDF5DataTypeID
    public  ::  GetHDF5FileName

#ifdef _GUI_
    public  :: HDF5InquireFile
    private ::      InquireSubGroup
    private ::      AddItemToHDFTree
    public  :: ConstructTreeViewImageLists
    public  :: HDF5ReadAttributes
    public  :: HDF5GetDimensions
#endif

    public  :: HDF5ReadGenericRealAttribute
    public  :: HDF5UpdateGenericRealAttribute
    
    !Management
    private ::  Ready
    private ::      LocateObjHDF5

    interface HDF5WriteData
        module procedure HDF5WriteDataR4_1D
        module procedure HDF5WriteDataR4_2D
        module procedure HDF5WriteDataR4_3D
        module procedure HDF5WriteDataR8_1D
        module procedure HDF5WriteDataR8_2D
        module procedure HDF5WriteDataR8_3D
        module procedure HDF5WriteDataI4_1D
        module procedure HDF5WriteDataI4_2D
        module procedure HDF5WriteDataI4_3D
    end interface

    interface HDF5ReadData
        module procedure HDF5ReadDataR4_1D
        module procedure HDF5ReadDataR4_2D
        module procedure HDF5ReadDataR4_3D
        module procedure HDF5ReadDataR8_1D
        module procedure HDF5ReadDataR8_2D
        module procedure HDF5ReadDataR8_3D
        module procedure HDF5ReadDataI4_1D
        module procedure HDF5ReadDataI4_2D
        module procedure HDF5ReadDataI4_3D
    end interface
    
    interface HDF5ReadWindow    
        module procedure HDF5ReadWindowR4_3D
        module procedure HDF5ReadWindowR4_2D
        module procedure HDF5ReadWindowR4_1D
        module procedure HDF5ReadWindowR8_3D
        module procedure HDF5ReadWindowR8_2D
        module procedure HDF5ReadWindowR8_1D
        module procedure HDF5ReadWindowI4_3D
        module procedure HDF5ReadWindowI4_2D
        module procedure HDF5ReadWindowI4_1D
    end interface
    
    interface HDF5WriteWindow
        module procedure HDF5WriteWindowR4_3D    
        module procedure HDF5WriteWindowR4_2D
        module procedure HDF5WriteWindowR4_1D
        module procedure HDF5WriteWindowR8_3D
        module procedure HDF5WriteWindowR8_2D
        module procedure HDF5WriteWindowR8_1D
        module procedure HDF5WriteWindowI4_3D
        module procedure HDF5WriteWindowI4_2D
        module procedure HDF5WriteWindowI4_1D
    end interface

    interface HDF5WriteGlobalAttribute
        module procedure HDF5WriteGlobalAttribute_Char
        module procedure HDF5WriteGlobalAttribute_Int
        module procedure HDF5WriteGlobalAttribute_R4
        module procedure HDF5WriteGlobalAttribute_R8   
    end interface

    interface HDF5ReadGlobalAttribute
        module procedure HDF5ReadGlobalAttribute_Char
        module procedure HDF5ReadGlobalAttribute_Int
        module procedure HDF5ReadGlobalAttribute_R4
        module procedure HDF5ReadGlobalAttribute_R8   
    end interface

    !Parameter
    integer, parameter                              :: HDF5_CREATE_     = 1
    integer, parameter                              :: HDF5_READ_       = 2
    integer, parameter                              :: HDF5_READWRITE_  = 3 
    
!    !Type definition
!    type T_Limits
!        integer                                     :: ILB, IUB
!        integer                                     :: JLB, JUB
!        integer                                     :: KLB, KUB
!    end type T_Limits
!
    !GUI
    type T_HDF5_DATA_ITEM
        sequence
        integer                                     :: ItemType
        character(Pathlength)                       :: FileName
        character(StringLength)                     :: GroupName
        character(StringLength)                     :: ItemName
        integer                                     :: Rank
        integer                                     :: NumType
        integer, dimension(7)                       :: Dims
    end type T_HDF5_DATA_ITEM
    
    type T_AuxMatrixes
        real(4),    dimension(:    ), pointer :: DataR4_1D
        real(4),    dimension(:,:  ), pointer :: DataR4_2D
        real(4),    dimension(:,:,:), pointer :: DataR4_3D
        real(8),    dimension(:    ), pointer :: DataR8_1D
        real(8),    dimension(:,:  ), pointer :: DataR8_2D
        real(8),    dimension(:,:,:), pointer :: DataR8_3D
        integer,    dimension(:    ), pointer :: DataI4_1D
        integer,    dimension(:,:  ), pointer :: DataI4_2D
        integer,    dimension(:,:,:), pointer :: DataI4_3D
        
        integer, dimension(7)                 :: dims_R4_1D
        integer, dimension(7)                 :: dims_R4_2D
        integer, dimension(7)                 :: dims_R4_3D
        integer, dimension(7)                 :: dims_R8_1D
        integer, dimension(7)                 :: dims_R8_2D
        integer, dimension(7)                 :: dims_R8_3D
        integer, dimension(7)                 :: dims_I4_1D
        integer, dimension(7)                 :: dims_I4_2D
        integer, dimension(7)                 :: dims_I4_3D
    end type T_AuxMatrixes    
    
    
    type T_HDF5
        integer                                     :: InstanceID
        integer (HID_T)                             :: FileID
        character(Pathlength)                       :: FileName, FileName2
        type    (T_Size3D)                          :: Limits3D, Limits
        type    (T_Size2D)                          :: Limits2D
        type    (T_Size1D)                          :: Limits1D
        type    (T_AuxMatrixes)                     :: AuxMatrixes
        type    (T_HDF5), pointer                   :: Next
    end type T_HDF5

    !Global Module Variables
    type (T_HDF5), pointer                          :: FirstHDF5
    logical                                         :: FirstPassage = .true.
    
    type (T_HDF5), pointer                          :: Me

    !GUI
    integer, parameter                              :: TypeFile = 1
    integer, parameter                              :: TypeVG   = 2
    integer, parameter                              :: TypeSDS  = 3
    integer, parameter                              :: TypeAttr = 4
    integer                                         :: IconFile
    integer                                         :: IconVG
    integer                                         :: IconSDS
    integer                                         :: IconSDS1
    integer                                         :: IconAttr
    integer                                         :: IconVG1

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructHDF5 (HDF5ID, FileName, Access, STAT)
#ifdef _GUI_
        !DEC$ IF DEFINED(_X86_)
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : '_ConstructHDF5@16' :: ConstructHDF5
        !DEC$ ELSE
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : 'ConstructHDF5'     :: ConstructHDF5
        !DEC$ ENDIF
#endif

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        character(len=*)                            :: FileName
        integer                                     :: Access
        integer(4), optional                        :: STAT

        !Local-----------------------------------------------------------------
        integer(4)                                  :: STAT_CALL
        integer                                     :: STAT_, ready_
        
        !Begin-----------------------------------------------------------------


        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (FirstPassage) then
            nullify (FirstHDF5)
            FirstPassage = .false.
        endif

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. OFF_ERR_) then

            !Allocates a new Instance
            call AllocateInstance

            !Initializes predefined datatypes
            call h5open_f (STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5 - ModuleHDF5 - ERR00'

            !Open the file
            if      (Access == HDF5_CREATE_) then
                call h5fcreate_f(trim(FileName), ACCESS_FLAGS = H5F_ACC_TRUNC_F,               &
                                 FILE_ID = Me%FileID, HDFERR = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5 - ModuleHDF5 - ERR01'
            elseif  (Access == HDF5_READ_) then
                call h5fopen_f (trim(FileName), ACCESS_FLAGS = H5F_ACC_RDONLY_F,               &
                                FILE_ID = Me%FileID, HDFERR = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5 - ModuleHDF5 - ERR02'
            elseif  (Access == HDF5_READWRITE_) then
                call h5fopen_f (trim(FileName), ACCESS_FLAGS = H5F_ACC_RDWR_F,                &
                                FILE_ID = Me%FileID, HDFERR = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5 - ModuleHDF5 - ERR03'
            else
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5 - ModuleHDF5 - ERR04'
            endif

            nullify(Me%AuxMatrixes%DataR4_1D, Me%AuxMatrixes%DataR4_2D, Me%AuxMatrixes%DataR4_3D, &
                    Me%AuxMatrixes%DataR8_1D, Me%AuxMatrixes%DataR8_2D, Me%AuxMatrixes%DataR8_3D, &
                    Me%AuxMatrixes%DataI4_1D, Me%AuxMatrixes%DataI4_2D, Me%AuxMatrixes%DataI4_3D)


            !Stores FileName
            Me%FileName  = trim(FileName)
            Me%FileName2 = trim(FileName)
            


            !Returns ID
            HDF5ID     = Me%InstanceID

        else

            stop 'ModuleHDF5 - ConstructHDF5 - ERR99' 

        endif

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructHDF5

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HDF5), pointer                      :: NewHDF5
        type (T_HDF5), pointer                      :: PreviousHDF5


        !Allocates new instance
        allocate (NewHDF5)
        nullify  (NewHDF5%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstHDF5)) then
            FirstHDF5   => NewHDF5
            Me          => NewHDF5
        else
            PreviousHDF5 => FirstHDF5
            Me           => FirstHDF5%Next
            do while (associated(Me))
                PreviousHDF5 => Me
                Me      => Me%Next
            enddo
            Me           => NewHDF5
            PreviousHDF5%Next => NewHDF5
        endif

        Me%InstanceID = RegisterNewInstance (mHDF5_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine HDF5CreateGroup  (HDF5ID, GroupName, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        character(len=*)                            :: GroupName
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_


        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then


            !Recursivly creates all Groups
            call CheckGroupExistence (Me%FileID, GroupName)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5CreateGroup

    !--------------------------------------------------------------------------

    subroutine HDF5SetLimits (HDF5ID, ILB, IUB, JLB, JUB, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        integer, optional                           :: ILB, IUB
        integer, optional                           :: JLB, JUB
        integer, optional                           :: KLB, KUB
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then


            if (present(ILB)) then
                Me%Limits%ILB = ILB
                Me%Limits3D%ILB = ILB
                Me%Limits2D%ILB = ILB
                Me%Limits1D%ILB = ILB                
            endif
            if (present(IUB)) then
                Me%Limits%IUB = IUB
                Me%Limits3D%IUB = IUB
                Me%Limits2D%IUB = IUB
                Me%Limits1D%IUB = IUB
            endif
            if (present(JLB)) then
                Me%Limits%JLB = JLB
                Me%Limits3D%JLB = JLB
                Me%Limits2D%JLB = JLB
            endif
            if (present(JUB)) then
                Me%Limits%JUB = JUB
                Me%Limits3D%JUB = JUB
                Me%Limits2D%JUB = JUB
            endif
            if (present(KLB)) then
                Me%Limits%KLB = KLB
                Me%Limits3D%KLB = KLB
            endif
            if (present(KUB)) then
                Me%Limits%KUB = KUB
                Me%Limits3D%KUB = KUB
            endif

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5SetLimits

    !--------------------------------------------------------------------------

    subroutine CreateMinMaxAttribute (ID, Units, Minimum, Maximum)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: ID
        character(len=*), optional                  :: Units
        real(4), optional                           :: Minimum, Maximum

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1, attr_id2, attr_id3
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: STAT_CALL
        INTEGER(HID_T)                              :: new_type_id

        !Creates data space for Minimum and Maximum attributes
        call h5screate_f (H5S_SCALAR_F, space_id, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR01'

        !Creates attribute Minimum
        call h5acreate_f(ID, "Minimum", H5T_NATIVE_REAL, space_id, attr_id1, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR02'

        !Creates attribute Maximum
        call h5acreate_f(ID, "Maximum", H5T_NATIVE_REAL, space_id, attr_id2, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR03'

        !Writes Miniumum
        if (present(Minimum)) then
            call h5awrite_f (attr_id1, H5T_NATIVE_REAL, Minimum, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR04'
        else
            call h5awrite_f (attr_id1, H5T_NATIVE_REAL, -FillValueReal, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR05'
        endif

        !Writes Maximum
        if (present(Maximum)) then
            call h5awrite_f (attr_id2, H5T_NATIVE_REAL, Maximum, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR06'
        else
            call h5awrite_f (attr_id2, H5T_NATIVE_REAL, FillValueReal, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR07'
        endif

        !Closes Minimum
        call h5aclose_f (attr_id1, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR08'

        !Closes Maximum
        call h5aclose_f (attr_id2, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR09'

        !Closes dataspaces
        call h5sclose_f (space_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR10'


        if (present(Units)) then

            !Creates data space for Units
            call h5screate_f   (H5S_SCALAR_F, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR11'

            !Copies Type
            call h5Tcopy_f     (H5T_NATIVE_CHARACTER, new_type_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR12'

            !Sets Size
            if (len_trim(Units)==0) then
                Units ="-"
            endif
            call h5Tset_size_f (new_type_id, int(len_trim(Units),8), STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR13'

            !Creates attribute
            call h5acreate_f   (ID, "Units", new_type_id, space_id, attr_id3, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR14'

            !Writes attribute
            call h5awrite_f    (attr_id3, new_type_id, Units, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR15'

            !Close type id
            call h5Tclose_f    (new_type_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR16'

            !Closes attribute
            call h5aclose_f    (attr_id3, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR17'

            !Closes dataspaces
            call h5sclose_f    (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CreateMinMaxAttribute - ModuleHDF5 - ERR18'
        endif

              
    end subroutine CreateMinMaxAttribute 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine CreateAverageRadiusAttribute (ID, Average, Radius)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: ID
        real(4)                                     :: Average, Radius

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1, attr_id2
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: STAT_CALL


        !Creates data space for Average and Radius attributes
        call h5screate_f (H5S_SCALAR_F, space_id, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CreateAverageRadiusAttribute - ModuleHDF5 - ERR10'

        !Creates attribute Average
        call h5acreate_f(ID, "Average", H5T_NATIVE_REAL, space_id, attr_id1, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateAverageRadiusAttribute - ModuleHDF5 - ERR20'

        !Creates attribute Radius
        call h5acreate_f(ID, "Radius" , H5T_NATIVE_REAL, space_id, attr_id2, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateAverageRadiusAttribute - ModuleHDF5 - ERR30'

        !Write Average
        call h5awrite_f (attr_id1,      H5T_NATIVE_REAL, Average, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateAverageRadiusAttribute - ModuleHDF5 - ERR40'

        !Writes Radius
        call h5awrite_f (attr_id2,      H5T_NATIVE_REAL, Radius, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateAverageRadiusAttribute - ModuleHDF5 - ERR50'

        !Closes Average
        call h5aclose_f (attr_id1, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateAverageRadiusAttribute - ModuleHDF5 - ERR60'

        !Closes Radius
        call h5aclose_f (attr_id2, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateAverageRadiusAttribute - ModuleHDF5 - ERR70'

        !Closes dataspaces
        call h5sclose_f (space_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateAverageRadiusAttribute - ModuleHDF5 - ERR80'
              
    end subroutine CreateAverageRadiusAttribute 

    !--------------------------------------------------------------------------

    subroutine UpdateMinMaxAttribute (ID, Minimum, Maximum)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: ID
        real(4)                                     :: Minimum, Maximum
        
        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: attr_id
        integer(HID_T)                              :: STAT_CALL
        integer(HSIZE_T), dimension(7)              :: dims
        real                                        :: OldValue = 0

        !Opens Minimum attribute
        call h5aopen_name_f (ID, "Minimum", attr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UpdateMinMaxAttribute - ModuleHDF5 - ERR01'

        !Reads Value
        call h5aread_f      (attr_id, H5T_NATIVE_REAL, OldValue, dims, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'UpdateMinMaxAttribute - ModuleHDF5 - ERR02'

        if (Minimum < OldValue) then
            call h5awrite_f (attr_id, H5T_NATIVE_REAL, Minimum, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateMinMaxAttribute - ModuleHDF5 - ERR03'
        endif

        !Closes attribute 
        call h5aclose_f     (attr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UpdateMinMaxAttribute - ModuleHDF5 - ERR04'


        !Opens Maximum attribute
        call h5aopen_name_f (ID, "Maximum", attr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UpdateMinMaxAttribute - ModuleHDF5 - ERR05'

        !Reads Value
        call h5aread_f      (attr_id, H5T_NATIVE_REAL, OldValue, dims, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'UpdateMinMaxAttribute - ModuleHDF5 - ERR06'

        if (Maximum > OldValue) then
            call h5awrite_f (attr_id, H5T_NATIVE_REAL, Maximum, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateMinMaxAttribute - ModuleHDF5 - ERR07'
        endif

        !Closes attribute 
        call h5aclose_f     (attr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UpdateMinMaxAttribute - ModuleHDF5 - ERR08'

    end subroutine UpdateMinMaxAttribute

    !--------------------------------------------------------------------------
    

    subroutine HDF5ReadGlobalAttribute_Char(HDF5ID, GroupName, AttributeName, Att_Char, STAT)

        !Arguments-------------------------------------------------------------
        integer                   , intent(in)      :: HDF5ID
        character(len=*)          , intent(in)      :: GroupName
        character(len=*)          , intent(in)      :: AttributeName
        character(len=*)          , intent(out)     :: Att_Char
        integer         , optional, intent(out)     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: gr_id, attr_id
        integer                                     :: STAT_, ready_, STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call CheckGroupExistence(Me%FileID, GroupName, .false.)

            !Opens the group
            call h5gopen_f (Me%FileID, trim(GroupName), gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Char - ModuleHDF5 - ERR10'

            !Reads Real Value
            call h5aopen_name_f     (gr_id, trim(AttributeName), attr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Char - ModuleHDF5 - ERR20'
            
            call h5aread_f          (attr_id, H5T_NATIVE_CHARACTER, Att_Char, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Char - ModuleHDF5 - ERR30'
            
            call h5aclose_f         (attr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Char - ModuleHDF5 - ERR40'

            !Closes the Group
            call h5gclose_f         (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Char - ModuleHDF5 - ERR50'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5ReadGlobalAttribute_Char

    !--------------------------------------------------------------------------
    

    subroutine HDF5ReadGlobalAttribute_Int(HDF5ID, GroupName, AttributeName, Att_Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                   , intent(in)      :: HDF5ID
        character(len=*)          , intent(in)      :: GroupName
        character(len=*)          , intent(in)      :: AttributeName
        integer                   , intent(out)     :: Att_Int
        integer         , optional, intent(out)     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: gr_id, attr_id
        integer                                     :: STAT_, ready_, STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call CheckGroupExistence(Me%FileID, GroupName, .false.)

            !Opens the group
            call h5gopen_f (Me%FileID, trim(GroupName), gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Int - ModuleHDF5 - ERR10'

            !Reads Real Value
            call h5aopen_name_f     (gr_id, trim(AttributeName), attr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Int - ModuleHDF5 - ERR20'
            
            call h5aread_f          (attr_id, H5T_NATIVE_INTEGER, Att_Int, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Int - ModuleHDF5 - ERR30'
            
            call h5aclose_f         (attr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Int - ModuleHDF5 - ERR40'

            !Closes the Group
            call h5gclose_f         (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_Int - ModuleHDF5 - ERR50'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5ReadGlobalAttribute_Int

    !--------------------------------------------------------------------------


    subroutine HDF5ReadGlobalAttribute_R4(HDF5ID, GroupName, AttributeName, Att_Real, STAT)

        !Arguments-------------------------------------------------------------
        integer                   , intent(in)      :: HDF5ID
        character(len=*)          , intent(in)      :: GroupName
        character(len=*)          , intent(in)      :: AttributeName
        real(4)                   , intent(out)     :: Att_Real
        integer         , optional, intent(out)     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: gr_id, attr_id
        integer                                     :: STAT_, ready_, STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call CheckGroupExistence(Me%FileID, GroupName, .false.)

            !Opens the group
            call h5gopen_f (Me%FileID, trim(GroupName), gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R4 - ModuleHDF5 - ERR10'

            !Reads Real Value
            call h5aopen_name_f     (gr_id, trim(AttributeName), attr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R4 - ModuleHDF5 - ERR20'
            
            call h5aread_f          (attr_id, H5T_NATIVE_REAL, Att_Real, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R4 - ModuleHDF5 - ERR30'
            
            call h5aclose_f         (attr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R4 - ModuleHDF5 - ERR40'

            !Closes the Group
            call h5gclose_f         (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R4 - ModuleHDF5 - ERR50'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5ReadGlobalAttribute_R4

    !--------------------------------------------------------------------------

    subroutine HDF5ReadGlobalAttribute_R8(HDF5ID, GroupName, AttributeName, Att_Real, STAT)

        !Arguments-------------------------------------------------------------
        integer                   , intent(in)      :: HDF5ID
        character(len=*)          , intent(in)      :: GroupName
        character(len=*)          , intent(in)      :: AttributeName
        real(8)                   , intent(out)     :: Att_Real
        integer         , optional, intent(out)     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: gr_id, attr_id
        integer                                     :: STAT_, ready_, STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call CheckGroupExistence(Me%FileID, GroupName, .false.)

            !Opens the group
            call h5gopen_f (Me%FileID, trim(GroupName), gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R8 - ModuleHDF5 - ERR10'

            !Reads Real Value
            call h5aopen_name_f     (gr_id, trim(AttributeName), attr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R8 - ModuleHDF5 - ERR20'
            
            call h5aread_f          (attr_id, H5T_NATIVE_DOUBLE, Att_Real, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R8 - ModuleHDF5 - ERR30'
            
            call h5aclose_f         (attr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R8 - ModuleHDF5 - ERR40'

            !Closes the Group
            call h5gclose_f         (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGlobalAttribute_R8 - ModuleHDF5 - ERR50'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5ReadGlobalAttribute_R8

    !--------------------------------------------------------------------------    

    subroutine HDF5WriteGlobalAttribute_Char(HDF5ID, GroupName, AttributeName, Att_Char, STAT)

        !Arguments-------------------------------------------------------------
        integer                   , intent(in)      :: HDF5ID
        character(len=*)          , intent(in)      :: GroupName
        character(len=*)          , intent(in)      :: AttributeName
        character(len=*)          , intent(in)      :: Att_Char
        integer         , optional, intent(out)     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: new_type_id, gr_id
        integer                                     :: STAT_, ready_, STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call CheckGroupExistence(Me%FileID, GroupName, .false.)

            !Opens the group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR00'


            !Creates data space for attribute
            call h5screate_f (H5S_SCALAR_F, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR01'

            !Copies Type
            call h5Tcopy_f     (H5T_NATIVE_CHARACTER, new_type_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR02'

            !Sets Size
            call h5Tset_size_f (new_type_id, int(len_trim(Att_Char),8), STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR03'

            !Creates attribute
            call h5acreate_f   (gr_id, trim(AttributeName), new_type_id, space_id, attr_id1, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR04'

            !Writes attribute
            call h5awrite_f    (attr_id1, new_type_id, Att_Char, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR05'

            !Closes attribute
            call h5aclose_f (attr_id1, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR13'

            !Closes dataspaces
            call h5sclose_f (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR14'
        
            !Closes group
            call h5gclose_f(gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Char - ModuleHDF5 - ERR14'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteGlobalAttribute_Char

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteGlobalAttribute_Int(HDF5ID, GroupName, AttributeName, &
                                            Att_Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                   , intent(in)      :: HDF5ID
        character(len=*)          , intent(in)      :: GroupName
        character(len=*)          , intent(in)      :: AttributeName
        integer                   , intent(in)      :: Att_Int
        integer         , optional, intent(out)     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: gr_id
        integer                                     :: STAT_, ready_, STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call CheckGroupExistence(Me%FileID, GroupName, .false.)

            !Opens the group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Int - ModuleHDF5 - ERR00'


            !Creates data space for attribute
            call h5screate_f (H5S_SCALAR_F, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Int - ModuleHDF5 - ERR01'

             !Creates attribute
            call h5acreate_f   (gr_id, trim(AttributeName), H5T_NATIVE_INTEGER, space_id, attr_id1, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Int - ModuleHDF5 - ERR06'

            !Writes attribute
            call h5awrite_f    (attr_id1, H5T_NATIVE_INTEGER, Att_Int, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Int - ModuleHDF5 - ERR07'

           !Closes attribute
            call h5aclose_f (attr_id1, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Int - ModuleHDF5 - ERR13'

            !Closes dataspaces
            call h5sclose_f (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Int - ModuleHDF5 - ERR14'
        
            !Closes group
            call h5gclose_f(gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_Int - ModuleHDF5 - ERR14'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteGlobalAttribute_Int
    
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine HDF5WriteGlobalAttribute_R4(HDF5ID, GroupName, AttributeName, Att_Real, STAT)

        !Arguments-------------------------------------------------------------
        integer                   , intent(in)      :: HDF5ID
        character(len=*)          , intent(in)      :: GroupName
        character(len=*)          , intent(in)      :: AttributeName
        real(4)                   , intent(in)      :: Att_Real
        integer         , optional, intent(out)     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: gr_id
        integer                                     :: STAT_, ready_, STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call CheckGroupExistence(Me%FileID, GroupName, .false.)

            !Opens the group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R4 - ModuleHDF5 - ERR00'


            !Creates data space for attribute
            call h5screate_f (H5S_SCALAR_F, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R4 - ModuleHDF5 - ERR01'

            !Creates attribute
            call h5acreate_f   (gr_id, trim(AttributeName), H5T_NATIVE_REAL, space_id, attr_id1, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R4 - ModuleHDF5 - ERR08'

            !Writes attribute
            call h5awrite_f    (attr_id1, H5T_NATIVE_REAL, Att_Real, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R4 - ModuleHDF5 - ERR09'

            !Closes attribute
            call h5aclose_f (attr_id1, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R4 - ModuleHDF5 - ERR13'

            !Closes dataspaces
            call h5sclose_f (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R4 - ModuleHDF5 - ERR14'
        
            !Closes group
            call h5gclose_f(gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R4 - ModuleHDF5 - ERR14'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteGlobalAttribute_R4

    !--------------------------------------------------------------------------
    
    subroutine HDF5WriteGlobalAttribute_R8(HDF5ID, GroupName, AttributeName, Att_Real, STAT)

        !Arguments-------------------------------------------------------------
        integer                   , intent(in)      :: HDF5ID
        character(len=*)          , intent(in)      :: GroupName
        character(len=*)          , intent(in)      :: AttributeName
        real(8)                   , intent(in)      :: Att_Real
        integer         , optional, intent(out)     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: space_id
        integer(HID_T)                              :: attr_id1
        integer(HSIZE_T), dimension(7)              :: dims
        integer(HID_T)                              :: gr_id
        integer                                     :: STAT_, ready_, STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call CheckGroupExistence(Me%FileID, GroupName, .false.)

            !Opens the group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R8 - ModuleHDF5 - ERR00'


            !Creates data space for attribute
            call h5screate_f (H5S_SCALAR_F, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R8 - ModuleHDF5 - ERR01'

            !Creates attribute
            call h5acreate_f   (gr_id, trim(AttributeName), H5T_NATIVE_DOUBLE, space_id, attr_id1, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R8 - ModuleHDF5 - ERR10'

            !Writes attribute
            call h5awrite_f    (attr_id1, H5T_NATIVE_DOUBLE, Att_Real, dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R8 - ModuleHDF5 - ERR11'

            !Closes attribute
            call h5aclose_f (attr_id1, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R8 - ModuleHDF5 - ERR13'

            !Closes dataspaces
            call h5sclose_f (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R8 - ModuleHDF5 - ERR14'
        
            !Closes group
            call h5gclose_f(gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteGlobalAttribute_R8 - ModuleHDF5 - ERR14'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteGlobalAttribute_R8

    !--------------------------------------------------------------------------
    
    subroutine HDF5WriteDataR4_1D (HDF5ID, GroupName, Name, Units,                   &
                                   Array1D, Average, Radius, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        real(4), dimension(:)      , pointer            :: Array1D
        real(4), optional                               :: Average
        real(4), optional                               :: Radius
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_REAL 

            Rank    = 1
#ifndef _STACK_LIMITS_
            !Minimum = minval(Array1D(Me%Limits%ILB:Me%Limits%IUB))
            !Maximum = maxval(Array1D(Me%Limits%ILB:Me%Limits%IUB))
            Minimum = minival(Array1D,Me%Limits1D)
            Maximum = maxival(Array1D,Me%Limits1D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,          &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))
                               
            AllocateMatrix = .false.
                               
            if (.not.Associated(Me%AuxMatrixes%DataR4_1D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R4_1D(1) /= dims(1)) then
                
                deallocate(Me%AuxMatrixes%DataR4_1D)
                nullify   (Me%AuxMatrixes%DataR4_1D)

                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR4_1D(Me%Limits%ILB:Me%Limits%IUB))
                Me%AuxMatrixes%dims_R4_1D(1) = dims(1)
            endif
            
            call SetMatrixValue( Me%AuxMatrixes%DataR4_1D, Me%Limits1D, Array1D)

            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataR4_1D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteDataR4 - ModuleHDF5 - ERR08'

            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            if (present(Average) .and. present(Radius))                                 &
                 call CreateAverageRadiusAttribute (dset_id, Average, Radius)

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataR4_1D

    !--------------------------------------------------------------------------

    subroutine HDF5WriteDataR4_2D (HDF5ID, GroupName, Name, Units,                   &
                                   Array2D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        real(4), dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_REAL 

            Rank    = 2
#ifndef _STACK_LIMITS_            
            !Minimum = minval(Array2D(Me%Limits%ILB:Me%Limits%IUB,          &
            !                         Me%Limits%JLB:Me%Limits%JUB))
            !Maximum = maxval(Array2D(Me%Limits%ILB:Me%Limits%IUB,          &
            !                         Me%Limits%JLB:Me%Limits%JUB))
            Minimum = minival(Array2D,Me%Limits2D)
            Maximum = maxival(Array2D,Me%Limits2D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif

            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,      &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))

            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataR4_2D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R4_2D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_R4_2D(2) /= dims (2)) then
                
                deallocate(Me%AuxMatrixes%DataR4_2D)
                nullify   (Me%AuxMatrixes%DataR4_2D)

                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR4_2D(Me%Limits%ILB:Me%Limits%IUB,          &
                                                Me%Limits%JLB:Me%Limits%JUB))
                Me%AuxMatrixes%dims_R4_2D(1:2) = dims(1:2)
            endif
        
            call SetMatrixValue(Me%AuxMatrixes%DataR4_2D, Me%Limits2D, Array2D)

            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataR4_2D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteDataR4 - ModuleHDF5 - ERR08a'

            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataR4_2D

    !--------------------------------------------------------------------------

    subroutine HDF5WriteDataR4_3D (HDF5ID, GroupName, Name, Units,                &
                                   Array3D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        real(4), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        logical                                         :: AllocateMatrix


        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
           
            NumType = H5T_NATIVE_REAL 

            Rank    = 3
#ifndef _STACK_LIMITS_
!            Minimum = minval(Array3D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB,          &
!                                     Me%Limits%KLB:Me%Limits%KUB))
!            Maximum = maxval(Array3D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB,          &
!                                     Me%Limits%KLB:Me%Limits%KUB))
            Minimum = minival(Array3D,Me%Limits3D)
            Maximum = maxival(Array3D,Me%Limits3D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif

            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1
            dims(3) = Me%Limits%KUB - Me%Limits%KLB + 1

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,      &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))

            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataR4_3D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R4_3D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_R4_3D(2) /= dims (2) .or.                      &
                     Me%AuxMatrixes%dims_R4_3D(3) /= dims (3)) then
                
                deallocate(Me%AuxMatrixes%DataR4_3D)
                nullify   (Me%AuxMatrixes%DataR4_3D)

                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR4_3D(Me%Limits%ILB:Me%Limits%IUB,          &
                                                Me%Limits%JLB:Me%Limits%JUB,            &
                                                Me%Limits%KLB:Me%Limits%KUB))
                Me%AuxMatrixes%dims_R4_3D(1:3) = dims(1:3)
            endif
        
            call SetMatrixValue(Me%AuxMatrixes%DataR4_3D, Me%Limits3D, Array3D)
 
            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataR4_3D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'HDF5WriteDataR4 - ModuleHDF5 - ERR08b'
            endif
            
            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataR4_3D

    !--------------------------------------------------------------------------

    subroutine HDF5WriteDataR8_1D (HDF5ID, GroupName, Name, Units,                   &
                                   Array1D, Average, Radius, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        real(8), dimension(:)      , pointer            :: Array1D
        real(8), optional                               :: Average
        real(8), optional                               :: Radius
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        real(4)                                         :: Average4, Radius4
        logical                                         :: AllocateMatrix
        


        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_DOUBLE 

            Rank    = 1
#ifndef _STACK_LIMITS_
!            Minimum = minval(Array1D(Me%Limits%ILB:Me%Limits%IUB))
!            Maximum = maxval(Array1D(Me%Limits%ILB:Me%Limits%IUB))
            Minimum = minival(Array1D,Me%Limits1D)
            Maximum = maxival(Array1D,Me%Limits1D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,      &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))

            AllocateMatrix = .false.
                               
            if (.not.Associated(Me%AuxMatrixes%DataR8_1D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R8_1D(1) /= dims(1)) then

                
                deallocate(Me%AuxMatrixes%DataR8_1D)
                nullify   (Me%AuxMatrixes%DataR8_1D)
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR8_1D(Me%Limits%ILB:Me%Limits%IUB))
                Me%AuxMatrixes%dims_R8_1D(1) = dims(1)
            endif
            
            call SetMatrixValue(Me%AuxMatrixes%DataR8_1D, Me%Limits1D, Array1D)

            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataR8_1D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteDataR8 - ModuleHDF5 - ERR08'

            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            if (present(Average) .and. present(Radius)) then
                Average4 = Average
                Radius4  = Radius
                call CreateAverageRadiusAttribute (dset_id, Average4, Radius4)
            endif

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataR8_1D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteDataR8_2D (HDF5ID, GroupName, Name, Units,                   &
                                   Array2D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        real(8), dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        logical                                         :: AllocateMatrix


        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_DOUBLE 

            Rank    = 2

#ifndef _STACK_LIMITS_
!            Minimum = minval(Array2D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB))
!            Maximum = maxval(Array2D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB))
            Minimum = minival(Array2D,Me%Limits2D)
            Maximum = maxival(Array2D,Me%Limits2D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif

            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,      &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))

            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataR8_2D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R8_2D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_R8_2D(2) /= dims (2)) then

                deallocate(Me%AuxMatrixes%DataR8_2D)
                nullify   (Me%AuxMatrixes%DataR8_2D)
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR8_2D(Me%Limits%ILB:Me%Limits%IUB,          &
                                                Me%Limits%JLB:Me%Limits%JUB))
                Me%AuxMatrixes%dims_R8_2D(1:2) = dims(1:2)
            endif
        
            call SetMatrixValue(Me%AuxMatrixes%DataR8_2D, Me%Limits2D, Array2D)

            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataR8_2D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteDataR8 - ModuleHDF5 - ERR08a'
            
            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataR8_2D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteDataR8_3D (HDF5ID, GroupName, Name, Units,                   &
                                   Array3D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        real(8), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_DOUBLE 

            Rank    = 3
#ifndef _STACK_LIMITS_
!            Minimum = minval(Array3D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB,          &
!                                     Me%Limits%KLB:Me%Limits%KUB))
!            Maximum = maxval(Array3D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB,          &
!                                     Me%Limits%KLB:Me%Limits%KUB))
            Minimum = minival(Array3D,Me%Limits3D)
            Maximum = maxival(Array3D,Me%Limits3D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif

            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1
            dims(3) = Me%Limits%KUB - Me%Limits%KLB + 1

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,      &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))

            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataR8_3D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R8_3D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_R8_3D(2) /= dims (2) .or.                      &
                     Me%AuxMatrixes%dims_R8_3D(3) /= dims (3)) then

                deallocate(Me%AuxMatrixes%DataR8_3D)
                nullify   (Me%AuxMatrixes%DataR8_3D)
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR8_3D(Me%Limits%ILB:Me%Limits%IUB,      &
                                                    Me%Limits%JLB:Me%Limits%JUB,    &
                                                    Me%Limits%KLB:Me%Limits%KUB))
                Me%AuxMatrixes%dims_R8_3D(1:3) = dims(1:3)
            endif

            call SetMatrixValue(Me%AuxMatrixes%DataR8_3D, Me%Limits3D, Array3D)

            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataR8_3D,                                 &
                              dims, STAT_CALL)

            if (STAT_CALL /= SUCCESS_) then
                stop 'HDF5WriteDataR8 - ModuleHDF5 - ERR08b'
            endif

            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataR8_3D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteDataI4_1D (HDF5ID, GroupName, Name, Units,                   &
                                   Array1D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        integer(HID_T), dimension(:)      , pointer     :: Array1D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_INTEGER 

            Rank    = 1
#ifndef _STACK_LIMITS_
!            Minimum = minval(Array1D(Me%Limits%ILB:Me%Limits%IUB))
!            Maximum = maxval(Array1D(Me%Limits%ILB:Me%Limits%IUB))
            Minimum = minival(Array1D,Me%Limits1D)
            Maximum = maxival(Array1D,Me%Limits1D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif

            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1


            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,      &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))

            AllocateMatrix = .false.
                               
            if (.not.Associated(Me%AuxMatrixes%DataI4_1D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_I4_1D(1) /= dims(1)) then

                deallocate(Me%AuxMatrixes%DataI4_1D)
                nullify   (Me%AuxMatrixes%DataI4_1D)
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataI4_1D(Me%Limits%ILB:Me%Limits%IUB))
                Me%AuxMatrixes%dims_I4_1D(1) = dims(1)
            endif
            
            call SetMatrixValue(Me%AuxMatrixes%DataI4_1D, Me%Limits1D, Array1D)

            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataI4_1D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteDataI4 - ModuleHDF5 - ERR08'
            
            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataI4_1D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteDataI4_2D (HDF5ID, GroupName, Name, Units,                   &
                                   Array2D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        integer(HID_T), dimension(:, :)   , pointer     :: Array2D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                         :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_INTEGER 

            Rank    = 2
#ifndef _STACK_LIMITS_
!            Minimum = minval(Array2D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB))
!            Maximum = maxval(Array2D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB))
            Minimum = minival(Array2D,Me%Limits2D)
            Maximum = maxival(Array2D,Me%Limits2D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,      &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))
                               
            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataI4_2D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_I4_2D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_I4_2D(2) /= dims (2)) then

                deallocate(Me%AuxMatrixes%DataI4_2D)
                nullify   (Me%AuxMatrixes%DataI4_2D)
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataI4_2D(Me%Limits%ILB:Me%Limits%IUB,  &
                                                    Me%Limits%JLB:Me%Limits%JUB))
                Me%AuxMatrixes%dims_I4_2D(1:2) = dims(1:2)
            endif
        
            call SetMatrixValue(Me%AuxMatrixes%DataI4_2D, Me%Limits2D, Array2D)
 

            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataI4_2D,                                 &
                              dims, STAT_CALL)                               
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteDataI4 - ModuleHDF5 - ERR08a'
            
            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataI4_2D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteDataI4_3D (HDF5ID, GroupName, Name, Units,                   &
                                   Array3D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        character(len=*)                                :: Units
        integer(HID_T), dimension(:, :, :), pointer     :: Array3D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                         :: STAT_CALL
        integer(HID_T)                                  :: dset_id, prp_id, gr_id
        character(StringLength)                         :: AuxChar
        real(4)                                         :: Minimum, Maximum
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_INTEGER 

            Rank    = 3
#ifndef _STACK_LIMITS_
!            Minimum = minval(Array3D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB,          &
!                                     Me%Limits%KLB:Me%Limits%KUB))
!            Maximum = maxval(Array3D(Me%Limits%ILB:Me%Limits%IUB,          &
!                                     Me%Limits%JLB:Me%Limits%JUB,          &
!                                     Me%Limits%KLB:Me%Limits%KUB))
            Minimum = minival(Array3D,Me%Limits3D)
            Maximum = maxival(Array3D,Me%Limits3D)
#else
            Minimum = FillValueReal
            Maximum = FillValueReal
#endif
           
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1
            dims(3) = Me%Limits%KUB - Me%Limits%KLB + 1

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif
            
            !Opens Group, Creates Dset, etc
            call PrepareWrite (Me%FileID, Rank, dims, space_id, prp_id, gr_id,      &
                               dset_id, NumType, GroupName, trim(adjustl(AuxChar)))

            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataI4_3D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_I4_3D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_I4_3D(2) /= dims (2) .or.                      &
                     Me%AuxMatrixes%dims_I4_3D(3) /= dims (3)) then

                deallocate(Me%AuxMatrixes%DataI4_3D)
                nullify   (Me%AuxMatrixes%DataI4_3D)
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataI4_3D(Me%Limits%ILB:Me%Limits%IUB,      &   
                                                    Me%Limits%JLB:Me%Limits%JUB,    & 
                                                    Me%Limits%KLB:Me%Limits%KUB))
                Me%AuxMatrixes%dims_I4_3D(1:3) = dims(1:3)
            endif
        
            call SetMatrixValue(Me%AuxMatrixes%DataI4_3D,Me%Limits3D,Array3D)
 

            !Writes the data to the file
            call h5dwrite_f  (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataI4_3D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'HDF5WriteDataI4 - ModuleHDF5 - ERR08b'
            endif
            
            !Creates attributes
            call CreateMinMaxAttribute (dset_id, Units, Minimum, Maximum)

            !Updates attributes of Group
            call UpdateMinMaxAttribute (gr_id, Minimum, Maximum)

            !Closes Group, Releases Dset, etc
            call FinishWrite (space_id, prp_id, gr_id, dset_id)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5WriteDataI4_3D
    
    !--------------------------------------------------------------------------

    subroutine PrepareWrite (FileID, Rank, dims, space_id, prp_id, gr_id, dset_id,       &
                             NumType, GroupName, ItemName)

        !Arguments-------------------------------------------------------------
        integer (HID_T)                             :: FileID, Rank, space_id
        integer (HID_T)                             :: prp_id, gr_id, dset_id, NumType
        integer (HSIZE_T), dimension(7)             :: dims
        character(len=*)                            :: GroupName
        character(len=*)                            :: ItemName

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: STAT_CALL

        !Creates a simple dataspace
        call h5screate_simple_f(Rank, dims, space_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) GroupName
            write(*,*) ItemName
            stop 'PrepareWrite - ModuleHDF5 - ERR01'
        endif

        !Creates a property list
        call h5pcreate_f (H5P_DATASET_CREATE_F, prp_id, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleHDF5 - ERR02'

        !Sets chunked
        call h5pset_chunk_f(prp_id, Rank, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleHDF5 - ERR03'

        !Sets the compression
        call h5pset_deflate_f(prp_id, 6, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleHDF5 - ERR04'

        !Verifies if group exists, if not create it
        call CheckGroupExistence (FileID, GroupName)

        !Opens the Group
        call h5gopen_f (FileID, GroupName, gr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleHDF5 - ERR05'

        !Creates the dataset with default properties
        call h5dcreate_f(gr_id, ItemName, NumType, space_id, dset_id,  STAT_CALL, prp_id) 
        if (STAT_CALL /= SUCCESS_)then
            write(*,*)trim(ItemName)
            write(*,*)trim(Me%FileName)
            stop 'PrepareWrite - ModuleHDF5 - ERR06'
        end if

    end subroutine PrepareWrite

    !--------------------------------------------------------------------------

    recursive subroutine CheckGroupExistence (FileID, GroupName, CreateMinMaxAttributes)

        !Arguments-------------------------------------------------------------
        integer (HID_T)                             :: FileID
        character(len=*)                            :: GroupName
        logical, optional, intent(in)               :: CreateMinMaxAttributes

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: nmembers    
        integer(HID_T)                              :: STAT_CALL
        character(StringLength)                     :: ParentGroupName
        integer  (HID_T)                            :: gr_id
        logical                                     :: lCreateMinMaxAttributes

        if(present(CreateMinMaxAttributes))then
            lCreateMinMaxAttributes = CreateMinMaxAttributes
        else
            lCreateMinMaxAttributes = .true.
        endif

        !Turns Error printing of
        call h5eset_auto_f  (0, STAT_CALL)

        call h5gn_members_f (FileID, trim(GroupName), nmembers, STAT_CALL)

        if (STAT_CALL == -1) then
            ParentGroupName = GroupName(1:index(GroupName, "/", .true.)-1)

            if (len_trim(ParentGroupName) > 0) then
                call h5gn_members_f (FileID, trim(ParentGroupName), nmembers, STAT_CALL)
                if (STAT_CALL == -1) then
                    call CheckGroupExistence (FileID, ParentGroupName)
                endif
            endif

            !Creates a new group with a given name
            call h5gcreate_f(FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckGroupExistence - ModuleHDF5 - ERR01'

            !write(*,*)'Created Group: ', trim(GroupName)
            if(lCreateMinMaxAttributes)then
                call CreateMinMaxAttribute (gr_id)   
            end if

            !Closes the group
            call h5gclose_f(gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckGroupExistence - ModuleHDF5 - ERR02'

        endif
      
        !Turns Error printing on
        call h5eset_auto_f  (1, STAT_CALL)


    end subroutine CheckGroupExistence

    !--------------------------------------------------------------------------
 
    subroutine FinishWrite (space_id, prp_id, gr_id, dset_id)

        !Arguments-------------------------------------------------------------
        integer (HID_T)                             :: space_id, prp_id, gr_id
        integer (HID_T)                             :: dset_id

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: STAT_CALL

        !Closes property list
        call h5pclose_f  (prp_id, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'FinishWrite - ModuleHDF5 - ERR03'

        !End access to the dataset
        call h5dclose_f  (dset_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FinishWrite - ModuleHDF5 - ERR02'

        !End access to the dataspace
        call h5sclose_f  (space_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FinishWrite - ModuleHDF5 - ERR04'

        !Closes group
        call h5gclose_f( gr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FinishWrite - ModuleHDF5 - ERR01'


    end subroutine FinishWrite

    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataR4_1D(HDF5ID, GroupName, Name,                          &
                                 Array1D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:)      , pointer            :: Array1D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix
        

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_REAL

            Rank    = 1
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1

            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_1D - ModuleHDF5 - ERR05'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_1D - ModuleHDF5 - ERR05'

            AllocateMatrix = .false.
                               
            if (.not.Associated(Me%AuxMatrixes%DataR4_1D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R4_1D(1) /= dims(1)) then
                
                deallocate(Me%AuxMatrixes%DataR4_1D)
                nullify   (Me%AuxMatrixes%DataR4_1D)                
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR4_1D(Me%Limits%ILB:Me%Limits%IUB))
                Me%AuxMatrixes%dims_R4_1D(1) = dims(1)
            endif
            
            !Read the data to the file
            call h5dread_f  (dset_id, NumType,                                         &
                             Me%AuxMatrixes%DataR4_1D,                                 &
                             dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_1D - ModuleHDF5 - ERR08'
            
            call SetMatrixValue(Array1D, Me%Limits1D, Me%AuxMatrixes%DataR4_1D)           

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_1D - ModuleHDF5 - ERR09'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_1D - ModuleHDF5 - ERR07'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataR4_1D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataR4_2D(HDF5ID, GroupName, Name,                          &
                                 Array2D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_REAL

            Rank    = 2
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1

            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_2D - ModuleHDF5 - ERR05'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_2D - ModuleHDF5 - ERR05'

                               
            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataR4_2D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R4_2D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_R4_2D(2) /= dims (2)) then
                
                deallocate(Me%AuxMatrixes%DataR4_2D)
                nullify   (Me%AuxMatrixes%DataR4_2D)                
                                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR4_2D(Me%Limits%ILB:Me%Limits%IUB,      &
                                                    Me%Limits%JLB:Me%Limits%JUB))
                Me%AuxMatrixes%dims_R4_2D(1:2) = dims(1:2)
            endif
        
 

            !Read the data to the file
            call h5dread_f  (dset_id, NumType,                                         &
                             Me%AuxMatrixes%DataR4_2D,                                 &
                             dims, STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_2D - ModuleHDF5 - ERR08a'
            
            call SetMatrixValue(Array2D, Me%Limits2D, Me%AuxMatrixes%DataR4_2D)

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_2D - ModuleHDF5 - ERR09'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_2D - ModuleHDF5 - ERR07'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataR4_2D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataR4_3D(HDF5ID, GroupName, Name,                          &
                                 Array3D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_REAL
            Rank    = 3
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1
            dims(3) = Me%Limits%KUB - Me%Limits%KLB + 1


            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'Group name ', GroupName, ' not found.'
                stop 'HDF5ReadDataR4_3D - ModuleHDF5 - ERR05'
            endif

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'Dataset ', AuxChar, ' not found.'
                stop 'HDF5ReadDataR4_3D - ModuleHDF5 - ERR05b'
            endif
            
            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataR4_3D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R4_3D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_R4_3D(2) /= dims (2) .or.                      &
                     Me%AuxMatrixes%dims_R4_3D(3) /= dims (3)) then

                
                deallocate(Me%AuxMatrixes%DataR4_3D)
                nullify   (Me%AuxMatrixes%DataR4_3D)                
                                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR4_3D(Me%Limits%ILB:Me%Limits%IUB,      &
                                                Me%Limits%JLB:Me%Limits%JUB,        &
                                                Me%Limits%KLB:Me%Limits%KUB))
                Me%AuxMatrixes%dims_R4_3D(1:3) = dims(1:3)
            endif
        
 

            !Read the data to the file
            call h5dread_f   (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataR4_3D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'HDF5ReadDataR4_3D - ModuleHDF5 - ERR08b'
            endif

            call SetMatrixValue(Array3D, Me%Limits3D, Me%AuxMatrixes%DataR4_3D)

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_3D - ModuleHDF5 - ERR09'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4_3D - ModuleHDF5 - ERR07'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataR4_3D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataR8_1D(HDF5ID, GroupName, Name,                          &
                                 Array1D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:)      , pointer            :: Array1D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_DOUBLE
            Rank    = 1
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1

            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_1D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_1D - ModuleHDF5 - ERR20'
            
            AllocateMatrix = .false.
                               
            if (.not.Associated(Me%AuxMatrixes%DataR8_1D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R8_1D(1) /= dims(1)) then
                
                deallocate(Me%AuxMatrixes%DataR8_1D)
                nullify   (Me%AuxMatrixes%DataR8_1D)                
                                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                
                nullify(Me%AuxMatrixes%DataR8_1D)

                allocate(Me%AuxMatrixes%DataR8_1D(Me%Limits%ILB:Me%Limits%IUB))
                Me%AuxMatrixes%dims_R8_1D(1) = dims(1)
            endif
            

            !Read the data to the file
            call h5dread_f  (dset_id, NumType,                                          &
                              Me%AuxMatrixes%DataR8_1D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_1D - ModuleHDF5 - ERR30'

            call SetMatrixValue(Array1D, Me%Limits1D, Me%AuxMatrixes%DataR8_1D)

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_1D - ModuleHDF5 - ERR40'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_1D - ModuleHDF5 - ERR50'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataR8_1D
    
    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataR8_2D(HDF5ID, GroupName, Name,                          &
                                 Array2D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_DOUBLE

            Rank    = 2
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1



            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_2D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_2D - ModuleHDF5 - ERR20'

                               
            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataR8_2D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R8_2D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_R8_2D(2) /= dims (2)) then
                
                deallocate(Me%AuxMatrixes%DataR8_2D)
                nullify   (Me%AuxMatrixes%DataR8_2D) 
                                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR8_2D(Me%Limits%ILB:Me%Limits%IUB,  &
                                                    Me%Limits%JLB:Me%Limits%JUB))
                Me%AuxMatrixes%dims_R8_2D(1:2) = dims(1:2)
            endif
        

            !Read the data to the file
            call h5dread_f  (dset_id, NumType,                                         &
                             Me%AuxMatrixes%DataR8_2D,                                 &
                             dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_2D - ModuleHDF5 - ERR30'

            call SetMatrixValue(Array2D, Me%Limits2D, Me%AuxMatrixes%DataR8_2D)
                             
            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_2D - ModuleHDF5 - ERR40'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_2D - ModuleHDF5 - ERR50'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataR8_2D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataR8_3D(HDF5ID, GroupName, Name,                          &
                                 Array3D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_DOUBLE

            Rank    = 3
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1
            dims(3) = Me%Limits%KUB - Me%Limits%KLB + 1

            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_3D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_3D - ModuleHDF5 - ERR20'

            
            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataR8_3D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_R8_3D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_R8_3D(2) /= dims (2) .or.                      &
                     Me%AuxMatrixes%dims_R8_3D(3) /= dims (3)) then
                
                deallocate(Me%AuxMatrixes%DataR8_3D)
                nullify   (Me%AuxMatrixes%DataR8_3D) 
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataR8_3D(Me%Limits%ILB:Me%Limits%IUB,          &
                                                Me%Limits%JLB:Me%Limits%JUB,            &
                                                Me%Limits%KLB:Me%Limits%KUB))
                Me%AuxMatrixes%dims_R8_3D(1:3) = dims(1:3)
            endif
        
            !Read the data to the file
            call h5dread_f   (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataR8_3D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_3D - ModuleHDF5 - ERR30'

            call SetMatrixValue(Array3D, Me%Limits3D, Me%AuxMatrixes%DataR8_3D)
                             
            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_3D - ModuleHDF5 - ERR40'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR8_3D - ModuleHDF5 - ERR50'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataR8_3D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataI4_1D(HDF5ID, GroupName, Name,                          &
                                 Array1D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer(4), dimension(:)      , pointer         :: Array1D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_INTEGER
            Rank    = 1
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1

            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_1D - ModuleHDF5 - ERR05a'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_1D - ModuleHDF5 - ERR05b'

            
            AllocateMatrix = .false.
                               
            if (.not.Associated(Me%AuxMatrixes%DataI4_1D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_I4_1D(1) /= dims(1)) then
                
                deallocate(Me%AuxMatrixes%DataI4_1D)
                nullify   (Me%AuxMatrixes%DataI4_1D) 
                                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataI4_1D(Me%Limits%ILB:Me%Limits%IUB))
                Me%AuxMatrixes%dims_I4_1D(1) = dims(1)
            endif
            

            !Read the data to the file
            call h5dread_f  (dset_id, NumType,                                          &
                             Me%AuxMatrixes%DataI4_1D,                                  &
                             dims, STAT_CALL)
                              
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_1D - ModuleHDF5 - ERR08'
            
            call SetMatrixValue(Array1D, Me%Limits1D, Me%AuxMatrixes%DataI4_1D)

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_1D - ModuleHDF5 - ERR09'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_1D - ModuleHDF5 - ERR07'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataI4_1D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataI4_2D(HDF5ID, GroupName, Name,                          &
                                 Array2D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer(4), dimension(:, :)   , pointer         :: Array2D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                         :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_INTEGER

            Rank    = 2
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1

            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_2D - ModuleHDF5 - ERR05c'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_2D - ModuleHDF5 - ERR05d'

                               
            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataI4_2D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_I4_2D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_I4_2D(2) /= dims (2)) then

                deallocate(Me%AuxMatrixes%DataI4_2D)
                nullify   (Me%AuxMatrixes%DataI4_2D)                 

                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataI4_2D(Me%Limits%ILB:Me%Limits%IUB,      &
                                                Me%Limits%JLB:Me%Limits%JUB))
                Me%AuxMatrixes%dims_I4_2D(1:2) = dims(1:2)
            endif
        
 

            !Read the data to the file
            call h5dread_f  (dset_id, NumType,                                          &
                             Me%AuxMatrixes%DataI4_2D,                                  &
                             dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_2D - ModuleHDF5 - ERR08a'

            call SetMatrixValue(Array2D,Me%Limits2D,Me%AuxMatrixes%DataI4_2D)
                             

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_2D - ModuleHDF5 - ERR09'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_2D - ModuleHDF5 - ERR07'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataI4_2D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadDataI4_3D(HDF5ID, GroupName, Name,                          &
                                 Array3D, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer(4), dimension(:, :, :), pointer         :: Array3D
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims
        integer(HID_T)                                  :: Rank
        integer(HID_T)                                  :: NumType
        integer(HID_T)                                         :: STAT_CALL
        integer(HID_T)                                  :: dset_id, gr_id
        character(StringLength)                         :: AuxChar
        logical                                         :: AllocateMatrix

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            NumType = H5T_NATIVE_INTEGER
            Rank    = 3
            dims(1) = Me%Limits%IUB - Me%Limits%ILB + 1
            dims(2) = Me%Limits%JUB - Me%Limits%JLB + 1
            dims(3) = Me%Limits%KUB - Me%Limits%KLB + 1

            !Opens the Group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_3D - ModuleHDF5 - ERR05e'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_3D - ModuleHDF5 - ERR05f'

            
            AllocateMatrix = .false.
                                   
            if (.not.Associated(Me%AuxMatrixes%DataI4_3D)) then
                
                AllocateMatrix = .true.
                
            else if (Me%AuxMatrixes%dims_I4_3D(1) /= dims (1) .or.                      &
                     Me%AuxMatrixes%dims_I4_3D(2) /= dims (2) .or.                      &
                     Me%AuxMatrixes%dims_I4_3D(3) /= dims (3)) then

                deallocate(Me%AuxMatrixes%DataI4_3D)
                nullify   (Me%AuxMatrixes%DataI4_3D)                 
                
                AllocateMatrix = .true.
                
            endif
            
            if (AllocateMatrix) then
                allocate(Me%AuxMatrixes%DataI4_3D(Me%Limits%ILB:Me%Limits%IUB,          &
                                                Me%Limits%JLB:Me%Limits%JUB,            &
                                                Me%Limits%KLB:Me%Limits%KUB))
                Me%AuxMatrixes%dims_I4_3D(1:3) = dims(1:3)
            endif
        
            !Read the data to the file
            call h5dread_f   (dset_id, NumType,                                         &
                              Me%AuxMatrixes%DataI4_3D,                                 &
                              dims, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'HDF5ReadDataI4_3D - ModuleHDF5 - ERR08b'
            endif
            
            call SetMatrixValue(Array3D,Me%Limits3D,Me%AuxMatrixes%DataI4_3D)
                                         
            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_3D - ModuleHDF5 - ERR09'

            !Closes group
            call h5gclose_f( gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataI4_3D - ModuleHDF5 - ERR07'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadDataI4_3D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadHyperSlab (HDF5ID, GroupName, ItemName,                 &
                                    lower_bound, upper_bound,                  &
                                    Array1D,   Array2D,   Array3D,             &
                                    TryToRead, OutputNumber, STAT)
#ifdef _GUI_
        !DEC$ IF DEFINED(_X86_)
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : '_HDF5ReadHyperSlab@40'  :: HDF5ReadHyperSlab
        !DEC$ ELSE
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : 'HDF5ReadHyperSlab'      :: HDF5ReadHyperSlab
        !DEC$ ENDIF
#endif
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: ItemName
        integer, dimension(3)                           :: lower_bound, upper_bound
        real(4), dimension(:)      , pointer, optional  :: Array1D
        real(4), dimension(:, :)   , pointer, optional  :: Array2D
        real(4), dimension(:, :, :), pointer, optional  :: Array3D
        logical,                  intent(IN), optional  :: TryToRead
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T), dimension(:),       pointer     :: iArray1D
        integer(HID_T), dimension(:, :),    pointer     :: iArray2D
        integer(HID_T), dimension(:, :, :), pointer     :: iArray3D
        logical                                         :: TryToRead_
        character(StringLength)                         :: AuxChar

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if (present(TryToRead)) then
                TryToRead_ = TryToRead
            else
                TryToRead_ = .false.
            endif

            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR01'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (ItemName, OutputNumber, AuxChar)
            else
                AuxChar = ItemName
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)

            if (TryToRead_ .and. STAT_CALL /= SUCCESS_) return

            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR02'

            !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR03'

            call h5tget_class_f(datatype_id, class_id,  STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR04a'

            if      (class_id == H5T_FLOAT_F) then
                NumType = H5T_NATIVE_REAL
            elseif  (class_id == H5T_INTEGER_F) then
                NumType = H5T_NATIVE_INTEGER
                if     (present(Array1D)) then
                    allocate (iArray1D(lower_bound(1):upper_bound(1)))
                elseif (present(Array2D)) then
                    allocate (iArray2D(lower_bound(1):upper_bound(1),                    &
                                       lower_bound(2):upper_bound(2)))
                elseif (present(Array3D)) then
                    allocate (iArray3D(lower_bound(1):upper_bound(1),                    &
                                       lower_bound(2):upper_bound(2),                    &
                                       lower_bound(3):upper_bound(3)))
                endif
            else
                stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR04b'
            endif

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR05'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR06'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            !if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR07'


            allocate (offset_in (rank))
            allocate (count_in  (rank))

            offset_in(1:rank) = lower_bound(1:rank) - 1
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR08'


            !Defines the memory dataspace
            if      (present(Array1D)) then
                rank_out    = 1
                dims_mem(1) = size(Array1D)
            elseif  (present(Array2D)) then
                rank_out = 2
                dims_mem(1) = size(Array2D, dim = 1)
                dims_mem(2) = size(Array2D, dim = 2)
            elseif  (present(Array3D)) then
                rank_out = 3
                dims_mem(1) = size(Array3D, dim = 1)
                dims_mem(2) = size(Array3D, dim = 2)
                dims_mem(3) = size(Array3D, dim = 3)
            endif
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR09'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR10'

            if (NumType == H5T_NATIVE_INTEGER) then
                if      (present(Array1D)) then
                    call h5dread_f (dset_id, NumType, iArray1D(lower_bound(1):upper_bound(1)),   &
                                    dims_mem, STAT_CALL, memspace_id, space_id)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR10a'
                    Array1D =  iArray1D
                    deallocate(iArray1D)
                elseif  (present(Array2D)) then
                    call h5dread_f (dset_id, NumType, iArray2D(lower_bound(1):upper_bound(1),    &
                                                               lower_bound(2):upper_bound(2)),   &
                                    dims_mem, STAT_CALL, memspace_id, space_id)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR10b'
                    Array2D =  iArray2D
                    deallocate(iArray2D)
                elseif  (present(Array3D)) then
                    call h5dread_f (dset_id, NumType, iArray3D(lower_bound(1):upper_bound(1),    &
                                                               lower_bound(2):upper_bound(2),    &
                                                               lower_bound(3):upper_bound(3)),   &
                                    dims_mem, STAT_CALL, memspace_id, space_id)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR10c'
                    Array3D =  iArray3D
                    deallocate(iArray3D)
                endif
            else
                if      (present(Array1D)) then
                    call h5dread_f (dset_id, NumType, Array1D(lower_bound(1):upper_bound(1)),   &
                                    dims_mem, STAT_CALL, memspace_id, space_id)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR10a'
                elseif  (present(Array2D)) then
                    call h5dread_f (dset_id, NumType, Array2D(lower_bound(1):upper_bound(1),    &
                                                              lower_bound(2):upper_bound(2)),   &
                                    dims_mem, STAT_CALL, memspace_id, space_id)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR10b'
                elseif  (present(Array3D)) then
                    call h5dread_f (dset_id, NumType, Array3D(lower_bound(1):upper_bound(1),    &
                                                              lower_bound(2):upper_bound(2),    &
                                                              lower_bound(3):upper_bound(3)),   &
                                    dims_mem, STAT_CALL, memspace_id, space_id)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR10c'
                endif
            endif


            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR05'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR05'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR05'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadHyperSlab - ModuleHDF5 - ERR06'


            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine HDF5ReadHyperSlab

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowR4_1D   (HDF5ID, GroupName, Name,                          &
                                    Array1D, OffSet1,                                   &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:)      , pointer            :: Array1D
        integer, optional                               :: OffSet1
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_REAL
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR80'
            endif
            

            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR80'


            !Defines the memory dataspace
            rank_out    = 1
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR90'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR100'

            call h5dread_f (dset_id, NumType, Array1D(lower_bound(1):upper_bound(1)),&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR110'

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR140'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR150'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR160'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_1D - ModuleHDF5 - ERR170'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowR4_1D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowR4_2D (HDF5ID, GroupName, Name,                            &
                                    Array2D, OffSet1, OffSet2,                          &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_REAL
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif


            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR80'
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR90'


            !Defines the memory dataspace
            rank_out    = 2
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR100'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR110'

            call h5dread_f (dset_id, NumType, Array2D(lower_bound(1):upper_bound(1), &
                                                      lower_bound(2):upper_bound(2)),&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR120'

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR140'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR150'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR160'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_2D - ModuleHDF5 - ERR170'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowR4_2D

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowR4_3D   (HDF5ID, GroupName, Name,                          &
                                      Array3D,                                          &
                                      OffSet1, OffSet2, OffSet3,                        &
                                      OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OffSet3                
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_REAL
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif


            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                if (present(OffSet3)) then
                    offset_in(3) = OffSet3
                else                
                    offset_in(3) = lower_bound(3) - 1    
                endif            
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR80'


            !Defines the memory dataspace
            rank_out    = 3
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            dims_mem(3) = upper_bound(3)-lower_bound(3) + 1
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR90'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR100'

            call h5dread_f (dset_id, NumType, Array3D(lower_bound(1):upper_bound(1), &
                                                      lower_bound(2):upper_bound(2), &
                                                      lower_bound(3):upper_bound(3)),&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR130'

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR140'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR150'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR160'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR4_3D - ModuleHDF5 - ERR170'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowR4_3D

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowR8_1D   (HDF5ID, GroupName, Name,                          &
                                    Array1D, OffSet1,                                   &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:)      , pointer            :: Array1D
        integer, optional                               :: OffSet1
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_DOUBLE
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR80'
            endif
            

            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR80'


            !Defines the memory dataspace
            rank_out    = 1
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR90'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR100'

            call h5dread_f (dset_id, NumType, Array1D(lower_bound(1):upper_bound(1)),&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR110'

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR140'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR150'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR160'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_1D - ModuleHDF5 - ERR170'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowR8_1D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowR8_2D (HDF5ID, GroupName, Name,                            &
                                    Array2D, OffSet1, OffSet2,                          &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_DOUBLE
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif


            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR80'
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR90'


            !Defines the memory dataspace
            rank_out    = 2
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR100'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR110'

            call h5dread_f (dset_id, NumType, Array2D(lower_bound(1):upper_bound(1), &
                                                      lower_bound(2):upper_bound(2)),&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR120'

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR140'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR150'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR160'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_2D - ModuleHDF5 - ERR170'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowR8_2D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowR8_3D   (HDF5ID, GroupName, Name,                          &
                                      Array3D,                                          &
                                      OffSet1, OffSet2, OffSet3,                        &
                                      OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OffSet3                
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_DOUBLE
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif


            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                if (present(OffSet3)) then
                    offset_in(3) = OffSet3
                else                
                    offset_in(3) = lower_bound(3) - 1    
                endif            
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR80'


            !Defines the memory dataspace
            rank_out    = 3
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            dims_mem(3) = upper_bound(3)-lower_bound(3) + 1
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR90'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR100'

            call h5dread_f (dset_id, NumType, Array3D(lower_bound(1):upper_bound(1), &
                                                      lower_bound(2):upper_bound(2), &
                                                      lower_bound(3):upper_bound(3)),&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR130'

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR140'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR150'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR160'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowR8_3D - ModuleHDF5 - ERR170'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowR8_3D

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowI4_1D (HDF5ID, GroupName, Name,                            &
                                    Array1D, OffSet1,                                   &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer, dimension(:)      , pointer            :: Array1D
        integer, optional                               :: OffSet1
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_INTEGER 
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR80'
            endif
            

            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR80'


            !Defines the memory dataspace
            rank_out    = 1
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR90'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR100'

            call h5dread_f (dset_id, NumType, Array1D(lower_bound(1):upper_bound(1)),&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR110'

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR140'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR150'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR160'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_1D - ModuleHDF5 - ERR170'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowI4_1D

    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowI4_2D (HDF5ID, GroupName, Name,                            &
                                    Array2D, OffSet1, OffSet2,                          &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer, dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        integer, dimension(:, :, :)   , pointer         :: Array3D
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            NumType = H5T_NATIVE_INTEGER 
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) AuxChar
                stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR20'
            endif
                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR60'

            if (rank >3) then
                stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR70'
            endif

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR80'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank ==3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif

            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >= 2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (rank == 3) then
                offset_in(3) = 0
            endif

            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR90'


            !Defines the memory dataspace
            rank_out    = 2
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            
            if (rank == 3) then
                rank_out    = 3
                dims_mem(3) = upper_bound(3)-lower_bound(3) + 1
            endif
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR100'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR120'
            
            if (rank ==3) then

                allocate(Array3D(lower_bound(1):upper_bound(1), lower_bound(2):upper_bound(2), &
                                 lower_bound(3):upper_bound(3)))

                call h5dread_f (dset_id, NumType, Array3D(lower_bound(1):upper_bound(1), &
                                                          lower_bound(2):upper_bound(2), &
                                                          lower_bound(3):upper_bound(3)), &
                                dims_mem, STAT_CALL, memspace_id, space_id)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR130'
                
                Array2D(lower_bound(1):upper_bound(1), lower_bound(2):upper_bound(2)) = & 
                Array3D(lower_bound(1):upper_bound(1), lower_bound(2):upper_bound(2), upper_bound(3))
                

                deallocate(Array3D)
                

            else                

                call h5dread_f (dset_id, NumType, Array2D(lower_bound(1):upper_bound(1), &
                                                          lower_bound(2):upper_bound(2)),&
                                dims_mem, STAT_CALL, memspace_id, space_id)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR140'

            endif

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR150'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR160'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR170'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_2D - ModuleHDF5 - ERR180'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowI4_2D

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine HDF5ReadWindowI4_3D   (HDF5ID, GroupName, Name,                          &
                                      Array3D,                                          &
                                      OffSet1, OffSet2, OffSet3,                        &
                                      OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer, dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OffSet3                
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank, rank_out
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_out
        integer(HSIZE_T ), dimension(:), allocatable    :: count_out
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_INTEGER 
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR10'


            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif


            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif


            allocate (offset_in (rank))
            allocate (count_in  (rank))

            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                if (present(OffSet3)) then
                    offset_in(3) = OffSet3
                else                
                    offset_in(3) = lower_bound(3) - 1    
                endif            
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR80'


            !Defines the memory dataspace
            rank_out    = 3
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            dims_mem(3) = upper_bound(3)-lower_bound(3) + 1
            
            call h5screate_simple_f (rank_out, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR90'


            !Define the memory hyperslab
            allocate (offset_out(rank_out))
            allocate (count_out (rank_out))
            offset_out = 0
            count_out  = count_in
            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, offset_out, count_out, &
                                        STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR100'

            call h5dread_f (dset_id, NumType, Array3D(lower_bound(1):upper_bound(1), &
                                                      lower_bound(2):upper_bound(2), &
                                                      lower_bound(3):upper_bound(3)),&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR130'

            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            deallocate (offset_out)
            deallocate (count_out )

            !Closes data space
            call h5sclose_f  (memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR140'

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR150'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR160'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR170'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5ReadWindowI4_3D
    
        !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowR4_3D   (HDF5ID, GroupName, Name,                         &
                                      Array3D,                                          &
                                      OffSet1, OffSet2, OffSet3,                        &
                                      OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OffSet3                
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_REAL
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR20'

            !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif

            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                if (present(OffSet3)) then
                    offset_in(3) = OffSet3
                else                
                    offset_in(3) = lower_bound(3) - 1    
                endif            
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR80'

            !Defines the memory dataspace
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            dims_mem(3) = upper_bound(3)-lower_bound(3) + 1
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR90'

            call h5dwrite_f (dset_id, NumType, Array3D,&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR100'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            
            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR110'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR120'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_3D - ModuleHDF5 - ERR130'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowR4_3D
        
        !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowR4_2D (HDF5ID, GroupName, Name,                    &
                                    Array2D, OffSet1, OffSet2,                   &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_REAL
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR20'

            !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif

            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR80'
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR90'

            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR100'
            
            call h5dwrite_f (dset_id, NumType, Array2D,&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR110'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR120'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR130'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_2D - ModuleHDF5 - ERR140'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowR4_2D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowR4_1D  (HDF5ID, GroupName, Name,                          &
                                    Array1D, OffSet1,                                   &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(4), dimension(:)      , pointer            :: Array1D
        integer, optional                               :: OffSet1
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_REAL
            
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR70'

            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR80'
            endif
            
            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR80'

            !Defines the memory dataspace
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR90'

            call h5dwrite_f (dset_id, NumType, Array1D,&
                             dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR100'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            
            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR110'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR120'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR4_1D - ModuleHDF5 - ERR130'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowR4_1D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowR8_3D  (HDF5ID, GroupName, Name,                          &
                                      Array3D,                                          &
                                      OffSet1, OffSet2, OffSet3,                        &
                                      OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OffSet3                
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_DOUBLE
            
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR20'

            !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif

            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                if (present(OffSet3)) then
                    offset_in(3) = OffSet3
                else                
                    offset_in(3) = lower_bound(3) - 1    
                endif            
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR80'

            !Defines the memory dataspace
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            dims_mem(3) = upper_bound(3)-lower_bound(3) + 1
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR90'

            call h5dwrite_f (dset_id, NumType, Array3D,&
                             dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR100'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            
            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR110'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR120'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_3D - ModuleHDF5 - ERR130'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowR8_3D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowR8_2D (HDF5ID, GroupName, Name,                            &
                                    Array2D, OffSet1, OffSet2,                          &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_DOUBLE
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR20'

            !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif

            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR80'
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR90'

            !Defines the memory dataspace
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR100'

            call h5dwrite_f (dset_id, NumType, Array2D,&
                             dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR110'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR120'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR130'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR140'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowR8_2D
    
    !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowR8_1D  (HDF5ID, GroupName, Name,                          &
                                    Array1D, OffSet1,                                   &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        real(8), dimension(:)      , pointer            :: Array1D
        integer, optional                               :: OffSet1
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_DOUBLE
            
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR20'

            !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR80'
            endif
            
            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR80'

            !Defines the memory dataspace
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR90'

            call h5dwrite_f (dset_id, NumType, Array1D,&
                             dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_2D - ModuleHDF5 - ERR100'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR110'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR120'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowR8_1D - ModuleHDF5 - ERR130'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowR8_1D    
        
    !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowI4_3D  (HDF5ID, GroupName, Name,                          &
                                      Array3D,                                          &
                                      OffSet1, OffSet2, OffSet3,                        &
                                      OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer, dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OffSet3                
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_INTEGER 
        
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR20'

                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank >=3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif

            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >=2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (Rank >=3) then
                if (present(OffSet3)) then
                    offset_in(3) = OffSet3
                else                
                    offset_in(3) = lower_bound(3) - 1    
                endif            
            endif
            
            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR80'

            !Defines the memory dataspace
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            dims_mem(3) = upper_bound(3)-lower_bound(3) + 1
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadWindowI4_3D - ModuleHDF5 - ERR90'

            call h5dwrite_f (dset_id, NumType, Array3D,&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR100'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )
            
            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR110'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR120'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR130'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowI4_3D

    !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowI4_2D(HDF5ID, GroupName, Name,                            &
                                    Array2D, OffSet1, OffSet2,                          &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer, dimension(:, :)   , pointer            :: Array2D
        integer, optional                               :: OffSet1
        integer, optional                               :: OffSet2
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            NumType = H5T_NATIVE_INTEGER 
            
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) AuxChar
                stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR20'
            endif
                !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR60'

            if (rank >3) then
                stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR70'
            endif

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR80'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                lower_bound(2) = Me%Limits%JLB
                upper_bound(2) = Me%Limits%JUB
            endif
            
            if (rank ==3) then
                lower_bound(3) = Me%Limits%KLB
                upper_bound(3) = Me%Limits%KUB
            endif

            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            if (Rank >= 2) then
                if (present(OffSet2)) then
                    offset_in(2) = OffSet2
                else                
                    offset_in(2) = lower_bound(2) - 1    
                endif            
            endif                

            if (rank == 3) then
                offset_in(3) = 0
            endif

            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR90'

            !Defines the memory dataspace
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            dims_mem(2) = upper_bound(2)-lower_bound(2) + 1
            
            if (rank == 3) then
                dims_mem(3) = 1
            endif
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR100'

            call h5dwrite_f (dset_id, NumType, Array2D,&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR110'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR120'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR130'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_2D - ModuleHDF5 - ERR140'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowI4_2D

    !--------------------------------------------------------------------------

    subroutine HDF5WriteWindowI4_1D(HDF5ID, GroupName, Name,                            &
                                    Array1D, OffSet1,                                   &
                                    OutputNumber, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        character(len=*)                                :: Name
        integer, dimension(:)      , pointer            :: Array1D
        integer, optional                               :: OffSet1
        integer, optional                               :: OutputNumber
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer, dimension(3)                           :: lower_bound, upper_bound
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HSIZE_T), dimension(7)                  :: dims_mem
        integer(HID_T)                                  :: dset_id, gr_id, datatype_id
        integer(HID_T)                                  :: space_id, rank
        integer(HID_T)                                  :: memspace_id, NumType, class_id
        integer(HSSIZE_T), dimension(:), allocatable    :: offset_in
        integer(HSIZE_T ), dimension(:), allocatable    :: count_in
        integer(HID_T)                                  :: STAT_CALL
        character(StringLength)                         :: AuxChar
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            NumType = H5T_NATIVE_INTEGER 
            
            !Opens the Group
            call h5gopen_f      (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR10'

            !Opens the DataSet
            if (present(OutputNumber)) then
                call ConstructDSName (Name, OutputNumber, AuxChar)
            else
                AuxChar = Name
            endif

            !Opens the Dataset
            call h5dopen_f      (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR20'

            !Gets the data type id
            call h5dget_type_f (dset_id, datatype_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR30'

            call h5tget_class_f(datatype_id, class_id, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR40'

            !Gets a handle of the dataspace
            call h5dget_space_f (dset_id, space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR50'

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR60'

            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR70'
            
            lower_bound(:) = FillValueInt
            upper_bound(:) = FillValueInt
            
            lower_bound(1) = Me%Limits%ILB
            upper_bound(1) = Me%Limits%IUB

            if (rank >=2) then
                stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR80'
            endif
            
            allocate (offset_in (rank))
            allocate (count_in  (rank))
            
            if (present(OffSet1)) then
                offset_in(1) = OffSet1
            else                
                offset_in(1) = lower_bound(1) - 1    
            endif            

            count_in (1:rank) = upper_bound(1:rank) - lower_bound(1:rank) + 1

            !Defines the hyperslab in the dataset
            call h5sselect_hyperslab_f        (space_id, H5S_SELECT_SET_F, offset_in, count_in, &
                                               STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR90'

            !Defines the memory dataspace
            dims_mem(1) = upper_bound(1)-lower_bound(1) + 1
            
            call h5screate_simple_f (rank, dims_mem, memspace_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR100'

            call h5dwrite_f (dset_id, NumType, Array1D,&
                            dims_mem, STAT_CALL, memspace_id, space_id)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_3D - ModuleHDF5 - ERR110'
            
            !Deallocates temporary matrixes
            deallocate (offset_in )
            deallocate (count_in  )

            !Closes data space
            call h5sclose_f  (space_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR120'

            !End access to the dataset
            call h5dclose_f  (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR130'

            !Closes group
            call h5gclose_f  (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5WriteWindowI4_1D - ModuleHDF5 - ERR140'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_        

    end subroutine HDF5WriteWindowI4_1D

    !--------------------------------------------------------------------------

    subroutine HDF5FlushMemory (HDF5ID, ErrorMessage, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        character(len=*), optional                  :: ErrorMessage
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: STAT_
        integer                                     :: ready_

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call h5fflush_f (Me%FileID, H5F_SCOPE_LOCAL_F, STAT_)
            if (STAT_ /= SUCCESS_) stop 'HDF5FlushMemory - ModuleHDF5 - ERR01'

            !Close and reopens file
            call h5fclose_f(Me%FileID,   HDFERR = STAT_)
            if (STAT_ /= SUCCESS_) stop 'HDF5FlushMemory - ModuleHDF5 - ERR02'
            
            if (trim(Me%FileName2)/= trim(Me%FileName)) then
                write(*,*) "FileName before",trim(Me%FileName)
                Me%FileName = trim(Me%FileName2)
            endif

            call h5fopen_f (trim(Me%FileName), ACCESS_FLAGS = H5F_ACC_RDWR_F,             &
                            FILE_ID = Me%FileID, HDFERR = STAT_)
            if (STAT_ /= SUCCESS_) then
                if (present(ErrorMessage)) then
                    write(*,*) trim(ErrorMessage)
                endif
                write(*,*) "FileName",trim(Me%FileName)
                write(*,*) "FileID", Me%FileID
                write(*,*) "HDF5ID", HDF5ID
                write(*,*) "STAT_", STAT_
                stop 'HDF5FlushMemory - ModuleHDF5 - ERR03'
            endif                

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine HDF5FlushMemory
    
    !--------------------------------------------------------------------------

    subroutine ConstructDSName (Name, OutputNumber, AuxChar)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Name
        integer                                     :: OutputNumber
        character(len=*)                            :: AuxChar

        !Local-----------------------------------------------------------------
        character(StringLength)                     :: AuxNum

        write(AuxNum, fmt=*)OutputNumber

        if     (OutputNumber < 10     ) then
            AuxChar = trim(adjustl(Name))//"_0000"//trim(adjustl(AuxNum))
        elseif (OutputNumber < 100    ) then
            AuxChar = trim(adjustl(Name))//"_000" //trim(adjustl(AuxNum))
        elseif (OutputNumber < 1000   ) then
            AuxChar = trim(adjustl(Name))//"_00"  //trim(adjustl(AuxNum))
        elseif (OutputNumber < 10000  ) then
            AuxChar = trim(adjustl(Name))//"_0"   //trim(adjustl(AuxNum))
        else
            AuxChar = trim(adjustl(Name))//"_"    //trim(adjustl(AuxNum))
        endif

    end subroutine ConstructDSName

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetHDF5FileID (HDF5ID, FileID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        integer                                     :: FileID    
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            FileID = Me%FileID
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine GetHDF5FileID

    !--------------------------------------------------------------------------
    
    

    subroutine GetHDF5FileName (HDF5ID, FileName, STAT)

        !Arguments-------------------------------------------------------------
        integer                  , intent(IN)       :: HDF5ID
        character(len=*)         , intent(OUT)      :: FileName    
        integer, optional        , intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            FileName = Me%FileName
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine GetHDF5FileName

    !--------------------------------------------------------------------------

    
    logical function GetHDF5FileOkToRead (FileName)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: FileName

        !Local-----------------------------------------------------------------
        integer                                     :: FileID, STAT_CALL    
        
        !Local-----------------------------------------------------------------
    
        call h5fopen_f (trim(FileName), ACCESS_FLAGS = H5F_ACC_RDONLY_F,                &
                        FILE_ID = FileID, HDFERR = STAT_CALL)    
                        
        if (STAT_CALL == SUCCESS_) then
        
            call h5fclose_f(FileID, HDFERR = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5FileOkToRead - ModuleHDF5 - ERR10'
            
            GetHDF5FileOkToRead = .true.
            
        else
            
            GetHDF5FileOkToRead = .false.
        
        endif                                                
                        
    end function GetHDF5FileOkToRead
    
    !--------------------------------------------------------------------------    

    subroutine GetHDF5FileAccess (HDF5_CREATE, HDF5_READ, HDF5_READWRITE)
#ifdef _GUI_
        !DEC$ IF DEFINED(_X86_)
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : '_GetHDF5FileAccess@12'  :: GetHDF5FileAccess
        !DEC$ ELSE
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : 'GetHDF5FileAccess'      :: GetHDF5FileAccess
        !DEC$ ENDIF
#endif

        !Arguments-------------------------------------------------------------
        integer, optional                               :: HDF5_CREATE
        integer, optional                               :: HDF5_READ
        integer, optional                               :: HDF5_READWRITE

        !Local-----------------------------------------------------------------

        if (present(HDF5_CREATE     )) HDF5_CREATE      = HDF5_CREATE_
        if (present(HDF5_READ       )) HDF5_READ        = HDF5_READ_
        if (present(HDF5_READWRITE  )) HDF5_READWRITE   = HDF5_READWRITE_

    end subroutine GetHDF5FileAccess
    
    !--------------------------------------------------------------------------

    subroutine GetHDF5GroupNumberOfItems (HDF5ID, GroupName, nItems, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        character(len=*)                            :: GroupName
        integer(HID_T), intent(OUT)                 :: nItems    
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer(HID_T)                              :: STAT_CALL
        integer(HID_T)                              :: gr_id
        

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call h5gopen_f      (Me%FileID, trim(adjustl(GroupName)), gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'GroupName= ', GroupName
                stop 'GetHDF5GroupNumberOfItems - ModuleHDF5 - ERR10'
            endif                
            
            call h5gn_members_f (gr_id, trim(adjustl(GroupName)), nItems, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'GroupName= ', GroupName
                stop 'GetHDF5GroupNumberOfItems - ModuleHDF5 - ERR20'
            endif

            call h5gclose_f     (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'GroupName= ', GroupName
                stop 'GetHDF5GroupNumberOfItems - ModuleHDF5 - ERR30'
            endif
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine GetHDF5GroupNumberOfItems
    
    !--------------------------------------------------------------------------


    subroutine GetHDF5ArrayDimensions (HDF5ID, GroupName, ItemName, OutputNumber, NDim, Imax, Jmax, Kmax, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        character(len=*)                            :: GroupName
        character(len=*)                            :: ItemName
        integer,  intent(IN ), optional             :: OutputNumber
        integer,  intent(OUT), optional             :: NDim
        integer,  intent(OUT), optional             :: Imax, Jmax, Kmax
        integer,  intent(OUT), optional             :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: gr_id, space_id, dset_id, rank
        integer(HSIZE_T), dimension(7)              :: dims, maxdims


        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        character(StringLength)                     :: ItemName_
        

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (ItemName, OutputNumber, ItemName_)
            else
                ItemName_ = ItemName
            endif
            
           !Opens the group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'Group name not present in the hdf5 input file', GroupName
                stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR10'
            endif
            
            !Opens the Dataset
            call h5dopen_f          (gr_id, ItemName_, dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR20'
            
            !Opens data space
            call h5dget_space_f     (dset_id, space_id,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR30'

           !Gets rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR35'

            !Gets dims
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims,  STAT_CALL)
            if (STAT_CALL < SUCCESS_) stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR40'

            !Closes data space
            call h5sclose_f         (space_id,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR50'

            !Closes data set
            call h5dclose_f         (dset_id,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR60'

            !Closes the Group
            call h5gclose_f         (gr_id,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR70'
            
            if (present(Imax)) then
                Imax = dims(1)
                if (rank <1) then
                    stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR80'
                endif
            endif                
            if (present(Jmax)) then
                Jmax = dims(2)
                if (rank <2) then
                    stop 'GetHDF5ArrayDimensions - ModuleHDF5 - ERR90'
                endif
            endif              
            if (present(Kmax)) then
                Kmax = dims(3)
                if (rank <3) then
                    Kmax = 1
                endif
            endif    
            
            if (present(NDim)) then
                NDim = rank    
            endif          

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine GetHDF5ArrayDimensions
    
    !--------------------------------------------------------------------------

    integer function GetHDF5ArrayDim (HDF5ID, GroupName, ItemName, OutputNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        character(len=*)                            :: GroupName
        character(len=*)                            :: ItemName
        integer,  intent(IN ), optional             :: OutputNumber
        integer,  intent(OUT), optional             :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: gr_id, space_id, dset_id, rank

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        character(StringLength)                     :: ItemName_
        

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !Creates the dataset with default properties
            if (present(OutputNumber)) then
                call ConstructDSName (ItemName, OutputNumber, ItemName_)
            else
                ItemName_ = ItemName
            endif
            
           !Opens the group
            call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDim - ModuleHDF5 - ERR10'

            !Opens the Dataset
            call h5dopen_f          (gr_id, ItemName_, dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDim - ModuleHDF5 - ERR20'
            
            !Opens data space
            call h5dget_space_f     (dset_id, space_id,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDim - ModuleHDF5 - ERR30'

           !Gets rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDim - ModuleHDF5 - ERR35'

            !Closes data space
            call h5sclose_f         (space_id,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDim - ModuleHDF5 - ERR50'

            !Closes data set
            call h5dclose_f         (dset_id,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDim - ModuleHDF5 - ERR60'

            !Closes the Group
            call h5gclose_f         (gr_id,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ArrayDim - ModuleHDF5 - ERR70'
            
            GetHDF5ArrayDim = rank
            
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end function GetHDF5ArrayDim 
    
    !--------------------------------------------------------------------------

         

    subroutine GetHDF5GroupID (HDF5ID, FatherGroupName, GroupPosition, GroupName, &
                               Units, Rank, Dimensions, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: FatherGroupName
        character(len=*)                                :: GroupName
        integer                                         :: GroupPosition
        character(len=*), optional                      :: Units
        integer, optional,               intent(OUT)    :: Rank
        integer, optional, dimension(7), intent(OUT)    :: Dimensions
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HID_T)                                  :: STAT_CALL, rank_
        integer(HID_T)                                  :: gr_id
        integer(HID_T)                                  :: dset_id
        integer(HID_T)                                  :: space_id
        !integer(HID_T)                                  :: class_id, size, datatype_id
        integer(HSIZE_T), dimension(7)                  :: dims, maxdims
        integer(HID_T)                                  :: GroupType
        integer(HID_T)                                  :: attr_id, type_id
        character(len=StringLength)                     :: Units_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
           
            call h5gopen_f       (Me%FileID, trim(adjustl(FatherGroupName)), gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'GetHDF5GroupID - ModuleHDF5 - ERR10'
            endif
                    
            call h5gget_obj_info_idx_f(Me%FileID, trim(adjustl(FatherGroupName)), GroupPosition-1, &   
                                       GroupName, GroupType, STAT_CALL)

            if (GroupType == H5G_DATASET_F) then

                !Opens data set
                call h5dopen_f      (gr_id, trim(adjustl(GroupName)), dset_id, STAT_CALL)

                !Opens data space
                call h5dget_space_f (dset_id, space_id, STAT_CALL)
    
                !Gets rank
                call h5sget_simple_extent_ndims_f (space_id, rank_, STAT_CALL)

                if(present(Rank))then
                    Rank = rank_
                endif

                !Gets datatype
                !call h5dget_type_f (dset_id, datatype_id,   STAT)
                !call h5tget_size_f (datatype_id, size,      STAT)
                !call h5tget_class_f(datatype_id, class_id,  STAT) 

                !call h5tclose_f    (datatype_id, STAT) 

                !Gets dims
                call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL)

                if(present(Dimensions))then
                    Dimensions = dims
                endif

                !Reads Units
                if(present(Units))then
                
                    call h5aopen_name_f     (dset_id, "Units", attr_id, STAT_CALL)
                    call h5Tcopy_f          (H5T_NATIVE_CHARACTER, type_id, STAT_CALL)
                    call h5Tset_size_f      (type_id, int(StringLength,8), STAT_CALL)
                    call h5aread_f          (attr_id, type_id, Units_, dims, STAT_CALL)
                    call h5aclose_f         (attr_id, STAT_CALL) 
                    call h5Tclose_f         (type_id, STAT_CALL)
                    
                    Units = Units_
                    
                endif


                !Closes data space
                call h5sclose_f     (space_id, STAT_CALL)

                !Closes data set
                call h5dclose_f     (dset_id, STAT_CALL)
            
            end if
            
            call h5gclose_f      (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5GroupID - ModuleHDF5 - ERR30'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine GetHDF5GroupID

    !--------------------------------------------------------------------------

    subroutine GetHDF5DataTypeID (HDF5ID, FatherGroupName, GroupPosition, GroupName, &
                                  DataType, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(*)                                    :: FatherGroupName
        character(*)                                    :: GroupName
        integer                                         :: GroupPosition
        integer, intent(OUT)                            :: DataType
        integer, optional, intent(OUT)                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: gr_id
        integer(HID_T)                                  :: dset_id
        integer(HID_T)                                  :: space_id
        integer(HID_T)                                  :: datatype_id
        integer(HID_T)                                  :: GroupType
        integer(HSIZE_T)                                :: size
        integer(HID_T)                                  :: class_id

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
           
            call h5gopen_f       (Me%FileID, adjustl(trim(FatherGroupName)), gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5DataTypeID - ModuleHDF5 - ERR01'
                    
            call h5gget_obj_info_idx_f(Me%FileID, adjustl(trim(FatherGroupName)), GroupPosition-1, &   
                                       GroupName, GroupType, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5DataTypeID - ModuleHDF5 - ERR02'
            
            if (GroupType == H5G_DATASET_F) then

                !Opens data set
                call h5dopen_f      (gr_id, adjustl(trim(GroupName)), dset_id, STAT_CALL)

                !Opens data space
                call h5dget_space_f (dset_id, space_id, STAT_CALL)
    
                !Gets datatype
                call h5dget_type_f (dset_id, datatype_id,   STAT)
                call h5tget_size_f (datatype_id, size,      STAT)
                call h5tget_class_f(datatype_id, class_id,  STAT) 

if1 :           if (class_id == H5T_FLOAT_F) then
if11 :              if (size == 8) then
                        dataType = H5T_NATIVE_DOUBLE
                    elseif (size == 4) then if11
                        dataType = H5T_NATIVE_REAL
                    endif if11
                elseif (class_id == H5T_INTEGER_F) then if1
                    dataType = H5T_NATIVE_INTEGER
                else if1
                    stop 'GetHDF5DataTypeID - ModModuleHDF5uleDDC - ERR03'
                endif if1

                call h5tclose_f    (datatype_id, STAT) 

                !Closes data space
                call h5sclose_f     (space_id, STAT_CALL)

                !Closes data set
                call h5dclose_f     (dset_id, STAT_CALL)
            
            end if
            
            call h5gclose_f      (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5DataTypeID - ModuleHDF5 - ERR04'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine GetHDF5DataTypeID

    !--------------------------------------------------------------------------

    subroutine GetHDF5GroupExist (HDF5ID, GroupName, Exist, nGroup, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: GroupName
        logical                                         :: Exist
        integer, optional                               :: nGroup
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(4)                                      :: STAT_CALL
        integer(HID_T)                                  :: nmembers

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
           

            !Turns Error printing of
            call h5eset_auto_f  (0, STAT_CALL)

            call h5gn_members_f (Me%FileID, trim(GroupName), nmembers, STAT_CALL)

            if (STAT_CALL == -1) then
                Exist = .false.
            else
                Exist = .true.
            endif

            if (present(nGroup)) nGroup = nmembers 
            
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine GetHDF5GroupExist

    !--------------------------------------------------------------------------
    
    subroutine GetHDF5DataSetExist (HDF5ID, DataSetName, Exist, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: DataSetName
        logical                                         :: Exist
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(4)                                      :: STAT_CALL
        integer(HID_T)                                  :: dset_id

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
           

            !Turns Error printing of
            call h5eset_auto_f  (0, STAT_CALL)

            call h5dopen_f(Me%FileID, trim(DataSetName), dset_id, STAT_CALL) 

            if (STAT_CALL == -1) then
                Exist = .false.
            else
                !Closes data set
                call h5dclose_f (dset_id, STAT)

                Exist = .true.
            endif
            
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine GetHDF5DataSetExist


    !--------------------------------------------------------------------------
 
     
    subroutine GetHDF5ObjectInfo (HDF5ID, FatherGroupName, GroupPosition, GroupName, &
                                  GroupType, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: HDF5ID
        character(len=*)                                :: FatherGroupName
        integer                                         :: GroupPosition
        character(len=*)                                :: GroupName
        integer(HID_T)                                  :: GroupType
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer(HID_T)                                  :: STAT_CALL
        integer(HID_T)                                  :: gr_id

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
       
            call h5gopen_f       (Me%FileID, trim(adjustl(FatherGroupName)), gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ObjectInfo - ModuleHDF5 - ERR10'
                    
            call h5gget_obj_info_idx_f(gr_id, FatherGroupName, GroupPosition-1, &   
                                       GroupName, GroupType, STAT_CALL)
            
            call h5gclose_f      (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ObjectInfo - ModuleHDF5 - ERR30'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine GetHDF5ObjectInfo

    !--------------------------------------------------------------------------

    subroutine HDF5ReadGenericRealAttribute (HDF5ID, GroupName, ItemName, ItemType, AttributeName, ValueReal, STAT)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: HDF5ID
        character(len=*)                            :: GroupName
        character(len=*)                            :: ItemName
        integer                                     :: ItemType
        character(len=*)                            :: AttributeName
        real                                        :: ValueReal
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer(HID_T)                              :: STAT_CALL        
        integer(HID_T)                              :: gr_id, attr_id, dset_id
        integer(HSIZE_T), dimension(7)              :: dims

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !If item is a group
            if (ItemType == TypeSDS) then

                !Opens the group
                call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR10'

                !Opens the Dataset
                call h5dopen_f          (gr_id, ItemName, dset_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR20'

                !Reads Real Value
                call h5aopen_name_f     (dset_id, trim(AttributeName), attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR30'
                
                call h5aread_f          (attr_id, H5T_NATIVE_REAL, ValueReal, dims, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR40'
                
                call h5aclose_f         (attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR50'

                !Closes the Dataset
                call h5dclose_f        (dset_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR60'
                

                !Closes the Group
                call h5gclose_f         (gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR70'

            endif


            if (ItemType == TypeVG) then

                !Opens the group
                call h5gopen_f (Me%FileID, trim(GroupName)//"/"//trim(ItemName), gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR80'

                !Reads Real Value
                call h5aopen_name_f     (dset_id, trim(AttributeName), attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR90'
                
                call h5aread_f          (attr_id, H5T_NATIVE_REAL, ValueReal, dims, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR90'
                call h5aclose_f         (attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR100'

                !Closes the Group
                call h5gclose_f         (gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadGenericRealAttribute - ModuleHDF5 - ERR110'

            endif
 
             STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_
        
    end subroutine HDF5ReadGenericRealAttribute

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine HDF5UpdateGenericRealAttribute(  HDF5ID, GroupName, ItemName,     &
                                                ItemType, AttributeName,         &
                                                ValueReal, STAT)

        !Arguments-------------------------------------------------------------
        integer(HID_T), intent(IN)                                  :: HDF5ID
        character(len=*), intent(IN)                                :: GroupName
        character(len=*) , intent(IN)                               :: ItemName
        integer, intent(IN)                                         :: ItemType
        character(len=*), intent(IN)                                :: AttributeName
        real, intent(IN)                                            :: ValueReal
        integer, optional, intent(OUT)                              :: STAT

        !Local-----------------------------------------------------------------
        integer                                                     :: STAT_, ready_
        integer(HID_T)                                              :: STAT_CALL        
        integer(HID_T)                                              :: gr_id
        integer(HID_T)                                              :: attr_id
        integer(HID_T)                                              :: dset_id
        integer(HSIZE_T), dimension(7)                              :: dims
        real                                                        ::OldValue

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !If item is a group
            if (ItemType == H5G_DATASET_F) then

                !Opens the group
                call h5gopen_f (Me%FileID, GroupName, gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR10'

                !Opens the Dataset
                call h5dopen_f          (gr_id, ItemName, dset_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR20'

                !Opens attribute
                call h5aopen_name_f     (dset_id, trim(AttributeName), attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR30'
                
                !Reads OldValue
                call h5aread_f      (attr_id, H5T_NATIVE_REAL, OldValue, dims, STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR40'
                
                if (trim(AttributeName) .EQ."Minimum") then
                
                    if (ValueReal < OldValue .OR. OldValue .EQ. null_real .OR. OldValue .EQ. null_int) then
                        call h5awrite_f (attr_id, H5T_NATIVE_REAL, ValueReal, dims, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR50a'
                    endif
                    
                elseif (trim(AttributeName) .EQ."Maximum") then
                
                    if (ValueReal > OldValue .OR. OldValue .EQ. null_real .OR. OldValue .EQ. null_int) then
                        call h5awrite_f (attr_id, H5T_NATIVE_REAL, ValueReal, dims, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR50b'
                    endif
                
                endif

                !Closes attribute 
                call h5aclose_f         (attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR60'

                !Closes the Dataset
                call h5dclose_f        (dset_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR70'
                
                !Closes the Group
                call h5gclose_f         (gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR80'

            endif

            if (ItemType == H5G_GROUP_F) then

                !Opens the group
                call h5gopen_f (Me%FileID, trim(GroupName)//"/"//trim(ItemName), gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR90'

                !Opens attribute
                call h5aopen_name_f     (dset_id, trim(AttributeName), attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR100'
                
                !Reads OldValue
                call h5aread_f      (attr_id, H5T_NATIVE_REAL, OldValue, dims, STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR110'
                
                if (trim(AttributeName) .EQ."Minimum") then
                
                    if (ValueReal < OldValue .OR. OldValue .EQ. null_real .OR. OldValue .EQ. null_int) then
                        call h5awrite_f (attr_id, H5T_NATIVE_REAL, ValueReal, dims, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR50a'
                    endif
                    
                elseif (trim(AttributeName) .EQ."Maximum") then
                
                    if (ValueReal > OldValue .OR. OldValue .EQ. null_real .OR. OldValue .EQ. null_int) then
                        call h5awrite_f (attr_id, H5T_NATIVE_REAL, ValueReal, dims, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR50b'
                    endif
                
                endif

                !Closes attribute 
                call h5aclose_f         (attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR130'

                !Closes the Dataset
                call h5aclose_f         (attr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR140'

                !Closes the Group
                call h5gclose_f         (gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5UpdateGenericRealAttribute - ModuleHDF5 - ERR150'

            endif
 
             STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_
        
    end subroutine HDF5UpdateGenericRealAttribute

    !--------------------------------------------------------------------------
      
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillHDF5 (HDF5ID, STAT)
#ifdef _GUI_
        !DEC$ IF DEFINED(_X86_)
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : '_KillHDF5@8'   :: KillHDF5
        !DEC$ ELSE
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : 'KillHDF5'      :: KillHDF5
        !DEC$ ENDIF
#endif

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_    
        integer(4)                                  :: STAT_CALL
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (HDF5ID, ready_)

        if (ready_ .NE. OFF_ERR_) then
            
            nUsers = DeassociateInstance(mHDF5_,  Me%InstanceID)

            if (nUsers == 0) then
                
                call h5fclose_f(Me%FileID, HDFERR = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillHDF5 - ModuleHDF5 - ERR01'
                
                call DeallocateAuxMatrixes

                !Deallocates Instance of HDF5
                call DeallocateInstance

                HDF5ID = 0
            
                STAT_ = SUCCESS_

            endif

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine KillHDF5

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HDF5), pointer                      :: AuxHDF5
        type (T_HDF5), pointer                      :: PreviousHDF5
        
        !Updates pointers
        if (Me%InstanceID == FirstHDF5%InstanceID) then
            FirstHDF5 => FirstHDF5%Next
        else
            PreviousHDF5 => FirstHDF5
            AuxHDF5      => FirstHDF5%Next
            do while (AuxHDF5%InstanceID /= Me%InstanceID)
                PreviousHDF5 => AuxHDF5
                AuxHDF5      => AuxHDF5%Next
            enddo

            !Now update linked list
            PreviousHDF5%Next => AuxHDF5%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------
    
    subroutine DeallocateAuxMatrixes
    
    
    !Begin------------------------------------------------------------------------------
    
        if (associated(Me%AuxMatrixes%DataR4_1D)) deallocate(Me%AuxMatrixes%DataR4_1D)
        if (associated(Me%AuxMatrixes%DataR4_2D)) deallocate(Me%AuxMatrixes%DataR4_2D)
        if (associated(Me%AuxMatrixes%DataR4_3D)) deallocate(Me%AuxMatrixes%DataR4_3D)

        if (associated(Me%AuxMatrixes%DataR8_1D)) deallocate(Me%AuxMatrixes%DataR8_1D)
        if (associated(Me%AuxMatrixes%DataR8_2D)) deallocate(Me%AuxMatrixes%DataR8_2D)
        if (associated(Me%AuxMatrixes%DataR8_3D)) deallocate(Me%AuxMatrixes%DataR8_3D)

        if (associated(Me%AuxMatrixes%DataI4_1D)) deallocate(Me%AuxMatrixes%DataI4_1D)
        if (associated(Me%AuxMatrixes%DataI4_2D)) deallocate(Me%AuxMatrixes%DataI4_2D)
        if (associated(Me%AuxMatrixes%DataI4_3D)) deallocate(Me%AuxMatrixes%DataI4_3D)

    end subroutine DeallocateAuxMatrixes    


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Ready (HDF5ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (HDF5ID > 0) then
            call LocateObjHDF5 (HDF5ID)
            ready_ = VerifyReadLock (mHDF5_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1
    
        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjHDF5 (HDF5ID)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5ID

        !Local-----------------------------------------------------------------

        Me => FirstHDF5
        do while (associated (Me))
            if (Me%InstanceID == HDF5ID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) then
            write(*,*) Me%FileName
            stop 'ModuleHDF5 - LocateObjHDF5 - ERR01'
        endif            

    end subroutine LocateObjHDF5

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GUI GU

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !For GUI only
#ifdef _GUI_
    subroutine HDF5InquireFile (FileName, hInstance, handleTV, STAT)
        !DEC$ IF DEFINED(_X86_)
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : '_HDF5InquireFile@16'  :: HDF5InquireFile
        !DEC$ ELSE
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : 'HDF5InquireFile'      :: HDF5InquireFile
        !DEC$ ENDIF

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: FileName
        integer                                     :: hInstance
        integer                                     :: handleTV
        integer                                     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: file_id
        integer                                     :: hParent, ret

        !Initializes predefined datatypes
        call h5open_f   (STAT)
        if (STAT /= SUCCESS_) return

        !Open the file as read only
        call h5fopen_f  (trim(FileName), H5F_ACC_RDONLY_F, file_id, STAT)
        if (STAT /= SUCCESS_) return

        !Creates tree view image list
        ret = ConstructTreeViewImageLists (hInstance, handleTV)

        !Adds the name of the file as root elemente
        hParent = AddItemToHDFTree(handleTV, 0, trim(FileName), trim(FileName), "/", TypeFile)

        !Iterates through all subgroups
        call InquireSubGroup (file_id, "/", 1, hParent, handleTV, trim(FileName))
                      
        !Closes the file
        call h5fclose_f (file_id, STAT)
        if (STAT /= SUCCESS_) return

    end subroutine HDF5InquireFile

    !--------------------------------------------------------------------------

    recursive subroutine InquireSubGroup (ID, GroupName, Level, hParent, handleTV, FileName)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: ID
        character(len=*)                            :: GroupName
        integer                                     :: Level
        integer                                     :: hParent
        integer                                     :: handleTV
        character(len=*)                            :: FileName

        !Local-----------------------------------------------------------------
        integer                                     :: nmembers
        character(StringLength)                     :: obj_name
        integer                                     :: obj_type, idx
        integer(HID_T)                              :: gr_id, dset_id
        integer(HID_T)                              :: space_id, datatype_id, class_id, size
        integer                                     :: STAT
        character(StringLength)                     :: NewGroupName
        integer(HSIZE_T), dimension(7)              :: dims, maxdims
        integer                                     :: NewParent
        integer                                     :: Rank    


        !Get the number of members in the Group
        call h5gn_members_f(ID, GroupName, nmembers, STAT)
        if (STAT /= SUCCESS_) return
    
        do idx = 1, nmembers

            !Gets information about the group
            call h5gget_obj_info_idx_f(ID, GroupName, idx-1, obj_name, obj_type, STAT)
            if (STAT /= SUCCESS_) return
!            write(*,*)("+", i=1,Level),trim(adjustl(obj_name))

            if     (obj_type == H5G_DATASET_F) then

                !Opens data set
                call h5dopen_f      (ID, trim(adjustl(obj_name)), dset_id, STAT)

                !Opens data space
                call h5dget_space_f (dset_id, space_id, STAT)
    
                !Gets rank
                call h5sget_simple_extent_ndims_f (space_id, rank, STAT)

                !Gets datatype
                call h5dget_type_f (dset_id, datatype_id,   STAT)
                call h5tget_size_f (datatype_id, size,      STAT)
                call h5tget_class_f(datatype_id, class_id,  STAT) 

                call h5tclose_f    (datatype_id, STAT) 

                !Gets dims
                call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT) 

                !Closes data space
                call h5sclose_f     (space_id, STAT)

                !Closes data set
                call h5dclose_f (dset_id, STAT)
                
                NewParent = AddItemToHDFTree(handleTV, hParent, obj_name, trim(FileName),  &
                                             GroupName, TypeSDS, Rank, class_id, size, Dims)

            elseif (obj_type == H5G_GROUP_F  ) then

                NewParent = AddItemToHDFTree(handleTV, hParent, obj_name, trim(FileName),      &
                                             GroupName, TypeVG)

            endif


            !Looks for attributes if object is an DATASET
!            if (obj_type == H5G_DATASET_F) then

!                call h5dopen_f          (ID, trim(adjustl(obj_name)), dset_id, STAT)

                !Reads Minimum
!                call h5aopen_name_f     (dset_id, "Minimum", attr_id, STAT) 
!                call h5aread_f          (attr_id, H5T_NATIVE_REAL, Minimum, dims, STAT)
!                call h5aclose_f         (attr_id, STAT) 

                !Reads Maximum
!                call h5aopen_name_f     (dset_id, "Maximum", attr_id, STAT) 
!                call h5aread_f          (attr_id, H5T_NATIVE_REAL, Maximum, dims, STAT)
!                call h5aclose_f         (attr_id, STAT) 

                !Reads Units
!                call h5aopen_name_f     (dset_id, "Units", attr_id, STAT)
!                call h5Tcopy_f          (H5T_NATIVE_CHARACTER, type_id, STAT)
!                call h5Tset_size_f      (type_id, StringLength, STAT)
!                call h5aread_f          (attr_id, type_id, Units, dims, STAT)
!                call h5aclose_f         (attr_id, STAT) 
!                call h5Tclose_f         (type_id, STAT)

!                write(*,*)("|-", i=1,Level),"Minimum :", Minimum
!                write(*,*)("|-", i=1,Level),"Maximum :", Maximum
!                write(*,*)("|-", i=1,Level),"Units   :", trim(adjustl(Units))

!                call h5dclose_f        (dset_id, STAT)

!            endif

            !Looks for attributes if object is an GROUP
!            if (obj_type == H5G_GROUP_F) then

!                call h5gopen_f          (ID, trim(adjustl(obj_name)), gr_id, STAT)

                !Reads Minimum
!                call h5aopen_name_f     (gr_id, "Minimum", attr_id, STAT) 
!                call h5aread_f          (attr_id, H5T_NATIVE_REAL, Minimum, dims, STAT)
!                call h5aclose_f         (attr_id, STAT) 

                !Reads Maximum
!                call h5aopen_name_f     (gr_id, "Maximum", attr_id, STAT) 
!                call h5aread_f          (attr_id, H5T_NATIVE_REAL, Maximum, dims, STAT)
!                call h5aclose_f         (attr_id, STAT) 

!                write(*,*)("|-", i=1,Level),"Minimum :", Minimum
!                write(*,*)("|-", i=1,Level),"Maximum :", Maximum

!                call h5gclose_f        (gr_id, STAT)

!            endif

            !Looks for futher subgroups
            if (obj_type == H5G_GROUP_F) then
                if (GroupName == "/") then
                    NewGroupName = GroupName//trim(adjustl(obj_name))
                else
                    NewGroupName = GroupName//"/"//trim(adjustl(obj_name))
                endif
                call h5gopen_f       (ID, trim(adjustl(NewGroupName)), gr_id, STAT)
                call InquireSubGroup (gr_id, trim(adjustl(NewGroupName)), Level + 1, NewParent, handleTV, trim(FileName))
                call h5gclose_f      (gr_id, STAT)
            endif
            
        enddo

    end subroutine InquireSubGroup

    !--------------------------------------------------------------------------

    integer function AddItemToHDFTree(hTreeViewControl, hParent, ItemName, FileName,     &
                                      GroupName, nType, Rank, NumType, size, Dims)


        !Arguments-------------------------------------------------------------
        integer, intent(in)                         :: hTreeViewControl
        integer, intent(in)                         :: hParent
        character(len=*)                            :: ItemName
        character(len=*)                            :: FileName
        character(len=*)                            :: GroupName
        integer, intent(in)                         :: nType
        integer, optional                           :: Rank
        integer, optional                           :: NumType
        integer, optional                           :: size
        integer(HSIZE_T), dimension(7), optional    :: dims


        !Local-----------------------------------------------------------------
        type (T_TVITEM)                             :: TVItem
        type (T_TVINSERTSTRUCT)                     :: TVInsertStruct
        type (T_HDF5_DATA_ITEM), pointer            :: HDF5_DataItem

        !Sets the mask
        TVItem%mask = IOR(TVIF_TEXT, IOR(TVIF_IMAGE, IOR(TVIF_SELECTEDIMAGE, TVIF_PARAM)))
 
        !Sets the text of the item
        TVItem%pszText    = loc(trim(ItemName)//""C)

        !Stores HDF Data
        allocate(HDF5_DataItem)
        HDF5_DataItem%ItemType  = nType
        HDF5_DataItem%FileName  = trim(FileName)
        HDF5_DataItem%GroupName = GroupName
        HDF5_DataItem%ItemName  = ItemName
        if (present(Rank)) then
            HDF5_DataItem%Rank = Rank
        else
            HDF5_DataItem%Rank = -1
        endif
        if (present(NumType)) then
            if      (NumType == H5T_INTEGER_F .and. size == 4) then
                HDF5_DataItem%NumType = 1
            elseif  (NumType == H5T_FLOAT_F   .and. size == 4) then
                HDF5_DataItem%NumType = 2
            elseif  (NumType == H5T_INTEGER_F .and. size == 8) then
                HDF5_DataItem%NumType = 3
            elseif  (NumType == H5T_FLOAT_F   .and. size == 8) then
                HDF5_DataItem%NumType = 4
            endif
        else
            HDF5_DataItem%NumType = -1
        endif
        if (present(Dims)) then
            HDF5_DataItem%Dims = Dims
        else
            HDF5_DataItem%Dims = -1
        endif
        TVItem%lParam           = loc(HDF5_DataItem)

        select case (nType)
        !File
        case (TypeFile)
            TVItem%iImage               = IconFile
            TVItem%iSelectedImage       = IconFile
            TVInsertStruct%hParent      = TVI_ROOT
            TVInsertStruct%hInsertAfter = TVI_FIRST

        !VG
        case (TypeVG)
            TVItem%iImage               = IconVG1
            TVItem%iSelectedImage       = IconVG
            TVInsertStruct%hParent      = hParent
            TVInsertStruct%hInsertAfter = TVI_LAST

        !SDS
        case (TypeSDS)
            TVItem%iImage               = IconSDS1
            TVItem%iSelectedImage       = IconSDS
            TVInsertStruct%hParent      = hParent
            TVInsertStruct%hInsertAfter = TVI_LAST

        !Attribute
        case (TypeAttr)
            TVItem%iImage               = IconAttr
            TVItem%iSelectedImage       = IconAttr
            TVInsertStruct%hParent      = hParent
            TVInsertStruct%hInsertAfter = TVI_LAST

        end select

        TVInsertStruct%Item = TVItem

        !Add the item to the tree view control
        AddItemToHDFTree = SendMessage(hTreeViewControl, TVM_INSERTITEM, 0, loc(TVInsertStruct))

    end function AddItemToHDFTree

    !--------------------------------------------------------------------------

    logical function ConstructTreeViewImageLists(ghInstance, hTreeViewControl)

        !Arguments-------------------------------------------------------------
        integer, intent(in)     :: ghInstance
        integer, intent(in)     :: hTreeViewControl

        !Local-----------------------------------------------------------------
        integer                 :: hIml     !handle to the image list
        integer                 :: hbmp     !handle to the bitmap list
        integer                 :: ret

        !Create the image list
        hIml = ImageList_Create(16, 16, ILC_COLOR, 6, 0)

        !Load icon 1 (File)
        hbmp             = LoadBitMap   (ghInstance, LOC("ICONFILE"C))
        IconFile         = ImageList_Add(hIml, hbmp, 0)
        ret              = DeleteObject (hbmp)

        !Load icon 2 (VG)
        hbmp             = LoadBitMap   (ghInstance, LOC("ICONVG"C))
        IconVG           = ImageList_Add(hIml, hbmp, 0)
        ret              = DeleteObject (hbmp)

        !Load icon 3 (SDS)
        hbmp             = LoadBitMap   (ghInstance, LOC("ICONSDS"C))
        IconSDS          = ImageList_Add(hIml, hbmp, 0)
        ret              = DeleteObject (hbmp)

        !Load icon 4 (Attribute)
        hbmp             = LoadBitMap   (ghInstance, LOC("ICONATTR"C))
        IconAttr         = ImageList_Add(hIml, hbmp, 0)
        ret              = DeleteObject (hbmp)

        !Load icon 5 (Unselected VG)
        hbmp             = LoadBitMap   (ghInstance, LOC("ICONVG1"C))
        IconVG1          = ImageList_Add(hIml, hbmp, 0)
        ret              = DeleteObject (hbmp)

        !Load icon 6 (Unselected SDS)
        hbmp             = LoadBitMap   (ghInstance, LOC("ICONSDS1"C))
        IconSDS1         = ImageList_Add(hIml, hbmp, 0)
        ret              = DeleteObject (hbmp)

        !Send a message to the TreeViewControl
        ret = SendMessage(hTreeViewControl, TVM_SETIMAGELIST, TVSIL_NORMAL, hIml)

        ConstructTreeViewImageLists = TRUE

    end function ConstructTreeViewImageLists

    !--------------------------------------------------------------------------

    subroutine HDF5ReadAttributes (FileName, GroupName, ItemName, ItemType, Minimum, Maximum, Units, STAT)
        !DEC$ IF DEFINED(_X86_)
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : '_HDF5ReadAttributes@32'  :: HDF5ReadAttributes
        !DEC$ ELSE
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : 'HDF5ReadAttributes'      :: HDF5ReadAttributes
        !DEC$ ENDIF

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: FileName
        character(len=*)                            :: GroupName
        character(len=*)                            :: ItemName
        integer                                     :: ItemType
        real                                        :: Minimum
        real                                        :: Maximum
        character(len=*)                            :: Units
        integer                                     :: STAT

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: file_id, gr_id, attr_id, type_id, dset_id
        integer(HSIZE_T), dimension(7)              :: dims

        !Initializes predefined datatypes
        call h5open_f   (STAT)
        if (STAT /= SUCCESS_) return

        !Open the file as read only
        call h5fopen_f  (trim(FileName), H5F_ACC_RDONLY_F, file_id, STAT)
        if (STAT /= SUCCESS_) return

        !If item is a group
        if (ItemType == TypeSDS) then

            !Opens the group
            call h5gopen_f (file_id, GroupName, gr_id, STAT)
            if (STAT /= SUCCESS_) return

            !Opens the Dataset
            call h5dopen_f          (gr_id, ItemName, dset_id, STAT)

            !Reads Minimum
            call h5aopen_name_f     (dset_id, "Minimum", attr_id, STAT) 
            call h5aread_f          (attr_id, H5T_NATIVE_REAL, Minimum, dims, STAT)
            call h5aclose_f         (attr_id, STAT) 

            !Reads Maximum
            call h5aopen_name_f     (dset_id, "Maximum", attr_id, STAT) 
            call h5aread_f          (attr_id, H5T_NATIVE_REAL, Maximum, dims, STAT)
            call h5aclose_f         (attr_id, STAT) 

            !Reads Units
            call h5aopen_name_f     (dset_id, "Units", attr_id, STAT)
            call h5Tcopy_f          (H5T_NATIVE_CHARACTER, type_id, STAT)
            call h5Tset_size_f      (type_id, Int8(StringLength), STAT)
            call h5aread_f          (attr_id, type_id, Units, dims, STAT)
            call h5aclose_f         (attr_id, STAT) 
            call h5Tclose_f         (type_id, STAT)

            !Closes the Dataset
            call h5dclose_f        (dset_id, STAT)

            !Closes the Group
            call h5gclose_f         (gr_id, STAT)

        endif


        if (ItemType == TypeVG) then

            !Opens the group
            call h5gopen_f (file_id, trim(GroupName)//"/"//trim(ItemName), gr_id, STAT)
            if (STAT /= SUCCESS_) return

            !Reads Minimum
            call h5aopen_name_f     (gr_id, "Minimum", attr_id, STAT) 
            call h5aread_f          (attr_id, H5T_NATIVE_REAL, Minimum, dims, STAT)
            call h5aclose_f         (attr_id, STAT) 

            !Reads Maximum
            call h5aopen_name_f     (gr_id, "Maximum", attr_id, STAT) 
            call h5aread_f          (attr_id, H5T_NATIVE_REAL, Maximum, dims, STAT)
            call h5aclose_f         (attr_id, STAT) 

            Units = ""

            !Closes the Group
            call h5gclose_f         (gr_id, STAT)

        endif
 
        
                      
        !Closes the file
        call h5fclose_f (file_id, STAT)
        if (STAT /= SUCCESS_) return

    end subroutine HDF5ReadAttributes

    !--------------------------------------------------------------------------

    subroutine HDF5GetDimensions (FileName, GroupName, ItemName, dims, STAT) 

        !DEC$ IF DEFINED(_X86_)
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : '_HDF5GetDimensions@20'  :: HDF5GetDimensions
        !DEC$ ELSE
        !DEC$ ATTRIBUTES STDCALL, REFERENCE, ALIAS : 'HDF5GetDimensions'      :: HDF5GetDimensions
        !DEC$ ENDIF

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: FileName
        character(len=*)                            :: GroupName
        character(len=*)                            :: ItemName
        integer  (HSIZE_T), dimension(7)            :: dims, maxdims
        integer                                     :: STAT


        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: file_id, gr_id
        integer(HID_T)                              :: dset_id, space_id

        !Initializes predefined datatypes
        call h5open_f   (STAT)
        if (STAT /= SUCCESS_) return

        !Open the file as read only
        call h5fopen_f  (trim(FileName), H5F_ACC_RDONLY_F, file_id, STAT)
        if (STAT /= SUCCESS_) return

        !Opens the group
        call h5gopen_f (file_id, GroupName, gr_id, STAT)
        if (STAT /= SUCCESS_) return

        !Opens the Dataset
        call h5dopen_f          (gr_id, ItemName, dset_id, STAT)

        !Opens data space
        call h5dget_space_f     (dset_id, space_id, STAT)

        !Gets dims
        call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT) 

        !Closes data space
        call h5sclose_f         (space_id, STAT)

        !Closes data set
        call h5dclose_f         (dset_id, STAT)

        !Closes the Group
        call h5gclose_f         (gr_id, STAT)

        !Closes the file
        call h5fclose_f (file_id, STAT)
        if (STAT /= SUCCESS_) return


    end subroutine HDF5GetDimensions


#endif

end module ModuleHDF5

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
