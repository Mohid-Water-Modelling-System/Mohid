!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : GridData 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to Read / Write Data Associated to a grid
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

Module ModuleGridData

    use ModuleGlobalData
    use ModuleTime              
    use ModuleEnterData 
    use ModuleDrawing        
    use ModuleHorizontalGrid,   only: GetHorizontalGridSize, GetHorizontalGrid,         &
                                      GetGridCoordType, GetGridAngle, GetGridZone,      &
                                      GetLatitudeLongitude, GetGridOrigin,              &
                                      UnGetHorizontalGrid, GetCoordTypeList,            &
#ifndef _NO_HDF5                                      
                                      WriteHorizontalGrid,                              &
#endif
                                      GetCheckDistortion,                               &
                                      GetGridLatitudeLongitude, GetZCoordinates,        &
                                      GetDDecompParameters
                                      
#ifndef _NO_HDF5
    use ModuleHDF5,             only: ConstructHDF5, HDF5ReadData, GetHDF5FileAccess,   &
                                      GetHDF5GroupNumberOfItems, HDF5SetLimits,         &
                                      HDF5WriteData, KillHDF5
#endif                                      
    use ModuleStopWatch,        only: StartWatch, StopWatch
    use ModuleFunctions,        only: CHUNK_J


    implicit none

    private

    !Subroutine----------------------------------------------------------------

    !Contructor
    public  :: ConstructGridData
    private ::      AllocateInstance
    private ::      ReadGridDataFile
#ifndef _NO_HDF5    
    private ::      ReadFileEvolution
#endif

    !Selector
    public  :: GetGridDataFileName
    public  :: GetGridData
    public  :: GetGridData2DReference
    public  :: GetMaximumValue
    public  :: GetMaxValueInPolygon
    public  :: GetMinimumValue
    public  :: GetFillValue
    public  :: GetGridDataType
    public  :: GetIsGridData3D
    public  :: GetGridDataEvolution
    public  :: UngetGridData

    !Modifier
    public  :: WriteGridData
    public  :: ModifyGridData

    !Destructor
    public  :: KillGridData
    private ::      DeallocateInstance
#ifndef _NO_HDF5
    
    private ::      WriteFileEvolution
#endif

    !Management
    private ::      Ready

    !Interface
    interface WriteGridData
        module procedure WriteGridData_v1
        module procedure WriteGridData_v2
    end interface

    interface GetGridData
        module procedure GetGridData2D
        module procedure GetGridData3D
    end interface

    interface UnGetGridData
        module procedure UnGetGridData2D
        module procedure UnGetGridData3D
    end interface

    interface ModifyGridData
        module procedure ModifyGridData2DIncrement
        module procedure ModifyNewMatrixGridData2D        
        module procedure ModifyConstantGridData2D
    end interface

    !Parameter-----------------------------------------------------------------
    character(LEN = StringLength), parameter :: BeginGridData2D   = '<BeginGridData2D>'
    character(LEN = StringLength), parameter :: EndGridData2D     = '<EndGridData2D>'
    character(LEN = StringLength), parameter :: BeginGridData3D   = '<BeginGridData3D>'
    character(LEN = StringLength), parameter :: EndGridData3D     = '<EndGridData3D>'

    character(StringLength),       parameter :: Char_Bathymetry   = 'Bathymetry'

    !Type----------------------------------------------------------------------

    type T_Evolution

        type(T_Time)                        :: Now
        logical                             :: Yes            = .false.
        character(LEN = PathLength)         :: File           = null_str
        integer                             :: OldInstants    = null_int
        real,    dimension(:,:),   pointer  :: TimeInstants   => null()
        real,    dimension(:,:,:), pointer  :: GridData2D     => null()
        integer                             :: ObjHDF5        = 0
        integer                             :: ObjTime        = 0
        character(Len = StringLength)       :: PropName       = null_str

    end type T_Evolution
    
    private :: T_DDecomp
    type       T_DDecomp        
        logical                                 :: MasterOrSlave = .false. 
        type (T_Size2D)                         :: HaloMap
    end type T_DDecomp
    
    

    type      T_GridData
        integer                         :: InstanceID          = null_int
        character(LEN = PathLength)     :: FileName            = null_str
        character(LEN = StringLength)   :: COMENT1 = '******'         
        character(LEN = StringLength)   :: COMENT2 = '******'          
        real, dimension(:,:), pointer   :: GridData2D          => null()
        real, dimension(:,:,:), pointer :: GridData3D          => null()
        real, dimension(:,:), pointer   :: GridData2Dreference => null()
        real                            :: DefaultValue        = null_real
        logical                         :: ConstantInSpace     = .false.   
        real                            :: MaximumValue        = null_real
        real                            :: MinimumValue        = null_real
        logical                         :: Is3D                = .false.
        logical                         :: ReadFile            = .false.
        type (T_Size2D)                 :: WorkSize, Size
        type (T_Size2D)                 :: GlobalWorkSize
        integer                         :: KLB                 = null_int
        integer                         :: KUB                 = null_int
        real                            :: FillValue           = null_real
        integer                         :: TypeZUV             = null_int
        type (T_Evolution)              :: Evolution
        type (T_DDecomp)    :: DDecomp
        integer                         :: ObjHorizontalGrid   = 0
        type (T_GridData), pointer      :: Next
    end type T_GridData

    !Global Module Variables
    type (T_GridData), pointer          :: FirstGridData
    type (T_GridData), pointer          :: Me

    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructGridData(GridDataID, HorizontalGridID, TimeID, FileName,        &
                                 KLB, KUB, DefaultValue, InMatrix3D, InMatrix2D,        &
                                 STAT)
                      
        !Arguments-------------------------------------------------------------
        integer                     , intent(INOUT) :: GridDataID
        integer                     , intent(INOUT) :: HorizontalGridID
        integer,            optional, intent(IN )   :: TimeID
        character(LEN = *), optional, intent(IN )   :: FileName
        integer,            optional, intent(IN )   :: KLB, KUB
        real,               optional, intent(IN )   :: DefaultValue
        real, dimension(:,:),   pointer, optional   :: InMatrix2D
        real, dimension(:,:,:), pointer, optional   :: InMatrix3D
        integer,            optional, intent(OUT)   :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: STAT_, ready_   
        integer                                     :: i, j, k
        logical                                     :: exists

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mGridData_)) then
            nullify (FirstGridData)
            call RegisterModule (mGridData_) 
        endif

        call Ready(GridDataID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            nullify (Me%GridData2D)
            nullify (Me%GridData3D)
            nullify (Me%GridData2Dreference)

            if (present(KLB) .and. present(KUB)) then
                Me%Is3D = .true.
                Me%KLB  = KLB
                Me%KUB  = KUB
            else
                Me%Is3D = .false.
            endif
            
            Me%ReadFile =.true.
            
            if (Me%Is3D) then
                if (present(InMatrix3D)) then
                    Me%ReadFile =.false.
                endif
            else
                if (present(InMatrix2D)) then
                    Me%ReadFile =.false.
                endif
            endif
            
            if (Me%ReadFile) then

                !Gets File name to read
                if (present(FileName)) then
                    Me%FileName = FileName    
                else
                    call ReadFileName('IN_GRID_DATA', Me%FileName, Message = "GridData File",  &
                                      STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleGridData - ERR01' 
                end if

                inquire(FILE = Me%FileName, EXIST = exists)
                if (.NOT. exists) then
                    write(*, *) 
                    write(*, *) 'GridData file especified does not exists: ', trim(Me%FileName)
                    stop 'ConstructGridData - ModuleGridData - ERR02' 
                end if

            endif

            !Associates Horizontal Grid
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)

            if (present(TimeID)) then
                Me%Evolution%ObjTime = AssociateInstance (mTIME_, TimeID)
            endif

            if (present(DefaultValue)) then
                Me%DefaultValue = DefaultValue
            else
                Me%DefaultValue = 0.
            endif

            !Gets Size
            call GetHorizontalGridSize (Me%ObjHorizontalGrid,                           &
                                        Size            = Me%Size,                      &
                                        WorkSize        = Me%WorkSize,                  &
                                        GlobalWorkSize  = Me%GlobalWorkSize,            &
                                        STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleGridData - ERR03' 
            
            call GetDDecompParameters(HorizontalGridID = Me%ObjHorizontalGrid,                 &
                                                  MasterOrSlave    = Me%DDecomp%MasterOrSlave, &
                                                  HaloMap          = Me%DDecomp%HaloMap,       &
                                                  STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleGridData - ERR04' 


            !Reads the GridData part
            if (Me%ReadFile) then
                call ReadGridDataFile
            else
                if (Me%Is3D) then
                    allocate(Me%GridData3D(Me%Size%ILB:Me%Size%IUB,                     &
                                           Me%Size%JLB:Me%Size%JUB,                     &
                                                   KLB:        KUB))

                    Me%GridData3D(:,:,:) = InMatrix3D(:,:,:)                                                               
                else
                    allocate(Me%GridData2D(Me%Size%ILB:Me%Size%IUB,                     &
                                           Me%Size%JLB:Me%Size%JUB))

                    Me%GridData2D(:,:)   = InMatrix2D(:,:)                
                endif
            endif

            !Searches for the maximum and minimum Depth / Frank May 99
            Me%MaximumValue = -1.0 * ABS(null_real)
            Me%MinimumValue =        ABS(null_real)


            if (associated(Me%GridData2D)) then

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%GridData2D(i, j) .GT. Me%MaximumValue)               &
                        Me%MaximumValue = Me%GridData2D(i, j)
                    if (Me%GridData2D(i, j) .LT. Me%MinimumValue .and.          &
                        Me%GridData2D(i, j) .GT. Me%FillValue)                  &
                        Me%MinimumValue = Me%GridData2D(i, j)
                end do       
                end do        

                if (Me%Evolution%Yes) then

                    Me%GridData2Dreference(:,:) = Me%GridData2D(:,:)

                    call ReadFileEvolution

                endif


            endif

            if (associated(Me%GridData3D)) then

                do k = Me%KLB,       Me%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%GridData3D(i, j, k) .GT. Me%MaximumValue)               &
                        Me%MaximumValue = Me%GridData3D(i, j, k)
                    if (Me%GridData3D(i, j, k) .LT. Me%MinimumValue .and.          &
                        Me%GridData3D(i, j, k) .GT. Me%FillValue)                  &
                        Me%MinimumValue = Me%GridData3D(i, j, k)
                end do
                end do
                end do

            endif


            !Returns ID
            GridDataID    = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            
            stop 'ModuleGridData - ConstructGridData - ERR05' 

        end if 


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructGridData

    !--------------------------------------------------------------------------

    subroutine AllocateInstance 

        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        type (T_GridData), pointer            :: NewObjGridData
        type (T_GridData), pointer            :: PreviousObjGridData


        !Allocates new instance
        nullify  (NewObjGridData)
        allocate (NewObjGridData)
        nullify  (NewObjGridData%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstGridData)) then
            FirstGridData            => NewObjGridData
            Me                       => NewObjGridData
        else
            PreviousObjGridData      => FirstGridData
            Me                       => FirstGridData%Next
            do while (associated(Me))
                PreviousObjGridData  => Me
                Me                   => Me%Next
            enddo
            Me                       => NewObjGridData
            PreviousObjGridData%Next => NewObjGridData
        endif

        Me%InstanceID = RegisterNewInstance (mGridData_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadGridDataFile

        !Arguments-------------------------------------------------------------                                                    
                                                                                                     
        !Local-----------------------------------------------------------------
        real                                        :: DefaultValue
        integer                                     :: ObjEnterData = 0
        integer                                     :: STAT_CALL
        integer                                     :: flag
        character(len=StringLength)                 :: Char_TypeZUV, Message

        !----------------------------------------------------------------------

        !Opens File
        call ConstructEnterData(ObjEnterData, Me%FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadGridDataFile - ModuleGridData - ERR10'

       

        !Gets FillValue
        call GetData                (Me%FillValue, ObjEnterData, flag,                   &
                                     keyword      = 'FILL_VALUE',                        &
                                     ClientModule = 'ModuleGridData',                    &
                                     default      = -99.0,                               &
                                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR20'

        !Gets FillValue
        call GetData                (Char_TypeZUV, ObjEnterData, flag,                   &
                                     keyword      = 'TYPE_ZUV',                          &
                                     ClientModule = 'ModuleGridData',                    &
                                     default      = "Z",                                 &
                                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR30'

        !Gets if the bathymetry can change in time
        call GetData                (Me%Evolution%Yes, ObjEnterData, flag,               &
                                     keyword      = 'EVOLUTION',                         &
                                     ClientModule = 'ModuleGridData',                    &
                                     default      =.false.,                              &
                                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR40'

        if (Me%Evolution%Yes .and. Me%Evolution%ObjTime == 0) stop 'ReadGridDataFile - ModuleGridData - ERR50'

        if (Me%Evolution%Yes .and. Me%Is3D) stop 'ReadGridDataFile - ModuleGridData - ERR60'


        if (Me%Evolution%Yes) then

            !Gets if the bathymetry can change in time
            call GetData            (Me%Evolution%File, ObjEnterData, flag,             &
                                     keyword      = 'EVOLUTION_FILE',                   &
                                     ClientModule = 'ModuleGridData',                   &
                                     default      ='******.***',                        &
                                     STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR70'

            if (flag == 0)  then
            
                call ReadFileName('EVOLUTION_FILE', Me%Evolution%File,                  &
                                   Message = Message, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'ReadGridDataFile - ModuleGridData - ERR80'

            endif


            call GetData            (Me%Evolution%PropName, ObjEnterData, flag,         &
                                     keyword      = 'PROPERTY_NAME',                    &
                                     ClientModule = 'ModuleGridData',                   &
                                     default      = trim(Char_Bathymetry),              &
                                     STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR90'


        endif


        call GetData                (Me%ConstantInSpace, ObjEnterData, flag,            &
                                     keyword      = 'CONSTANT_IN_SPACE',                &
                                     ClientModule = 'ModuleGridData',                   &
                                     default      =  .false.,                           &
                                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR100'

        DefaultValue = Me%DefaultValue

        call GetData                (Me%DefaultValue, ObjEnterData, flag,               &
                                     keyword      = 'DEFAULT_VALUE',                    &
                                     ClientModule = 'ModuleGridData',                   &
                                     default      =  DefaultValue,                      &
                                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR110'        

                
        select case(Char_TypeZUV)

            case("Z", "z")

                Me%TypeZUV = TypeZ_

            case("U", "u")

                Me%TypeZUV = TypeU_
                Me%WorkSize%JUB        = Me%WorkSize%JUB + 1

            case("V", "v")

                Me%TypeZUV = TypeV_
                Me%WorkSize%IUB        = Me%WorkSize%IUB + 1
           
            case default
                
                write(*,*)'Invalid type ZUV in grid data '//trim(Me%FileName)
                stop 'ReadGridDataFile - ModuleGridData - ERR120'

        end select

        if (Me%Is3D) then
            allocate(Me%GridData3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%KLB:Me%KUB), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR180a'
            
            Me%GridData3D(:,:,:) = Me%DefaultValue        
        else
            allocate(Me%GridData2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%GridData2D(:,:  ) = Me%DefaultValue
        endif        
        

        if (.not.Me%ConstantInSpace) then
            !Looks for data block
            call ReadFromBlocks(ObjEnterData)
        endif
            
        call KillEnterData(ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridDataFile - ModuleGridData - ERR280'

        !----------------------------------------------------------------------

    end subroutine ReadGridDataFile

    !--------------------------------------------------------------------------

    subroutine ReadFromBlocks(ObjEnterData)

        !Arguments-------------------------------------------------------------                                                    
        integer                                     :: ObjEnterData
        !Local-----------------------------------------------------------------
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        integer                                     :: STAT_CALL
        integer                                     :: flag
        integer                                     :: line
        integer                                     :: i, j, k, l, ii, jj
        real, dimension(:), allocatable             :: Aux
        real                                        :: AuxValue
        logical                                     :: Start3DFrom2D 

        !----------------------------------------------------------------------

        !Looks for data block
        if (.not. Me%Is3D) then
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                      &
                                        BeginGridData2D, EndGridData2D, BlockFound,      &
                                        FirstLine = FirstLine, LastLine = LastLine,      &
                                        STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadFromBlocks - ModuleGridData - ERR110'

            if (.not. BlockFound) then

                call SetError (WARNING_, KEYWORD_, 'Block <BeginGridData2D>, <EndGridData2D> not found', OFF)
                call SetError (WARNING_, KEYWORD_, 'Are you using the old format <BeginBathymetry>, <EndBathymetry>?', OFF)
                call SetError (WARNING_, KEYWORD_, 'File : '//trim(adjustl(Me%FileName)), OFF)

                !Tries to read 2D Data (New Format)
                call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                            '<BeginBathymetry>', '<EndBathymetry>', BlockFound,  &
                                            FirstLine = FirstLine, LastLine = LastLine,          &
                                            STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ReadFromBlocks - ModuleGridData - ERR120'
            endif

        else

            Start3DFrom2D = .false.
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                      &
                                        BeginGridData3D, EndGridData3D, BlockFound,      &
                                        FirstLine = FirstLine, LastLine = LastLine,      &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFromBlocks - ModuleGridData - ERR130'
            !Verifies if there is a 2D block to initialize 3D filed (constant in vertical)
            if (.not. BlockFound) then
                call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                      &
                                            BeginGridData2D, EndGridData2D, BlockFound,      &
                                            FirstLine = FirstLine, LastLine = LastLine,      &
                                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFromBlocks - ModuleGridData - ERR140'
                if (BlockFound) Start3DFrom2D = .true.
            endif
        endif
        
BF:     if (BlockFound) then 
            
Is3D:       if (.not. Me%Is3D) then

                if (Me%Evolution%Yes) then
                    allocate(Me%GridData2Dreference(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                    Me%GridData2Dreference(:,:) = Me%DefaultValue
                endif

                line = FirstLine + 1

                allocate  (Aux(3))

                call GetData(Aux, ObjEnterData, flag, Buffer_Line  = Line, STAT = STAT_CALL)

                if (.not.(STAT_CALL == SUCCESS_ .or. STAT_CALL == SIZE_ERR_))            &
                    stop 'ReadFromBlocks - ModuleGridData - ERR150'

Coln1:          if      (flag == 3) then 

                    do  l = Line, LastLine - 1

                        call GetData(Aux, ObjEnterData,                                  &
                                     flag, Buffer_Line  = l, STAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_)                                       &
                            stop 'ReadFromBlocks - ModuleGridData - ERR160'

                        i = int(Aux(1))
                        j = int(Aux(2))

                        if (Me%DDecomp%MasterOrSlave) then
                            if (i>= Me%DDecomp%HaloMap%ILB .and.            &
                                i<= Me%DDecomp%HaloMap%IUB+1) then
                                ii = i + 1 - Me%DDecomp%HaloMap%ILB
                            else
                                cycle
                            endif                                
                        else
                            ii = i
                        endif    
                                
                        if (Me%DDecomp%MasterOrSlave) then
                            if (j>= Me%DDecomp%HaloMap%JLB .and. &
                                j<= Me%DDecomp%HaloMap%JUB+1) then
                                jj = j + 1 - Me%DDecomp%HaloMap%JLB
                            else
                                cycle
                            endif                                
                        else
                            jj = j
                        endif    
                        

                        Me%GridData2D(ii, jj) = Aux(3)

                    enddo

                else if (flag == 1) then Coln1

                    line = FirstLine

                    do i = Me%GlobalWorkSize%ILB, Me%GlobalWorkSize%IUB
                    do j = Me%GlobalWorkSize%JLB, Me%GlobalWorkSize%JUB
                        line = line+1
                
                        !Reached last line before end?
                        if (line .EQ. LastLine) then
                            write(*,*) 'Error in File=', trim(Me%FileName)
                            write(*,*) 'Error in Line=', line
                            write(*,*) 'Error reading GridData2D'
                            stop       'ReadFromBlocks - ModuleGridData - ERR170'
                        end if

                        call GetData(AuxValue, ObjEnterData,  flag,                     &
                                     Buffer_Line  = Line,                               &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFromBlocks - ModuleGridData - ERR180'
                        

                        if (Me%DDecomp%MasterOrSlave) then
                            if (i>= Me%DDecomp%HaloMap%ILB .and. &
                                i<= Me%DDecomp%HaloMap%IUB+1) then
                                ii = i + 1 - Me%DDecomp%HaloMap%ILB
                            else
                                cycle
                            endif                                
                        else
                            ii = i
                        endif    
                                
                        if (Me%DDecomp%MasterOrSlave) then
                            if (j>= Me%DDecomp%HaloMap%JLB .and. &
                                j<= Me%DDecomp%HaloMap%JUB+1) then
                                jj = j + 1 - Me%DDecomp%HaloMap%JLB
                            else
                                cycle
                            endif                                
                        else
                            jj = j
                        endif                            

                        Me%GridData2D(ii,jj) = AuxValue

                    enddo
                    enddo

                endif Coln1

                deallocate(Aux)

            else Is3D

                if (Start3DFrom2D) then

                    line = FirstLine

                    allocate(Aux(1))
                    do i = Me%GlobalWorkSize%ILB, Me%GlobalWorkSize%IUB
                    do j = Me%GlobalWorkSize%JLB, Me%GlobalWorkSize%JUB
                        line = line+1
                
                        !Reached last line before end?
                        if (line .EQ. LastLine) then
                            write(*,*) 'Error in File=', trim(Me%FileName)
                            write(*,*) 'Error in Line=', line
                            write(*,*) 'Error reading GridData2D'
                            stop       'ReadFromBlocks - ModuleGridData - ERR190'
                        end if

                        call GetData(Aux, ObjEnterData,  flag,         & 
                                     Buffer_Line  = Line,              &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFromBlocks - ModuleGridData - ERR200'
                        
                        if (Me%DDecomp%MasterOrSlave) then
                            if (i>= Me%DDecomp%HaloMap%ILB .and. &
                                i<= Me%DDecomp%HaloMap%IUB+1) then
                                ii = i + 1 - Me%DDecomp%HaloMap%ILB
                            else
                                cycle
                            endif                                
                        else
                            ii = i
                        endif    
                                
                        if (Me%DDecomp%MasterOrSlave) then
                            if (j>= Me%DDecomp%HaloMap%JLB .and. &
                                j<= Me%DDecomp%HaloMap%JUB+1) then
                                jj = j + 1 - Me%DDecomp%HaloMap%JLB
                            else
                                cycle
                            endif                                
                        else
                            jj = j
                        endif                            

                        Me%GridData3D(ii, jj,:) = Aux(1)
                        
                    enddo
                    enddo

                else


                    line = FirstLine + 1
                    allocate(Aux(4))


                    call GetData(Aux, ObjEnterData, flag, Buffer_Line  = Line, STAT = STAT_CALL)

                    if (.not.(STAT_CALL == SUCCESS_ .or. STAT_CALL == SIZE_ERR_))            &
                        stop 'ReadFromBlocks - ModuleGridData - ERR210'


Coln:               if (flag == 4)  then
                        do  l = Line, LastLine - 1

                            call GetData(Aux, ObjEnterData, flag, Buffer_Line  = l, STAT = STAT_CALL)

                            if (STAT_CALL /= SUCCESS_)                                       &
                                stop 'ReadFromBlocks - ModuleGridData - ERR220'

                            i = int(Aux(1))
                            j = int(Aux(2))
                            k = int(Aux(3))
                            
                            if (Me%DDecomp%MasterOrSlave) then
                                if (i>= Me%DDecomp%HaloMap%ILB .and. &
                                    i<= Me%DDecomp%HaloMap%IUB+1) then
                                    ii = i + 1 - Me%DDecomp%HaloMap%ILB
                                else
                                    cycle
                                endif                                
                            else
                                ii = i
                            endif    
                                    
                            if (Me%DDecomp%MasterOrSlave) then
                                if (j>= Me%DDecomp%HaloMap%JLB .and. &
                                    j<= Me%DDecomp%HaloMap%JUB+1) then
                                    jj = j + 1 - Me%DDecomp%HaloMap%JLB
                                else
                                    cycle
                                endif                                
                            else
                                jj = j
                            endif                            

                            Me%GridData3D(ii, jj, k) = Aux(4)
                       
                        enddo

                        deallocate(Aux)

                    else if (flag == 3) then Coln

                        deallocate(Aux) 

                        allocate  (Aux(3))

                        do  l = Line, LastLine - 1

                            call GetData(Aux, ObjEnterData,                                  &
                                         flag, Buffer_Line  = l, STAT = STAT_CALL)

                            if (STAT_CALL /= SUCCESS_)                                       &
                                stop 'ReadFromBlocks - ModuleGridData - ERR230'

                            i = int(Aux(1))
                            j = int(Aux(2))
                            
                            if (Me%DDecomp%MasterOrSlave) then
                                if (i>= Me%DDecomp%HaloMap%ILB .and. &
                                    i<= Me%DDecomp%HaloMap%IUB+1) then
                                    ii = i + 1 - Me%DDecomp%HaloMap%ILB
                                else
                                    cycle
                                endif                                
                            else
                                ii = i
                            endif    
                                    
                            if (Me%DDecomp%MasterOrSlave) then
                                if (j>= Me%DDecomp%HaloMap%JLB .and. &
                                    j<= Me%DDecomp%HaloMap%JUB+1) then
                                    jj = j + 1 - Me%DDecomp%HaloMap%JLB
                                else
                                    cycle
                                endif                                
                            else
                                jj = j
                            endif                            

                            !Next line crashes in debug.... replaced by do loop - Frank
                            !Me%GridData3D(i, j, Me%KLB : Me%KUB) = Aux(3)
                            do k = Me%KLB, Me%KUB
                                Me%GridData3D(ii, jj, k) = Aux(3)
                            enddo

                        enddo

                        deallocate(Aux)

                    else if (flag == 1) then Coln


                        line = FirstLine

                        do i = Me%GlobalWorkSize%ILB, Me%GlobalWorkSize%IUB
                        do j = Me%GlobalWorkSize%JLB, Me%GlobalWorkSize%JUB
                        do k = Me%KLB,       Me%KUB
                            line = line+1
                
                            !Reached last line before end?
                            if (line .EQ. LastLine) then
                                write(*,*) 
                                write(*,*) 'Error reading GridData3D'
                                stop       'ReadFromBlocks - ModuleGridData - ERR240'
                            end if

                            call GetData(AuxValue, ObjEnterData,  flag,                 & 
                                         Buffer_Line  = Line,                           &
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReadFromBlocks - ModuleGridData - ERR250'
                            
                            if (Me%DDecomp%MasterOrSlave) then
                                if (i>= Me%DDecomp%HaloMap%ILB .and. &
                                    i<= Me%DDecomp%HaloMap%IUB+1) then
                                    ii = i + 1 - Me%DDecomp%HaloMap%ILB
                                else
                                    cycle
                                endif                                
                            else
                                ii = i
                            endif    
                                    
                            if (Me%DDecomp%MasterOrSlave) then
                                if (j>= Me%DDecomp%HaloMap%JLB .and. &
                                    j<= Me%DDecomp%HaloMap%JUB+1) then
                                    jj = j + 1 - Me%DDecomp%HaloMap%JLB
                                else
                                    cycle
                                endif                                
                            else
                                jj = j
                            endif                            

                            Me%GridData3D(ii,jj,k) = AuxValue                            

                        enddo
                        enddo
                        enddo

                    endif Coln

                endif


            endif Is3D

        else BF
            
            write(*,*)'Invalid Grid Data File'
            write(*,*)'File :',trim(adjustl(Me%FileName))
            stop 'ReadFromBlocks - ModuleGridData - ERR260'

        endif BF

        call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFromBlocks - ModuleGridData - ERR270'


        !----------------------------------------------------------------------

    end subroutine ReadFromBlocks

#ifndef _NO_HDF5    
    subroutine ReadFileEvolution

        !Arguments-------------------------------------------------------------                                                    
                                                                                                     
        !External--------------------------------------------------------------
        real, dimension(:,:  ), pointer         :: Field2D
        integer                                 :: STAT_CALL
        integer                                 :: HDF5_READ
        logical                                 :: exist

        !Local-----------------------------------------------------------------
        integer                                 :: NumberOfInstants, i

        !Begin-----------------------------------------------------------------

        Me%Evolution%OldInstants = 0

        inquire (file=trim(Me%Evolution%File), exist = exist)

ex:     if (exist) then
      
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
        
            call ConstructHDF5 (Me%Evolution%ObjHDF5, trim(Me%Evolution%File), HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFileEvolution - ModuleGridData - ERR01'

            call GetHDF5GroupNumberOfItems(Me%Evolution%ObjHDF5, "/Time",                       &
                                           NumberOfInstants, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFileEvolution - ModuleGridData - ERR02'

            call GetComputeCurrentTime(Me%Evolution%ObjTime, Me%Evolution%Now, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFileEvolution - ModuleGridData - ERR03'

            do i = 1, NumberOfInstants

                if (HDF5TimeInstant(i) > Me%Evolution%Now) exit
            
            enddo

            Me%Evolution%OldInstants = i - 1

old:        if (Me%Evolution%OldInstants > 0) then

                allocate (Me%Evolution%TimeInstants(6,Me%Evolution%OldInstants))

                allocate (Me%Evolution%GridData2D(Me%Size%ILB:Me%Size%IUB,               &
                                                  Me%Size%JLB:Me%Size%JUB,               &
                                                  Me%Evolution%OldInstants))

                Me%Evolution%GridData2D(:,:,:) = FillValueReal

                allocate (Field2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB)) 

                Field2D (:,:) = FillValueReal
            
                do i = 1, Me%Evolution%OldInstants

                   call ExtractDate(HDF5TimeInstant(i), Year   = Me%Evolution%TimeInstants(1,i),    &
                                                        Month  = Me%Evolution%TimeInstants(2,i),    &
                                                        Day    = Me%Evolution%TimeInstants(3,i),    &
                                                        Hour   = Me%Evolution%TimeInstants(4,i),    & 
                                                        Minute = Me%Evolution%TimeInstants(5,i),    &
                                                        Second = Me%Evolution%TimeInstants(6,i))

                    
            
                    call ReadHDF5Values2D (i, Field2D)

                    !The bathymetry actualize
                    if (HDF5TimeInstant(i) == Me%Evolution%Now) Me%GridData2D(:,:) = Field2D(:,:)

                enddo

                deallocate (Field2D) 

                nullify(Field2D)


            endif old



            call KillHDF5(Me%Evolution%ObjHDF5, STAT = STAT_CALL)
        
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFileEvolution - ModuleGridData - ERR04'

        endif ex



    end subroutine ReadFileEvolution
#endif
        !----------------------------------------------------------------------



    !--------------------------------------------------------------------------
#ifndef _NO_HDF5
    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%Evolution%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%Evolution%ObjHDF5,                     &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleGridData - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    
    !--------------------------------------------------------------------------


    subroutine ReadHDF5Values2D (Instant, Field)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        real, dimension(:,:), pointer           :: Field
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%Evolution%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,           &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleGridData - ERR01'
        
        
        call HDF5ReadData(Me%Evolution%ObjHDF5, "/Results/"//trim(Me%Evolution%PropName),    &
                          trim(Me%Evolution%PropName),                                       &
                          Array2D = Field, OutputNumber = Instant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleGridData - ERR02'


    end subroutine ReadHDF5Values2D
#endif
   
    
    !--------------------------------------------------------------------------

    

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetGridData2D(GridDataID, GridData2D, STAT)     

        !Arguments---------------------------------------------------------------
        integer                                     :: GridDataID
        real, pointer, dimension(:,:)               :: GridData2D
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (.not. associated(Me%GridData2D)) then
                stop 'GetGridData2D - ModuleGridData - ERR01'
            endif

            call Read_Lock(mGridData_, Me%InstanceID)
            GridData2D => Me%GridData2D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetGridData2D

    !--------------------------------------------------------------------------

    subroutine GetGridData3D(GridDataID, GridData3D, STAT)     

        !Arguments---------------------------------------------------------------
        integer                                     :: GridDataID
        real, pointer, dimension(:,:,:)             :: GridData3D
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (.not. associated(Me%GridData3D)) then
                stop 'GetGridData3D - ModuleGridData - ERR01'
            endif

            call Read_Lock(mGridData_, Me%InstanceID)
            GridData3D => Me%GridData3D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridData3D

    !--------------------------------------------------------------------------
    subroutine GetGridData2DReference(GridDataID, GridData2Dreference, STAT)     

        !Arguments---------------------------------------------------------------
        integer                                     :: GridDataID
        real, pointer, dimension(:,:)               :: GridData2Dreference
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (.not. associated(Me%GridData2Dreference)) then
                stop 'GetGridData2Dreference - ModuleGridData - ERR01'
            endif

            call Read_Lock(mGridData_, Me%InstanceID)
            GridData2Dreference => Me%GridData2Dreference

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetGridData2Dreference


    !--------------------------------------------------------------------------

    subroutine GetGridDataFileName(GridDataID, FileName, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        character(LEN = *), intent(OUT)             :: FileName
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            FileName = Me%FileName

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridDataFileName

    !--------------------------------------------------------------------------

    subroutine GetMaximumValue(GridDataID, MaximumValue, STAT)     

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        real,              intent(OUT)              :: MaximumValue
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            MaximumValue = Me%MaximumValue

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetMaximumValue

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine GetMaxValueInPolygon(GridDataID, MaximumValue, Polygon, STAT)     

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        real,              intent(OUT)              :: MaximumValue
        type(T_Polygon),                 pointer    :: Polygon
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        real,  dimension(:,:), pointer              :: CoordX, CoordY
        type (T_PointF),   pointer                  :: Point
        integer                                     :: ready_, i, j, STAT_CALL        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call GetZCoordinates(Me%ObjHorizontalGrid, CoordX, CoordY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetMaxValueInPolygon - ModuleGridData - ERR10'

            allocate(Point)
            
            MaximumValue = FillValueReal
            
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
                Point%X = CoordX(i, j)
                Point%Y = CoordY(i, j)                    
                
                if (Me%GridData2D(i, j) .GT. MaximumValue .and. IsPointInsidePolygon(Point, Polygon)) &
                    MaximumValue = Me%GridData2D(i, j)
                    
            end do       
            end do     
            
            deallocate(Point)
            
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetMaxValueInPolygon - ModuleGridData - ERR20'
                        
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetMaxValueInPolygon - ModuleGridData - ERR30'

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetMaxValueInPolygon

    !--------------------------------------------------------------------------    

    subroutine GetMinimumValue(GridDataID, MinimumValue, STAT)     

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        real,              intent(OUT)              :: MinimumValue
        integer, optional, intent(OUT)              :: STAT
             

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            MinimumValue = Me%MinimumValue

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetMinimumValue

    !--------------------------------------------------------------------------
    
    subroutine GetGridDataType(GridDataID, TypeZUV, STAT)     

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        integer,           intent(OUT)              :: TypeZUV
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            TypeZUV = Me%TypeZUV

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridDataType
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine GetIsGridData3D(GridDataID, Is3D, STAT)     

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        logical,           intent(OUT)              :: Is3D
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Is3D = Me%Is3D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetIsGridData3D
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine GetFillValue(GridDataID, FillValue, STAT)     

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        real,              intent(OUT)              :: FillValue
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            FillValue = Me%FillValue

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetFillValue

    !--------------------------------------------------------------------------

    
    subroutine GetGridDataEvolution(GridDataID, Evolution, STAT)     

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        logical,           intent(OUT)              :: Evolution
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Evolution = Me%Evolution%Yes

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridDataEvolution

    !--------------------------------------------------------------------------


    subroutine UngetGridData2D(GridDataID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        real, pointer, dimension(:,:)               :: Array
        integer, optional, intent (OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            nullify(Array)

            call Read_UnLock(mGridData_, Me%InstanceID, "UngetGridData")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetGridData2D

    !--------------------------------------------------------------------------

    subroutine UngetGridData3D(GridDataID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        real, pointer, dimension(:,:,:)             :: Array
        integer, optional, intent (OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            nullify(Array)

            call Read_UnLock(mGridData_, Me%InstanceID, "UngetGridData")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetGridData3D


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyGridData2DIncrement(GridDataID, Increment2D, Add, STAT)     

        !Arguments---------------------------------------------------------------
        integer                                     :: GridDataID
        real, pointer, dimension(:,:)               :: Increment2D
        logical,    intent(IN)                      :: Add
        integer, optional, intent(OUT)              :: STAT
        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: i, j, STAT_
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (Me%Evolution%Yes) then

                if (MonitorPerformance) then
                    call StartWatch ("ModuleGridData", "ModifyGridData2DIncrement")
                endif

                CHUNK = CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
                if (Add) then
                    !$OMP PARALLEL PRIVATE(i,j)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j=Me%WorkSize%JLB,Me%WorkSize%JUB
                    do i=Me%WorkSize%ILB,Me%WorkSize%IUB
                        Me%GridData2D(i, j) =  Me%GridData2D(i, j) + Increment2D(i, j)
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL
                else
                    !$OMP PARALLEL PRIVATE(i,j)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j=Me%WorkSize%JLB,Me%WorkSize%JUB
                    do i=Me%WorkSize%ILB,Me%WorkSize%IUB
                        Me%GridData2D(i, j) =  Me%GridData2D(i, j) - Increment2D(i, j)
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL
                endif

                if (MonitorPerformance) then
                    call StopWatch ("ModuleGridData", "ModifyGridData2DIncrement")
                endif

                STAT_ = SUCCESS_

            else

                STAT_ = UNKNOWN_

            endif

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine ModifyGridData2DIncrement

    !----------------------------------------------------------------------------

    subroutine ModifyNewMatrixGridData2D(GridDataID, NewGridData2D, STAT)     

        !Arguments---------------------------------------------------------------
        integer                                     :: GridDataID
        real, pointer, dimension(:,:)               :: NewGridData2D
        integer, optional, intent(OUT)              :: STAT
        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: i, j, STAT_
        !T integer                                     :: CHUNK

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (MonitorPerformance) then
                call StartWatch ("ModuleGridData", "ModifyNewMatrixGridData2D")
            endif

            !ACanas: Paralelization not tested because subrotine not used in MOHID Water
            
            !T CHUNK = CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
            !T !$OMP PARALLEL PRIVATE(i,j)
            !T !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j=Me%WorkSize%JLB,Me%WorkSize%JUB
            do i=Me%WorkSize%ILB,Me%WorkSize%IUB
                Me%GridData2D(i, j) =  NewGridData2D(i, j)
            enddo
            enddo
            !T !$OMP END DO
            !T !$OMP END PARALLEL
            
            if (MonitorPerformance) then
                call StopWatch ("ModuleGridData", "ModifyNewMatrixGridData2D")
            endif 
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine ModifyNewMatrixGridData2D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine ModifyConstantGridData2D(GridDataID, ConstantValue, STAT)     

        !Arguments---------------------------------------------------------------
        integer                                     :: GridDataID
        real                                        :: ConstantValue
        integer, optional, intent(OUT)              :: STAT
        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: i, j, STAT_
        !T integer                                     :: CHUNK

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (MonitorPerformance) then
                call StartWatch ("ModuleGridData", "ModifyConstantGridData2D")
            endif

            !ACanas: Paralelization not tested because subrotine not used in MOHID Water

            !T CHUNK = CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
            !T !$OMP PARALLEL PRIVATE(i,j)
            !T !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j=Me%WorkSize%JLB,Me%WorkSize%JUB
            do i=Me%WorkSize%ILB,Me%WorkSize%IUB
                Me%GridData2D(i, j) =  ConstantValue
            enddo
            enddo
            !T !$OMP END DO
            !T !$OMP END PARALLEL
            
            if (MonitorPerformance) then
                call StopWatch ("ModuleGridData", "ModifyConstantGridData2D")
            endif
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine ModifyConstantGridData2D

    !----------------------------------------------------------------------------

    subroutine WriteGridData_v1 (FileName,                                          &
                                 COMENT1, COMENT2,                                  &
                                 HorizontalGridID,                                  &
                                 FillValue, Overwrite,                              &
                                 GridData2D_Real, GridData2D_Int,                   &
                                 GridData3D_Real, GridData3D_Int,                   &
                                 KLB, KUB,                                          &
                                 STAT)                                                                           
        
        !Arguments---------------------------------------------------------------
        character(LEN = *), intent(IN)              :: FileName
        character(LEN = *), intent(IN)              :: COMENT1          
        character(LEN = *), intent(IN)              :: COMENT2   
        integer                                     :: HorizontalGridID       
        real                                        :: FillValue
        logical, intent(IN)                         :: Overwrite
        real,    pointer, optional, dimension(:,:)  :: GridData2D_Real
        integer, pointer, optional, dimension(:,:)  :: GridData2D_Int
        real,    pointer, optional, dimension(:,:,:):: GridData3D_Real
        integer, pointer, optional, dimension(:,:,:):: GridData3D_Int
        integer, optional, intent(IN)               :: KLB, KUB
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        real                                        :: Xorig      
        real                                        :: Yorig       
        real                                        :: GRID_ANGLE  
        integer                                     :: Zone        
        real                                        :: Longitude
        real                                        :: Latitude    
        integer                                     :: CoordType  
        real, pointer,           dimension(:  )     :: XX, YY
        real, pointer,           dimension(:,:)     :: XX_IE, YY_IE
        real, pointer,           dimension(:,:)     :: LatConn, LongConn
        type (T_Size2D)                             :: WorkSize

        integer                                     :: Unit
        integer                                     :: STAT_CALL
        integer                                     :: UTM_, MIL_PORT_
        integer                                     :: PORTUGUESE_UTM_ZONE_, SIMPLE_GEOG_
        logical                                     :: exist, DistortionYes
        integer                                     :: i, j, k

        !------------------------------------------------------------------------

        !Gets Origin
        call GetGridOrigin(HorizontalGridID, Xorig, Yorig, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR10'
        
        !Gets GridAngle
        call GetGridAngle(HorizontalGridID, GRID_ANGLE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR20'

        !Gets GridZone
        call GetGridZone (HorizontalGridID, Zone, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR30'

        !Gets Lat/lon
        call GetLatitudeLongitude (HorizontalGridID, Latitude, Longitude, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR40'
        
        !Gets Lat/lon
        call GetGridCoordType (HorizontalGridID, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR50'
        
        !Gets XX/ YY and XX_IE, YY_IE
        call GetHorizontalGrid(HorizontalGridID, XX = XX, YY = YY, XX_IE = XX_IE, YY_IE = YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR60'

        call GetGridLatitudeLongitude(HorizontalGridID,                     &
                                      GridLatitudeConn  = LatConn,          &
                                      GridLongitudeConn = LongConn,         &
                                      STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR70'


        !Gets WorkSize
        call GetHorizontalGridSize(HorizontalGridID, WorkSize = WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR80'

        call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR90'

        if(.not. Overwrite)then
            !Verifies if file exists
            inquire (file = FileName, exist = exist)
            if (exist) then
                write(*,* )'Cannot write Grid Data. File Exists already!'
                write(*, *)trim(adjustl(FileName))
                stop 'WriteGridData_v1 - ModuleGridData - ERR100' 
            endif

            !Opens file
            open (Unit, file = FileName, status = 'new', IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR110'
        else
            !Opens file
            open (Unit, file = FileName, status = 'replace', IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) then 
                write(*,*) 'Was not possible to write to file ', FileName
                write(*,*) 'OPEN returned IOSTAT = ', STAT_CALL
                stop 'WriteGridData_v1 - ModuleGridData - ERR120'
            endif
        end if

        call WriteDataLine(unit, 'COMENT1', COMENT1)
        call WriteDataLine(unit, 'COMENT2', COMENT2)
        write(unit, *)                                                                                                      
        write(unit, *)                                                                                                      
                       
        write(unit, *) 'ILB_IUB    :', WorkSize%ILB, WorkSize%IUB
        write(unit, *) 'JLB_JUB    :', WorkSize%JLB, WorkSize%JUB
        
        write(unit, *) 'COORD_TIP  :', CoordType

        call GetCoordTypeList(UTM = UTM_, MIL_PORT = MIL_PORT_, &
                              PORTUGUESE_UTM_ZONE = PORTUGUESE_UTM_ZONE_, &
                              SIMPLE_GEOG = SIMPLE_GEOG_)
        if ((CoordType .EQ. UTM_) .OR. (CoordType .EQ. MIL_PORT_)) then
            if (CoordType .EQ. MIL_PORT_) then
                write(unit, *)  'ZONE       :', PORTUGUESE_UTM_ZONE_
            else
                write(unit, *)  'ZONE       :', Zone
            endif
        end if

        write(unit, *) 'ORIGIN     :', Xorig, Yorig
        write(unit, *) 'GRID_ANGLE :', GRID_ANGLE
        write(unit, *) 'LATITUDE   :', Latitude           
        write(unit, *) 'LONGITUDE  :', Longitude             
        write(unit, *) 'FILL_VALUE :', FillValue
        write(unit, *) 
        write(unit, *) 

        call GetCheckDistortion(HorizontalGridID, DistortionYes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR130'

        if (DistortionYes) then

            write(unit, *) "<CornersXY>"

            if (CoordType .EQ. SIMPLE_GEOG_) then
                
                do i = WorkSize%ILB, WorkSize%IUB+1
                do j = WorkSize%JLB, WorkSize%JUB+1 
                    write(unit, *) LongConn(i,j), LatConn(i,j)
                end do
                end do

            else

                do i = WorkSize%ILB, WorkSize%IUB+1
                do j = WorkSize%JLB, WorkSize%JUB+1 
                    write(unit, *) XX_IE(i,j),YY_IE(i,j)
                end do
                end do

            end if

            write(unit, *) "<"//backslash//"CornersXY>"             
!            write(unit, *) "<\CornersXY>"


        else

            write(unit, *) "<BeginXX>"
            do j = WorkSize%JLB, WorkSize%JUB+1 
                write(unit, *) XX(j)                       
            end do
            write(unit, *) "<EndXX>"


            write(unit, *) "<BeginYY>"
            do i = WorkSize%ILB, WorkSize%IUB+1
                write(unit, *) YY(i)  
            end do
            write(unit, *) "<EndYY>"

        endif

        if ((present(GridData3D_Real) .or. present(GridData3D_Int)) .and.               &
            (.not. present(KLB) .or. .not. present(KUB))) then
            if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR135'
        endif
   
        if (present(GridData2D_Real)) then
            write(unit, *)trim(BeginGridData2D)   
            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
                write(unit, *) GridData2D_Real(i, j)
            end do
            end do
            write(unit, *)trim(EndGridData2D)
        elseif(present(GridData2D_Int)) then
            write(unit, *)trim(BeginGridData2D)   
            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
                write(unit, *) GridData2D_Int(i, j)
            end do
            end do
            write(unit, *)trim(EndGridData2D)
        elseif(present(GridData3D_Int)) then
            write(unit, *)trim(BeginGridData3D)   
            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
            do k = KLB         , KUB   
                write(unit, *) GridData3D_Int(i, j, k)
            enddo
            end do
            end do
            write(unit, *)trim(EndGridData3D)
        elseif(present(GridData3D_Real)) then
            write(unit, *)trim(BeginGridData3D)   
            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
            do k = KLB         , KUB   
                write(unit, *) GridData3D_Real(i, j, k)
            enddo
            end do
            end do
            write(unit, *)trim(EndGridData3D)
        else
            write(unit, *)trim(BeginGridData2D)   
            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
                write(unit, *) Me%GridData2D(i, j)
            end do
            end do
            write(unit, *)trim(EndGridData2D)
        end if

        !Closes Files             
        call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR140'

        !Ungets XX, YY
        call UngetHorizontalGrid (HorizontalGridID, XX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR150'

        call UngetHorizontalGrid (HorizontalGridID, YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR160'

        call UngetHorizontalGrid(HorizontalGridID, XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR170'

        call UngetHorizontalGrid(HorizontalGridID, YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR180'

        call UngetHorizontalGrid(HorizontalGridID, LatConn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR190'

        call UngetHorizontalGrid(HorizontalGridID, LongConn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v1 - ModuleGridData - ERR200'



        if (present(STAT)) STAT = SUCCESS_

        !------------------------------------------------------------------------

    end subroutine WriteGridData_v1

    !----------------------------------------------------------------------------

    subroutine WriteGridData_v2  (FileName,                                         &
                                  XX, YY,                                           &
                                  ConnectionX, ConnectionY,                         &
                                  COMENT1, COMENT2,                                 &
                                  WorkSize,                                         &  
                                  CoordType,                                        &
                                  Xorig, Yorig,                                     &
                                  Zone,GRID_ANGLE,                                  &
                                  Latitude, Longitude,                              &
                                  FillValue, Overwrite,                             &
                                  GridData2D_Real, GridData2D_Int,                  &
                                  Datum, ProjType, SP1, SP2,                                  &
                                  STAT)                                                                           
        
        !Arguments---------------------------------------------------------------
        character(LEN = *), intent(IN)              :: FileName
        character(LEN = *), intent(IN)              :: COMENT1          
        character(LEN = *), intent(IN)              :: COMENT2          
        real, pointer, optional,     dimension(:  ) :: XX
        real, pointer, optional,     dimension(:  ) :: YY
        real, pointer, optional,     dimension(:,:) :: ConnectionX
        real, pointer, optional,     dimension(:,:) :: ConnectionY
        real, intent(IN)                            :: Xorig      
        real, intent(IN)                            :: Yorig       
        real, intent(IN)                            :: Latitude    
        real, intent(IN)                            :: Longitude
        real, intent(IN)                            :: FillValue
        real, intent(IN)                            :: GRID_ANGLE  
        integer, intent(IN)                         :: CoordType  
        integer, intent(IN)                         :: Zone        
        type (T_Size2D), intent(IN)                 :: WorkSize
        real,    pointer, optional, dimension(:,:)  :: GridData2D_Real
        integer, pointer, optional, dimension(:,:)  :: GridData2D_Int
        logical, intent(IN)                         :: Overwrite
        integer, optional, intent(IN)               :: Datum
        integer, optional, intent(IN)               :: ProjType
        real   , optional, intent(IN)               :: SP1, SP2
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: Unit
        integer                                     :: STAT_CALL
        integer                                     :: UTM_, MIL_PORT_, PORTUGUESE_UTM_ZONE_
        logical                                     :: exist
        integer                                     :: i, j

        !------------------------------------------------------------------------

        call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v2 - ModuleGridData - ERR10'

        if(.not. Overwrite)then
            !Verifies if file exists
            inquire (file = FileName, exist = exist)
            if (exist) then
                write(*,* )'Cannot write Grid Data. File Exists already!'
                write(*, *)trim(adjustl(FileName))
                stop 'WriteGridData_v2 - ModuleGridData - ERR20' 
            endif

            !Opens file
            open (Unit, file = FileName, status = 'new', IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v2 - ModuleGridData - ERR30'
        else
            !Opens file
            open (Unit, file = FileName, status = 'replace', IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) then
                write(*,*)"Error trying to open the output file: "//trim(FileName)
                write(*,*)"Please check if the path is correct."
                stop 'WriteGridData_v2 - ModuleGridData - ERR40'
            end if
        end if

        call WriteDataLine(unit, 'COMENT1', COMENT1)
        call WriteDataLine(unit, 'COMENT2', COMENT2)

        write(unit, *)                                                                                                      
        write(unit, *)                                                                                                      
                       
        write(unit, *) 'ILB_IUB                :', WorkSize%ILB, WorkSize%IUB
        write(unit, *) 'JLB_JUB                :', WorkSize%JLB, WorkSize%JUB
        
        write(unit, *) 'COORD_TIP              :', CoordType

        call GetCoordTypeList(UTM = UTM_, MIL_PORT = MIL_PORT_, &
                              PORTUGUESE_UTM_ZONE = PORTUGUESE_UTM_ZONE_)
        if ((CoordType .EQ. UTM_) .OR. (CoordType .EQ. MIL_PORT_)) then
            if (CoordType .EQ. MIL_PORT_) then
                write(unit, *)  'ZONE                  :', PORTUGUESE_UTM_ZONE_
            else
                write(unit, *)  'ZONE                  :', Zone
            endif
        end if

        if (present(Datum))     write(unit, *) 'DATUM                  :', Datum
        if (present(ProjType))  write(unit, *) 'PROJ_TYPE              :', ProjType
        if (present(SP1  ))     write(unit, *) 'SP1                    :', SP1
        if (present(SP2  ))     write(unit, *) 'SP2                    :', SP2

       
        write(unit, *) 'ORIGIN                 :', Xorig, Yorig
        write(unit, *) 'GRID_ANGLE             :', GRID_ANGLE
        write(unit, *) 'LATITUDE               :', Latitude           
        write(unit, *) 'LONGITUDE              :', Longitude             
        write(unit, *) 'FILL_VALUE             :', FillValue
        write(unit, *) 
        write(unit, *)

        if(present(XX))then

            write(unit, *) "<BeginXX>"
            do j = WorkSize%JLB, WorkSize%JUB+1 
                write(unit, *) XX(j)                       
            end do
            write(unit, *) "<EndXX>"

        end if

        if(present(XX))then

            write(unit, *) "<BeginYY>"
            do i = WorkSize%ILB, WorkSize%IUB+1
                write(unit, *) YY(i)  
            end do
            write(unit, *) "<EndYY>"

        end if

        if(present(ConnectionX) .and. present(ConnectionY))then
            write(unit, *) "<CornersXY>"
            do i = WorkSize%ILB, WorkSize%IUB+1
            do j = WorkSize%JLB, WorkSize%JUB+1 
                write(unit, *) ConnectionX(i,j), ConnectionY(i,j)
            end do
            end do
            write(unit, *) "<"//backslash//"CornersXY>"
!            write(unit, *) "<\CornersXY>"
            write(unit, *)
        end if

        if (present(GridData2D_Real)) then
            write(unit, *)trim(BeginGridData2D)   
            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
                write(unit, *) GridData2D_Real(i, j)
            end do
            end do
            write(unit, *)trim(EndGridData2D)
        elseif(present(GridData2D_Int)) then
            write(unit, *)trim(BeginGridData2D)   
            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
                write(unit, *) GridData2D_Int(i, j)
            end do
            end do
            write(unit, *)trim(EndGridData2D)
        else
            write(unit, *)trim(BeginGridData2D)   
            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
                write(unit, *) Me%GridData2D(i, j)
            end do
            end do
            write(unit, *)trim(EndGridData2D)
        end if

        !Closes Files             
        call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridData_v2 - ModuleGridData - ERR50'

        if (present(STAT)) STAT = SUCCESS_

        !------------------------------------------------------------------------

    end subroutine WriteGridData_v2

#ifndef _NO_HDF5
    subroutine WriteFileEvolution

        !Arguments-------------------------------------------------------------                                                    
                                                                                                     
        !Local-----------------------------------------------------------------
        real,    dimension(:,:  ), pointer      :: Array2D
        integer, dimension(:,:  ), pointer      :: MappingPoints2D
        real,    dimension(:    ), pointer      :: TimePtr
        real,    dimension(6    ), target       :: AuxTime
        integer                                 :: STAT_CALL
        integer                                 :: HDF5_CREATE
        integer                                 :: ILB, IUB, JLB, JUB, i, j
        integer                                 :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        !----------------------------------------------------------------------


        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 

        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 


        !Begin-----------------------------------------------------------------

      
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
    
        call ConstructHDF5 (Me%Evolution%ObjHDF5, trim(Me%Evolution%File), HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR01'

        call GetComputeCurrentTime(Me%Evolution%ObjTime, Me%Evolution%Now, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR02'


        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%Evolution%ObjHDF5,             &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR03'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%Evolution%ObjHDF5, WorkILB, WorkIUB, WorkJLB,           &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR04'


        !Writes the Grid
        call HDF5WriteData   (Me%Evolution%ObjHDF5, "/Grid",                             &
                              trim(Me%Evolution%PropName),"m",                           &
                              Array2D = Me%GridData2Dreference,                          &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR05'

        allocate(MappingPoints2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        do i = Me%WorkSize%ILB, Me%WorkSize%IUB  
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB  

            if (Me%GridData2Dreference(i,j) < -50.) then
            
                MappingPoints2D(i, j) = 0

            else

                MappingPoints2D(i, j) = 1

            endif

        enddo
        enddo

        call HDF5WriteData   (Me%Evolution%ObjHDF5, "/Grid", "MappingPoints2D", "-", &
                              Array2D = MappingPoints2D,                             &
                              STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR06'

        deallocate(MappingPoints2D)

d1:     do i = 1, Me%Evolution%OldInstants + 1

        
            if ( i <= Me%Evolution%OldInstants ) then

                TimePtr => Me%Evolution%TimeInstants(:,i)
                
                Array2D => Me%Evolution%GridData2D(:,:,i)

            else

                !Writes current time
                call ExtractDate   (Me%Evolution%Now, AuxTime(1), AuxTime(2), AuxTime(3),&
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime

                Array2D => Me%GridData2D

            endif

            call HDF5SetLimits  (Me%Evolution%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR07'

            call HDF5WriteData  (Me%Evolution%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS", &
                                 Array1D = TimePtr, OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR08'

            !Writes GridData instant i
            call HDF5SetLimits  (Me%Evolution%ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                                 WorkJUB, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR09'

            call HDF5WriteData  (Me%Evolution%ObjHDF5, "/Results/"//trim(Me%Evolution%PropName),   &
                                trim(Me%Evolution%PropName),                             &
                                 "-", Array2D = Array2D,                                 &
                                OutputNumber = i, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR10'

        enddo d1

        if (Me%Evolution%OldInstants > 0) then

            deallocate(Me%Evolution%TimeInstants)
            
            deallocate(Me%Evolution%GridData2D)

        endif

        call KillHDF5(Me%Evolution%ObjHDF5, STAT = STAT_CALL)
        
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteFileEvolution - ModuleGridData - ERR11'

    end subroutine WriteFileEvolution
#endif
        !----------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillGridData(GridDataID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GridDataID
        integer, optional, intent(OUT)              :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GridDataID, ready_)    
        
cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mGridData_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%Evolution%Yes) then
                    call WriteFileEvolution
                endif

                if (Me%Evolution%ObjTime /= 0) nUsers = DeassociateInstance (mTIME_, Me%Evolution%ObjTime)


                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillGridData - ModuleGridData - ERR02'
                   
                if (associated(Me%GridData2D)) then
                    deallocate(Me%GridData2D)
                    nullify   (Me%GridData2D)
                    if (Me%Evolution%Yes) then
                        deallocate(Me%GridData2Dreference)
                        nullify   (Me%GridData2Dreference)
                    endif
                endif

                if (associated(Me%GridData3D)) then
                    deallocate(Me%GridData3D)
                    nullify   (Me%GridData3D)
                endif

                call DeallocateInstance
                GridDataID = 0
                STAT_        = SUCCESS_

            end if 
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillGridData

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_GridData), pointer          :: AuxObjGridData
        type (T_GridData), pointer          :: PreviousObjGridData

        !Updates pointers
        if (Me%InstanceID == FirstGridData%InstanceID) then
            FirstGridData => FirstGridData%Next
        else
            PreviousObjGridData => FirstGridData
            AuxObjGridData      => FirstGridData%Next
            do while (AuxObjGridData%InstanceID /= Me%InstanceID)
                PreviousObjGridData => AuxObjGridData
                AuxObjGridData      => AuxObjGridData%Next
            enddo

            !Now update linked list
            PreviousObjGridData%Next => AuxObjGridData%Next

        endif
            
        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

    end subroutine DeallocateInstance


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Ready (ObjGridData_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGridData_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjGridData_ID > 0) then
            call LocateObjGridData(ObjGridData_ID)
            ready_ = VerifyReadLock (mGridData_,  Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjGridData (ObjGridDataID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGridDataID

        !Local-----------------------------------------------------------------

        Me => FirstGridData
        do while (associated (Me))
            if (Me%InstanceID == ObjGridDataID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleGridData - LocateObjGridData - ERR01'

    end subroutine LocateObjGridData

    !--------------------------------------------------------------------------


end module ModuleGridData

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
