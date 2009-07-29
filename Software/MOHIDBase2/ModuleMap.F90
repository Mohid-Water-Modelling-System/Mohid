!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : ModuleMap 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to which defines grid points into several categories (see below)
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
! This module "maps" the following points/faces:
!
! WaterPoints3D         - All points where:
!                           - WaterPoints2D == 1 AND k >= kfloorZ
!
! LandPoints3D          - All points where:
!                           - LandPoints2D  == 1
!                           - WaterPoints2D == 1 AND k < kfloorZ
!
! OpenPoints3D          - All points where 
!                           - ComputeFaces3D%U(i, j, k) + ComputeFaces3D%U(i, j+1, k) + 
!                             ComputeFaces3D%V(i, j, k) + ComputeFaces3D%V(i+1, j, k) + 
!                             ComputeFaces3D%W(i, j, k) + ComputeFaces3D%W(i, j, k+1) > 0 
!
! ComputeFaces3D%U      - All faces where:
!                           - ComputeFace2D%U == 1 AND k >= kfloorU
!
! ComputeFaces3D%V      - All faces where:
!                           - ComputeFace2D%V == 1 AND k >= kfloorV
!
! ComputeFaces3D%W      - All faces in (k+1) where
!                           - ComputeFaces3D%U(i, j, k) + ComputeFaces3D%U(i, j+1, k) + 
!                             ComputeFaces3D%V(i, j, k) + ComputeFaces3D%V(i+1, j, k) > 0 
!                           - expect where BoudaryPoints2D(i, j) == 1
!
! LandBoundaryFaces3D   - All faces which have:
!                           - one side of the face a WaterPoints3D
!                           - and one the other side a LandPoints3D
! ImposedNormalFaces    - All faces which are:
!                           - use to defined the boundary conditions normal to the boundary
!
! ImposedTangentialFaces- All faces which are:
!                           - use to defined the boundary conditions tangential to the boundary
!
! WetFaces%U          - All faces which have (in U direction):
!                           - at least one side of the face a OpenPoints3D = 1
!
! WetFaces%V          - All faces which have (in V direction):
!                           - at least one side of the face a OpenPoints3D = 1
!
! History
!   - Frank     Jun 99     : Creation
!   - Frank     Jul 99     : LandBoundaryFaces3D
!   - Frank     Sep 00     : UnCoveredFaces removed
!
!
module ModuleMap

    use ModuleGlobalData
    use ModuleTime                 
    use ModuleGridData,         only: GetGridData, UngetGridData,               & 
                                      GetGridDataFileName, WriteGridData          
    use ModuleHorizontalMap,    only: GetWaterPoints2D, GetLandPoints2D,        &
                                      UnGetHorizontalMap, UpdateComputeFaces2D, &
                                      GetComputeFaces2D, GetBoundaries,         &
                                      GetExteriorBoundaryFaces, GetBoundaryFaces     
    use ModuleGeometry,         only: GetGeometrySize, GetGeometryKFloor,       &
                                      UnGetGeometry, GetGeometryMinWaterColumn
    use ModuleHorizontalGrid,   only: GetHorizontalGrid, GetGridCoordType,      &
                                      GetGridAngle, GetLatitudeLongitude,       &
                                      GetGridZone, GetGridOrigin,               & 
                                      GetCheckDistortion
    use ModuleFunctions,        only : SetMatrixValue, Chunk_J
    use ModuleStopWatch,        only : StartWatch, StopWatch         


    implicit none
    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructMap
    private ::      AllocateInstance
    private ::      AllocateVariables
    private ::      ConstructLandBoundaryFaces3D
    private ::      CheckIsolatedCell
    private ::      CorrectIsolatedCells


    !Selector
    public  :: GetLandPoints3D
    public  :: GetWaterPoints3D
    public  :: GetOpenPoints3D
    public  :: GetComputeFaces3D
    public  :: GetLandBoundaryFaces3D
    public  :: GetImposedNormalFaces
    public  :: GetImposedTangentialFaces
    public  :: GetWetFaces


    public  :: UngetMap



    !Modifier
    public  :: UpdateComputeFaces3D
    private ::      UpdateOpenPoints3D
    private ::      UpdateImposedValues

    public  :: SetComputesFaces3D

    !Destructor
    public  :: KillMap
    private ::      DeallocateInstance



    !Management
    private ::      Ready
    private ::          LocateObjMap

    !Interfaces----------------------------------------------------------------

    private :: UpdateSoilComputeFaces3D
    private :: UpdateMovBoundCompFaces3D

    interface  UpdateComputeFaces3D
        module procedure UpdateSoilComputeFaces3D
        module procedure UpdateMovBoundCompFaces3D
        module procedure UpdateSedimentCompFaces3D
        module procedure UpdateUnsaturatedComputeFaces3D
    end interface  UpdateComputeFaces3D

    private :: UngetMap3D
    interface  UngetMap
        module procedure UngetMap3D
    end interface  UngetMap


    !Types---------------------------------------------------------------------

    private :: T_3D_INT
    type       T_3D_INT
        integer, pointer, dimension(:,:,:)  :: U, V, W
    end type T_3D_INT


    private :: T_2D_INT
    type       T_2D_INT
        integer, pointer, dimension(:,:,:)  :: U, V
    end type T_2D_INT

    private :: T_MapCell
    type       T_MapCell
        integer                             :: Icell
        integer                             :: Jcell
        type(T_MapCell),            pointer :: Next                 => null()
    end type T_MapCell

    private :: T_Map
    type       T_Map
        private 
        integer                             :: InstanceID
        type(T_Size3D  )                    :: Size
        type(T_Size3D  )                    :: WorkSize
        type(T_Time    )                    :: ActualTime
        type(T_3D_INT  )                    :: ComputeFaces3D 
        type(T_3D_INT  )                    :: LandBoundaryFaces3D
        type(T_2D_INT  )                    :: ImposedNormalFaces
        type(T_2D_INT  )                    :: ImposedTangentialFaces
        type(T_2D_INT  )                    :: WetFaces
        integer, pointer, dimension(:,:,:)  :: WaterPoints3D, LandPoints3D, OpenPoints3D

        !Instance of ModuleTime
        integer                             :: ObjTime              = 0

        !Instance of ModuleHorizontalMap
        integer                             :: ObjHorizontalMap     = 0
       
        !Instance of Geometry
        integer                             :: ObjGeometry          = 0
        
        !Collection of instances                
        type(T_Map), pointer                :: Next
    end type T_Map 
    
        !Global Module Variables
    type (T_Map),    pointer                :: Me, FirstObjMap

    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructMap(Map_ID,                                           &
                            GeometryID,                                       &
                            HorizontalMapID,                                  &
                            TimeID, GridDataID, HorizontalGridID, STAT)  

        !Arguments-------------------------------------------------------------
        integer                             :: Map_ID
        integer                             :: HorizontalMapID
        integer                             :: GeometryID
        integer, optional                   :: TimeID
        integer, optional                   :: GridDataID
        integer, optional                   :: HorizontalGridID
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: ready_         
        integer, pointer, dimension(:, :)   :: WaterPoints2D, LandPoints2D, KFloorZ
                    
        !Local-----------------------------------------------------------------
        integer                             :: i, j, k, KFundo
        integer                             :: STAT_
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        logical                             :: IsolatedCell = .false.
        type(T_MapCell), pointer            :: FirstIsolatedCell

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_


        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mMap_)) then
            nullify (FirstObjMap)
            call RegisterModule (mMap_) 
        endif
        
        call Ready (Map_ID, ready_)

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance()
            
            !Associates External Instances

            if (present(TimeID)) then
                Me%ObjTime       = AssociateInstance (mTIME_,           TimeID          )
            endif

            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )

            !Gets the size from the Geometry
            call GetGeometrySize(Me%ObjGeometry,                                     &
                                   Size        = Me%Size,                            &
                                   WorkSize    = Me%WorkSize,                        &
                                   STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMap - ModuleMap - ERR01'

            !Auxiliar variables for the do loops
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB

            !Recieves the WaterPoints2D
            call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMap - ModuleMap - ERR02'

            !Recieves the LandPoints2D
            call GetLandPoints2D(Me%ObjHorizontalMap, LandPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMap - ModuleMap - ERR03'

            !Recieves KFloor
            call GetGeometryKFloor(Me%ObjGeometry, Z = KFloorZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMap - ModuleMap - ERR04'

            !Allocates Variables
            call AllocateVariables()

            !Constructs Land/Water points3D
do1:        do j = JLB, JUB
do2:        do i = ILB, IUB 

                if  (WaterPoints2D(i, j) == 1) then
                    KFUNDO = KFloorZ(i, j)
                    do K = KLB, KFUNDO - 1
                        Me%LandPoints3D(i, j, k)      = 1
                    enddo
                    do K = KFUNDO, KUB
                        Me%WaterPoints3D(i, j, k)     = 1
                    end do
                endif

                if  (LandPoints2D(i, j) == 1) then
                    do K = KLB, KUB
                        Me%LandPoints3D(i, j, k)      = 1
                    end do
                end if

            end do do2
            end do do1

            !This subroutine searches for the coastline
            call ConstructLandBoundaryFaces3D()

            if (present(GridDataID) .AND. present(HorizontalGridID)) then

                !Checks the existence of isolated cells which may become endlessly hotter
                !(cycle the cells and check the WaterPoints)
                do j = JLB, JUB
                do i = ILB, IUB 

                    if  (WaterPoints2D(i, j) == 1) then
                        KFUNDO = KFloorZ(i, j)
                        do k = KFUNDO, KUB !only waterpoints are checked

                            call CheckIsolatedCell(i, j, k, FirstIsolatedCell,      & 
                                 IsolatedCell) 

                        end do

                    endif

                end do 
                end do 

                if (IsolatedCell) then
                    
                    write(*,*)
                    write(*,*) 'Isolated cells in horizontal are not allowed.'
                    write(*,*)

                    !Construct a new bathymetry with no isolated cells
                    call CorrectIsolatedCells(GridDataID, HorizontalGridID,         &
                         FirstIsolatedCell)
               
                endif

            endif

            call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)           
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMap - ModuleMap - ERR06'

            call UnGetHorizontalMap(Me%ObjHorizontalMap, LandPoints2D, STAT = STAT_CALL)         
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMap - ModuleMap - ERR07'
                   
            call UnGetGeometry     (Me%ObjGeometry, KFloorZ, STAT = STAT_CALL)       
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMap - ModuleMap - ERR08'


            STAT_ = SUCCESS_
            
            !Returns ID
            Map_ID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleMap - ConstructMap - ERR10' 

        end if cd0

        if (present(STAT)) STAT = STAT_


        !----------------------------------------------------------------------

    end subroutine ConstructMap

    !--------------------------------------------------------------------------

    subroutine AllocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Map), pointer           :: NewMap
        type (T_Map), pointer           :: PreviousMap


        !Allocates new instance
        allocate (NewMap)
        nullify  (NewMap%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjMap)) then
            FirstObjMap      => NewMap
            Me               => NewMap
        else
            PreviousMap      => FirstObjMap
            Me               => FirstObjMap%Next
            do while (associated(Me))
                PreviousMap  => Me
                Me           => Me%Next
            enddo
            Me               => NewMap
            PreviousMap%Next => NewMap
        endif

        Me%InstanceID = RegisterNewInstance (mMAP_)

    end subroutine AllocateInstance    
    
    !--------------------------------------------------------------------------

    subroutine AllocateVariables() 

        !External--------------------------------------------------------------

        !Local----------------------------------------------------------------
        integer                         :: ILB, IUB
        integer                         :: JLB, JUB
        integer                         :: KLB, KUB
        integer                         :: STAT

        !----------------------------------------------------------------------
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        allocate(Me%WaterPoints3D       (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR01'
        allocate(Me%LandPoints3D        (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR02'
        allocate(Me%OpenPoints3D        (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR03'

        allocate(Me%ComputeFaces3D%U    (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR04'
        allocate(Me%ComputeFaces3D%V    (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR05'
        allocate(Me%ComputeFaces3D%W    (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR06'

        allocate(Me%LandBoundaryFaces3D%U (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR07'
        allocate(Me%LandBoundaryFaces3D%V (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR08'

        allocate(Me%ImposedNormalFaces%U    (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR09'

        allocate(Me%ImposedNormalFaces%V    (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR10'

        allocate(Me%ImposedTangentialFaces%U(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR11'

        allocate(Me%ImposedTangentialFaces%V(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR12'

        allocate(Me%WetFaces%U(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR13'

        allocate(Me%WetFaces%V(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT)
        if (STAT /= SUCCESS_) stop 'AllocateVariables - ModuleMap - ERR14'


        !By default all values are zero
        Me%WaterPoints3D                    = 0
        Me%LandPoints3D                     = 0
        Me%OpenPoints3D                     = 0
                                                
        Me%ComputeFaces3D%U                 = 0
        Me%ComputeFaces3D%V                 = 0
        Me%ComputeFaces3D%W                 = 0
                                                
        Me%LandBoundaryFaces3D%U            = 0
        Me%LandBoundaryFaces3D%V            = 0


        Me%ImposedNormalFaces%U(:,:,:)      = 0
        Me%ImposedNormalFaces%V(:,:,:)      = 0
                                                
        Me%ImposedTangentialFaces%U(:,:,:)  = 0
        Me%ImposedTangentialFaces%V(:,:,:)  = 0

        Me%WetFaces%U(:,:,:)                = 0
        Me%WetFaces%V(:,:,:)                = 0


        !----------------------------------------------------------------------

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine ConstructLandBoundaryFaces3D()   

        !Arguments--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                         :: i, j, k
        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (    (    (Me%WaterPoints3D(i, j-1, k) == 1) .and.        &   !Left Water
                         (Me%LandPoints3D (i, j, k)   == 1)   )  .or.    &   !Right Land
                    (    (Me%WaterPoints3D(i, j, k)   == 1) .and.        &   !Right Water
                         (Me%LandPoints3D (i, j-1, k) == 1)   )       )  &   !Left Land

                Me%LandBoundaryFaces3D%U(i, j, k) = 1


            if (    (    (Me%WaterPoints3D(i-1, j, k) == 1) .and.        &   !Down Water
                         (Me%LandPoints3D (i, j, k)   == 1)   )  .or.    &   !Up Land
                    (    (Me%WaterPoints3D(i, j, k)   == 1) .and.        &   !Up Water
                         (Me%LandPoints3D (i-1, j, k) == 1)   )       )  &   !Down Land

                Me%LandBoundaryFaces3D%V(i, j, k) = 1

        enddo
        enddo
        enddo

        !----------------------------------------------------------------------

    end subroutine ConstructLandBoundaryFaces3D

        !----------------------------------------------------------------------

    subroutine CheckIsolatedCell(i, j, k, FirstIsolatedCell, IsolatedCell)  

        !Arguments-------------------------------------------------------------
        type(T_MapCell), pointer            :: FirstIsolatedCell
        logical                             :: IsolatedCell
        integer                             :: i, j, k

        !External--------------------------------------------------------------
                    
        !Local-----------------------------------------------------------------
        type(T_MapCell), pointer            :: NewIsolatedCell
        type(T_MapCell), pointer            :: ObjIsolatedCell, PreviousIsolatedCell

        !----------------------------------------------------------------------

        if (Me%WaterPoints3D(i,j-1,k) == 0) then
            if (Me%WaterPoints3D(i,j+1,k) == 0) then
                if (Me%WaterPoints3D(i-1,j,k) == 0) then
                    if (Me%WaterPoints3D(i+1,j,k) == 0) then

900                     FORMAT(i4,1x,i4,1x,i4)
                        write(*,*)'Isolated cell in horizontal (i,j,k):'
                        write(*,fmt = 900) i, j, k

                        if (.not. IsolatedCell) then
                            IsolatedCell = .true.
                        endif

                        allocate (NewIsolatedCell)
                        nullify (NewIsolatedCell%Next)

                        if (.not. associated(FirstIsolatedCell)) then
                            FirstIsolatedCell => NewIsolatedCell
                            ObjIsolatedCell => NewIsolatedCell
                        else
                            PreviousIsolatedCell => FirstIsolatedCell
                            ObjIsolatedCell => FirstIsolatedCell%Next
                            do while (associated(ObjIsolatedCell))
                                PreviousIsolatedCell => ObjIsolatedCell
                                ObjIsolatedCell => ObjIsolatedCell%Next
                            enddo
                            ObjIsolatedCell => NewIsolatedCell
                            PreviousIsolatedCell%Next => NewIsolatedCell
                        end if

                        ObjIsolatedCell%Icell = i
                        ObjIsolatedCell%Jcell = j

                    endif
                endif
            endif
        endif

    end subroutine CheckIsolatedCell

        !----------------------------------------------------------------------

    subroutine CorrectIsolatedCells(GridDataID, HorizontalGridID, FirstIsolatedCell)  

        !Arguments-------------------------------------------------------------
        integer                             :: GridDataID
        integer                             :: HorizontalGridID
        type(T_MapCell), pointer            :: FirstIsolatedCell

        !External--------------------------------------------------------------
        real, dimension(:, :), pointer      :: Bathymetry, NewBathymetry
        character(len=StringLength)         :: BathymetryFile
                    
        !Local-----------------------------------------------------------------
        integer                             :: i, j
        type(T_MapCell), pointer            :: ObjIsolatedCell
        logical                             :: Distortion
        character(len=StringLength)         :: Comment1, Comment2
        integer                             :: LengthWithoutExt
        type (T_Size2D)                     :: Size2D
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------

        call GetGridData(GridDataID, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CorrectIsolatedCells - ModuleMap - ERR01'

        nullify (NewBathymetry)
        allocate(NewBathymetry(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        !Copy bathymetry array
        NewBathymetry(:,:) = Bathymetry(:,:)

        ObjIsolatedCell => FirstIsolatedCell
        do while (associated(ObjIsolatedCell))
            i = ObjIsolatedCell%Icell
            j = ObjIsolatedCell%Jcell
            NewBathymetry(i,j) = max(NewBathymetry(i,j+1),                          &
                NewBathymetry(i,j-1), NewBathymetry(i-1,j),                         &
                NewBathymetry(i+1,j))
            ObjIsolatedCell => ObjIsolatedCell%Next
            write(*,*)'i              = ', i
            write(*,*)'j              = ', j
            write(*,*)'Bathymetry     = ', Bathymetry(i, j)
            write(*,*)'New Bathymetry = ',NewBathymetry(i, j)
        enddo


        !Writes new Bathymetry
        call GetGridDataFileName(GridDataID, FileName = BathymetryFile,       & 
             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CorrectIsolatedCells - ModuleMap - ERR11'

        BathymetryFile  = adjustl(BathymetryFile)
        Comment1        = "Automatic Generated Grid Data File"
        Comment2        = "Based On "//trim(BathymetryFile)
        LengthWithoutExt= len_trim(BathymetryFile) - 4
        BathymetryFile  = BathymetryFile(1:LengthWithoutExt)//"_"//".new2"

        call GetCheckDistortion (HorizontalGridID, Distortion, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CorrectIsolatedCells - ModuleMap - ERR12'

        Size2D%ILB = Me%WorkSize%ILB
        Size2D%IUB = Me%WorkSize%IUB
        Size2D%JLB = Me%WorkSize%JLB
        Size2D%JUB = Me%WorkSize%JUB
                
        call WriteGridData(FileName         = BathymetryFile,                           &
                           COMENT1          = Comment1,                                 &      
                           COMENT2          = Comment2,                                 &      
                           HorizontalGridID = HorizontalGridID,                         &
                           FillValue        = -99.,                                     &
                           Overwrite        = ON,                                       &
                           GridData2D_Real  = NewBathymetry,                            &
                           STAT             = STAT_CALL)        
                           
        if (STAT_CALL /= SUCCESS_) stop 'CorrectIsolatedCells - ModuleMap - ERR14'

        deallocate(NewBathymetry)
        nullify   (NewBathymetry)

        !Disposes pointer to the Bathymetry
        call UngetGridData(GridDataID, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CorrectIsolatedCells - ModuleMap - ERR15'

        !Displays Message to inform
        write(*,*)'A new Bathymetry has been created, with isolated cells removed'
        write(*,*)'Modify the file Nomfich.dat and Re-run the model'
        write(*,*)'New Bathymetry file : ', trim(BathymetryFile)
        stop

    end subroutine CorrectIsolatedCells


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !Updates the compute faces to taking in consideration a moving boundary condition
    subroutine UpdateMovBoundCompFaces3D(Map_ID, SurfaceElevation, ActualTime, STAT)   
        
        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        real, dimension(:, :), pointer              :: SurfaceElevation
        type (T_Time)                               :: ActualTime
        integer, optional, intent(OUT)              :: STAT      
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_          
        integer                                     :: STAT_CALL
        integer, dimension(:, :), pointer           :: ComputeFaces2DU, ComputeFaces2DV
        integer, dimension(:, :), pointer           :: KFloorU, KFloorV
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB, i, j, k
        real                                        :: MinWaterColumn
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then


            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB

            !Actualize the time
            Me%ActualTime       = ActualTime

            !Initializes (U, V, W)
            Me%ComputeFaces3D%U = 0
            Me%ComputeFaces3D%V = 0
            Me%ComputeFaces3D%W = 0

            !Gets Min WaterColumn
            call GetGeometryMinWaterColumn(Me%ObjGeometry, MinWaterColumn, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateMovBoundCompFaces3D - ModuleMap - ERR02'

            !Updates Compute Faces 2D before updating the 3D ones
            call UpdateComputeFaces2D(Me%ObjHorizontalMap, SurfaceElevation,             &
                                      MinWaterColumn, ActualTime,                        &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateMovBoundCompFaces3D - ModuleMap - ERR03'

            !Gets ComputeFaces2D
            call GetComputeFaces2D(Me%ObjHorizontalMap, ComputeFaces2DU,                 &
                                   ComputeFaces2DV, ActualTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateMovBoundCompFaces3D - ModuleMap - ERR04'


            !Gets KFloorU, KFloorV
            call GetGeometryKFloor(Me%ObjGeometry, U = KFloorU, V = KFloorV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateMovBoundCompFaces3D - ModuleMap - ERR05'

            if (MonitorPerformance) call StartWatch ("ModuleMap", "UpdateMovBoundCompFaces3D")

            CHUNK = CHUNK_J(JLB,JUB)
            !$OMP PARALLEL PRIVATE(I,J)

            !ComputeFacesU / ComputeFacesV
            do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB
                if (ComputeFaces2DU(i, j) == 1 .and. k >= KFloorU(i, j))                &
                    Me%ComputeFaces3D%U(i, j, k) = 1

                if (ComputeFaces2DV(i, j) == 1 .and. k >= KFloorV(i, j))                &
                    Me%ComputeFaces3D%V(i, j, k) = 1

            enddo
            enddo
            !$OMP END DO
            enddo

            !$OMP END PARALLEL

            if (MonitorPerformance) call StopWatch ("ModuleMap", "UpdateMovBoundCompFaces3D")

            !Updates all ComputefacesW
            call UpdateComputeFacesW()
        
            !Updates all points which have at least one Computeface
            call UpdateOpenPoints3D()

            !Updates all the imposed points
            call UpdateImposedValues()

            call UnGetHorizontalMap(Me%ObjHorizontalMap, ComputeFaces2DU, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateMovBoundCompFaces3D - ModuleMap - ERR07'

            call UnGetHorizontalMap(Me%ObjHorizontalMap, ComputeFaces2DV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateMovBoundCompFaces3D - ModuleMap - ERR08'

            call UnGetGeometry(Me%ObjGeometry, KFloorU, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateMovBoundCompFaces3D - ModuleMap - ERR09'

            call UnGetGeometry(Me%ObjGeometry, KFloorV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateMovBoundCompFaces3D - ModuleMap - ERR10'


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UpdateMovBoundCompFaces3D

    !--------------------------------------------------------------------------
    
    subroutine UpdateSoilComputeFaces3D(Map_ID, Prop3D, PropMin, BoundaryType,          & 
                                        Imposed, STAT)   
        
        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        real,    dimension(:,:,:), pointer          :: Prop3D, PropMin
        integer, dimension(:,:,:), pointer          :: BoundaryType
        integer,           Intent(IN )              :: Imposed
        integer, optional, intent(OUT)              :: STAT      

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_          
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k


        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB


            !Initializes (U, V, W)
            Me%ComputeFaces3D%U = 0
            Me%ComputeFaces3D%V = 0
            Me%ComputeFaces3D%W = 0


            !ComputeFacesU / ComputeFacesV
            do k = KLB  , KUB
            do j = JLB+1, JUB
            do i = ILB  , IUB

                if (Me%WaterPoints3D(i, j, k) == 1 .and.                                & 
                    Me%WaterPoints3D(i, j-1, k) == 1) then
                
                    if (.not.((Prop3D      (i, j, k)   <= PropMin (i, j, k) .and.       &
                               Prop3D      (i, j-1, k) <= PropMin(i, j-1, k))   .or.    &
                              (BoundaryType(i, j, k) == Imposed           .and.         &
                               BoundaryType(i, j-1, k) == Imposed))) then
                        Me%ComputeFaces3D%U(i, j, k) = 1
                    endif
                endif

            enddo
            enddo
            enddo

            do k = KLB  , KUB
            do j = JLB  , JUB
            do i = ILB+1, IUB

                if (Me%WaterPoints3D(i, j, k) == 1 .and.                                &
                    Me%WaterPoints3D(i-1, j, k) == 1) then
                
                    if (.not.((Prop3D      (i, j, k)   <= PropMin(i, j, k)  .and.       &
                               Prop3D      (i-1, j, k) <= PropMin(i-1, j, k))   .or.    &
                              (BoundaryType(i, j, k) == Imposed          .and.          &
                               BoundaryType(i-1, j, k) == Imposed))) then
                        Me%ComputeFaces3D%V(i, j, k) = 1
                    endif
                endif

            enddo
            enddo
            enddo


            do k = KLB+1, KUB
            do j = JLB  , JUB
            do i = ILB  , IUB
                !
                if (Me%WaterPoints3D(i, j, k) == 1 .and. Me%WaterPoints3D(i, j, k-1) == 1) then
                
                    if (.not.((Prop3D      (i, j, k) <= PropMin(i, j, k)  .and.          &
                               Prop3D      (i, j, k-1) <= PropMin(i, j, k-1))   .or.     &
                              (BoundaryType(i, j, k) == Imposed          .and.           &
                               BoundaryType(i, j, k-1) == Imposed))) then
                        Me%ComputeFaces3D%W(i, j, k) = 1

                    endif
                endif

            enddo
            enddo
            enddo
       
            !Updates all points which have at least one Computeface
            call UpdateOpenPoints3D()

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UpdateSoilComputeFaces3D

    !--------------------------------------------------------------------------

    subroutine UpdateUnsaturatedComputeFaces3D (Map_ID, DummyR, DummyI, STAT)   

        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        integer                                     :: DummyI
        real                                        :: DummyR
        integer, optional, intent(OUT)              :: STAT      

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_          
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k, STAT_CALL
        integer, dimension(:,:), pointer            :: KFloorZ


        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)

        !Dummy existe to distinguish UpdateUnsaturatedComputeFaces3D interface from other
        !interfaces
        DummyR = 1.
        DummyI = 1

        if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB


            !Initializes (U, V, W)
            Me%ComputeFaces3D%U = 0
            Me%ComputeFaces3D%V = 0
            Me%ComputeFaces3D%W = 0

            !Gets KFloorZ
            call GetGeometryKFloor(Me%ObjGeometry, Z = KFloorZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateComputeFacesW - ModuleMap - ERR10'


            !ComputeFacesU / ComputeFacesV            
            do j = JLB+1, JUB
            do i = ILB  , IUB
                if (Me%WaterPoints3D(i, j, KUB) == 1) then
                    do k = KFloorZ(i,j), KUB
                        if (Me%WaterPoints3D(i, j, k) == 1 .and. Me%WaterPoints3D(i, j-1, k) == 1) then
                            Me%ComputeFaces3D%U(i, j  , k) = 1
                        endif
                    enddo
                endif
            enddo
            enddo
            
            do j = JLB  , JUB
            do i = ILB+1, IUB
                if (Me%WaterPoints3D(i, j, KUB) == 1) then
                    do k = KFloorZ(i,j), KUB
                        if (Me%WaterPoints3D(i, j, k) == 1 .and. Me%WaterPoints3D(i-1, j, k) == 1) then
                            Me%ComputeFaces3D%V(i, j, k) = 1
                        endif
                    enddo
                endif
            enddo
            enddo

            do j = JLB  , JUB
            do i = ILB  , IUB
                if (Me%WaterPoints3D(i, j, KUB) == 1) then
                    do k = KFloorZ(i,j)+1, KUB
                        Me%ComputeFaces3D%W(i, j, k) = 1
                    enddo
                endif
            enddo
            enddo

            call UnGetGeometry     (Me%ObjGeometry, KFloorZ, STAT = STAT_CALL)       
            if (STAT_CALL /= SUCCESS_) stop 'UpdateComputeFacesW - ModuleMap - ERR20'

            !Updates all points which have at least one Computeface
            call UpdateOpenPoints3D()

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    endsubroutine UpdateUnsaturatedComputeFaces3D

    !--------------------------------------------------------------------------


    subroutine UpdateSedimentCompFaces3D(Map_ID, STAT)   
        
        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        integer, optional, intent(OUT)              :: STAT      

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_          
        integer                                     :: i, j, k

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            do k = Me%WorkSize%KLB  , Me%WorkSize%KUB
            do j = Me%WorkSize%JLB  , Me%WorkSize%JUB
            do i = Me%WorkSize%ILB  , Me%WorkSize%IUB

                if (Me%WaterPoints3D(i, j, k) == 1) then

                    Me%ComputeFaces3D%U(i, j, k) = 0
                    Me%ComputeFaces3D%V(i, j, k) = 0
                    Me%ComputeFaces3D%W(i, j, k) = 1

                end if

            enddo
            enddo
            enddo

            call UpdateOpenPoints3D()

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UpdateSedimentCompFaces3D

    !--------------------------------------------------------------------------

    subroutine UpdateComputeFacesW()

        !Arguments-------------------------------------------------------------
                    
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer            :: KFloorZ
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k, STAT_CALL
        integer                                     :: Chunk

        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        !Recieves KFloor
        call GetGeometryKFloor(Me%ObjGeometry, Z = KFloorZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UpdateComputeFacesW - ModuleMap - ERR10'

        if (MonitorPerformance) call StartWatch ("ModuleMap", "UpdateComputeFacesW")

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(i,j,k)

        !If the surface cell have at least one horizontal computeface, so the upper face
        !will also be a compute face
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ComputeFaces3D%U(i, j, KUB) + Me%ComputeFaces3D%U(i, j+1, KUB) +    &
                Me%ComputeFaces3D%V(i, j, KUB) + Me%ComputeFaces3D%V(i+1, j, KUB) > 0) then
                !The face close to the surface isnt ComputeFace    
                do k = KFloorZ(i,j) + 1, KUB
                    Me%ComputeFaces3D%W(i, j, k) = 1
                end do
            endif
        enddo
        enddo
        !$OMP END DO     
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleMap", "UpdateComputeFacesW")

        call UnGetGeometry     (Me%ObjGeometry, KFloorZ, STAT = STAT_CALL)       
        if (STAT_CALL /= SUCCESS_) stop 'UpdateComputeFacesW - ModuleMap - ERR20'


    end subroutine UpdateComputeFacesW

    !--------------------------------------------------------------------------

    subroutine CleanBoundary()

        !Arguments-------------------------------------------------------------
                            
        !Local-----------------------------------------------------------------
        integer, dimension(:, :), pointer           :: BoundaryPoints2D
        integer, dimension(:, :), pointer           :: ExteriorFaceU, ExteriorFaceV
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        !Gets BoundaryPoints2D
        call GetBoundaries(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CleanBoundary - ModuleMap - ERR01'  

        !Gets ExteriorPoints2D
        call GetExteriorBoundaryFaces(Me%ObjHorizontalMap, ExteriorFaceU,            &
                                      ExteriorFaceV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CleanBoundary - ModuleMap - ERR02'  

        !Eliminate ComputeFaces U along boundary
        do k = KLB,   KUB
        do j = JLB+1, JUB
        do i = ILB,   IUB
            if ((BoundaryPoints2D(i, j) == 1) .and. (BoundaryPoints2D(i, j-1) == 1))then
                 Me%ComputeFaces3D%U(I, J, K) = 0
            end if
        enddo
        enddo
        enddo

        !Eliminate ComputeFaces V along boundary
        do k = KLB,   KUB
        do j = JLB,   JUB   
        do i = ILB+1, IUB
            if ((BoundaryPoints2D(i,   j) == 1) .and. (BoundaryPoints2D(i-1, j) == 1))then
                 Me%ComputeFaces3D%V(I, J, K) = 0
            end if
        enddo
        enddo
        enddo


        !Eliminate Exterior Faces U 
        do k = KLB,   KUB
        do j = JLB,   JUB + 1
        do i = ILB,   IUB
            if (ExteriorFaceU(i, j) == 1)then
                 Me%ComputeFaces3D%U(I, J, K) = 0
            endif
        enddo
        enddo
        enddo

        !Eliminate Exterior Faces V 
        do k = KLB,   KUB
        do j = JLB,   JUB   
        do i = ILB,   IUB + 1
            if (ExteriorFaceV(i, j) == 1)then
                 Me%ComputeFaces3D%V(I, J, K) = 0
            endif
        enddo
        enddo
        enddo

        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CleanBoundary - ModuleMap - ERR03'  

        call UnGetHorizontalMap(Me%ObjHorizontalMap, ExteriorFaceU,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CleanBoundary - ModuleMap - ERR04'  

        call UnGetHorizontalMap(Me%ObjHorizontalMap, ExteriorFaceV,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CleanBoundary - ModuleMap - ERR05'  

    end subroutine CleanBoundary

    !--------------------------------------------------------------------------

    subroutine UpdateOpenPoints3D()   
        
        !Arguments-------------------------------------------------------------
                    
        !Local-----------------------------------------------------------------
        integer :: ILB, IUB, JLB, JUB, KLB, KUB, i, j, k
        integer :: CHUNK
                                
        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        if (MonitorPerformance) call StartWatch ("ModuleMap", "UpdateOpenPoints3D")

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(k,i,j)

        !Searches for open
        do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ComputeFaces3D%U(i, j, k) + Me%ComputeFaces3D%U(i, j+1, k) +  & 
                Me%ComputeFaces3D%V(i, j, k) + Me%ComputeFaces3D%V(i+1, j, k) +  & 
                Me%ComputeFaces3D%W(i, j, k) + Me%ComputeFaces3D%W(i, j, k+1) > 0 ) then
                Me%OpenPoints3D(i, j, k) = 1
            else
                Me%OpenPoints3D(i, j, k) = 0
            endif
        enddo
        enddo
        !$OMP END DO

        !obtain the U WetFaces
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB+1, JUB
        do i = ILB, IUB
            if (Me%OpenPoints3D(i, j, k) == 1 .or. Me%OpenPoints3D(i, j-1, k) == 1) then
                Me%WetFaces%U(i, j, k) = 1
            else
                Me%WetFaces%U(i, j, k) = 0       
            endif
        
            if (Me%OpenPoints3D(i, JLB, k) == 1) then
                Me%WetFaces%U(i, JLB, k) = 1
            else
                Me%WetFaces%U(i, JLB, k) = 0       
            endif
        enddo
        enddo
        !$OMP END DO

        !obtain the V WetFaces
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB+1, IUB
            if (Me%OpenPoints3D(i, j, k) == 1 .or. Me%OpenPoints3D(i-1, j, k) == 1) then
                Me%WetFaces%V(i, j, k) = 1
            else
                Me%WetFaces%V(i, j, k) = 0       
            endif
        enddo
        
        if (Me%OpenPoints3D(ILB, j, k) == 1) then
            Me%WetFaces%V(ILB, j, k) = 1
        else
            Me%WetFaces%V(ILB, j, k) = 0       
        endif
        enddo
        !$OMP END DO
    
        enddo
        
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleMap", "UpdateOpenPoints3D")

        !----------------------------------------------------------------------

    end subroutine UpdateOpenPoints3D

    !--------------------------------------------------------------------------

    subroutine UpdateImposedValues() 

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer, dimension(:,:), pointer    :: BoundaryPoints2D, BoundaryFacesU
        integer, dimension(:,:), pointer    :: BoundaryFacesV, KFloorU, KFloorV
        integer, dimension(:,:), pointer    :: ExteriorFacesU, ExteriorFacesV
        integer                             :: STATUS
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: i, j, k, kbottom
        integer                             :: CHUNK
        !logical                             :: ErrorOcurred1, ErrorOcurred2

        !-----------------------------------------------------------------------


        !Gets BoundaryPoints2D
        call GetBoundaries(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STATUS)
        if (STATUS /= SUCCESS_)stop 'UpdateImposedValues - ModuleMap - ERR01'


        !Gets BoundaryFaces
        call GetBoundaryFaces(Me%ObjHorizontalMap, BoundaryFacesU, BoundaryFacesV, STAT = STATUS)
        if (STATUS /= SUCCESS_)stop 'UpdateImposedValues - ModuleMap - ERR02'

        !Gets KFloorU, KFloorV
        call GetGeometryKFloor(Me%ObjGeometry, U = KFloorU, V = KFloorV, STAT = STATUS)
        if (STATUS /= SUCCESS_)stop 'UpdateImposedValues - ModuleMap - ERR03'

        call GetExteriorBoundaryFaces(Me%ObjHorizontalMap,                           &
                                      BoundaryPointsFaceU = ExteriorFacesU,              &
                                      BoundaryPointsFaceV = ExteriorFacesV,              &
                                      STAT = STATUS)
        if (STATUS /= SUCCESS_)stop 'UpdateImposedValues - ModuleMap - ERR04'


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        Me%ImposedNormalFaces%U(:,:,:) = 0
        Me%ImposedNormalFaces%V(:,:,:) = 0

        Me%ImposedTangentialFaces%U(:,:,:) = 0
        Me%ImposedTangentialFaces%V(:,:,:) = 0

        if (MonitorPerformance) call StartWatch ("ModuleMap", "UpdateImposedValues")

        !ErrorOcurred1 = .false.
        !ErrorOcurred2 = .false.

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(kbottom,k,i,j)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do1 :   do j = JLB, JUB + 1
do2 :   do i = ILB, IUB
                   
if1:        if (BoundaryPoints2D(i, j    ) == 1 .and.                                    & 
                BoundaryPoints2D(i, j - 1) == 1)  then 

                kbottom = KFloorU(i, j)

dok1:           do k = kbottom, KUB

                    Me%ImposedTangentialFaces%U(i, j, k) = 1            

                enddo dok1

            endif if1


if45 :      if (ExteriorFacesU (i, j) == 1 .and.                                         &
                (BoundaryFacesU (i-1, j)  == 1 .or. BoundaryFacesU (i+1, j)  == 1 ))  then 
                
                kbottom = KFloorU(i, j)

dok45:          do k = kbottom, KUB

                        Me%ImposedTangentialFaces%U(i, j, k) = 1            

                                
                enddo dok45

            endif if45



if2:        if (BoundaryFacesU (i, j)  == 1) then
                
dok2:           do k = KLB, KUB
                    
if3:                if(Me%ComputeFaces3D%U(i, j, k)  == 1) then 

                        if  (BoundaryPoints2D (i, j) == BoundaryPoints2D (i, j-1)) then

                            !$OMP CRITICAL (ERR05)
                            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR05")        
                            !$OMP END CRITICAL (ERR05)
                        endif

                        !if  (BoundaryPoints2D (i, j) == BoundaryPoints2D (i, j-1)) then
                            !ErrorOcurred1 = .true.
                            !exit do1
                        !endif

                        if  (BoundaryPoints2D (i, j    )  == 1)                          &
                              Me%ImposedNormalFaces%U(i, j + 1, k)  = 1

                        if  (BoundaryPoints2D (i, j - 1)  == 1)                          &
                              Me%ImposedNormalFaces%U(i, j - 1, k)  = 1                    

                    endif if3

                enddo dok2 

            endif if2

        end do do2
        end do do1
        !$OMP END DO

        !if (ErrorOcurred1)                                                              &
        !    call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR05")        

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :   do j = JLB, JUB
do4 :   do i = ILB, IUB + 1
                   
if4:        if (BoundaryPoints2D(i    , j) == 1 .and.                                    & 
                BoundaryPoints2D(i - 1, j) == 1)  then 

                kbottom = KFloorV(i, j)

dok3:           do k = kbottom, KUB
                    
                    Me%ImposedTangentialFaces%V(i, j, k) = 1            
                    
                enddo dok3

            endif if4



if46 :      if (ExteriorFacesV (i, j) == 1 .and.                                         &
                (BoundaryFacesV (i, j-1)  == 1 .or. BoundaryFacesV (i, j+1)  == 1 ))  then 
                
                kbottom = KFloorV(i, j)

dok46:          do k = kbottom, KUB

                        Me%ImposedTangentialFaces%V(i, j, k) = 1            

                                
                enddo dok46

            endif if46


if6 :       if (BoundaryFacesV (i, j)  == 1) then

dok4:           do k = KLB, KUB
                    
if5 :                if(Me%ComputeFaces3D%V(i, j, k)  == 1) then 

                        if  (BoundaryPoints2D (i, j) == BoundaryPoints2D (i-1, j)) then
                            !ErrorOcurred1 = .true.
                            !exit do3
                            !$OMP CRITICAL (ERR06)
                            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR06")        
                            !$OMP END CRITICAL (ERR06)
                        endif

                        if  (BoundaryPoints2D (i, j    )  == 1)                          &
                              Me%ImposedNormalFaces%V(i + 1, j, k)  = 1

                        if  (BoundaryPoints2D (i - 1, j)  == 1)                          &
                              Me%ImposedNormalFaces%V(i - 1, j, k)  = 1                    

                    endif if5

                enddo dok4

            endif if6

        end do do4
        end do do3
        !$OMP END DO

        !if (ErrorOcurred1)                                                       &
        !    call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR06")        

do5:    do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do6 :   do j = JLB, JUB
do7 :   do i = ILB, IUB + 1

            if (Me%ImposedTangentialFaces%V(i, j, k) == 1 .and.                      &
                Me%ComputeFaces3D%V        (i, j, k) == 1) then
                !ErrorOcurred1 = .true.
                !exit do6
                
                !$OMP CRITICAL (ERR07)
                call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR07")        
                !$OMP END CRITICAL (ERR07)
            endif

            if (Me%ImposedNormalFaces%V(i, j, k) == 1 .and.                          &
                Me%ComputeFaces3D%V    (i, j, k) == 1) then
                !ErrorOcurred2 = .true.
                !exit do6
                
                !$OMP CRITICAL (ERR08)
                call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR08")        
                !$OMP END CRITICAL (ERR08)
            endif

        end do do7
        end do do6      
        !$OMP END DO NOWAIT

        !if (ErrorOcurred1)                                                       &
        !    call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR07")        

        !if (ErrorOcurred2)                                                       &
        !    call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR08")        
        
        end do do5


do8 :   do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do9 :   do j = JLB, JUB + 1
do10:   do i = ILB, IUB

            if (Me%ImposedTangentialFaces%U(i, j, k) == 1 .and.                      &
                Me%ComputeFaces3D%U        (i, j, k) == 1) then
                !ErrorOcurred1 = .true.
                !exit do9
                
                !$OMP CRITICAL (ERR09)
                call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR09")        
                !$OMP END CRITICAL (ERR09)
            endif

            if (Me%ImposedNormalFaces%U(i, j, k) == 1 .and.                          &
                Me%ComputeFaces3D%U    (i, j, k) == 1) then
                !ErrorOcurred1 = .true.
                !exit do9
                
                !$OMP CRITICAL (ERR10)
                call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR10")        
                !$OMP END CRITICAL (ERR10)
            endif

        end do do10
        end do do9
        !$OMP END DO
        
        !if (ErrorOcurred1)                                                           &
        !    call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR09")        

        !if (ErrorOcurred2)                                                           &
        !    call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR10")        
        
        end do do8
        

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleMap", "UpdateImposedValues")

        !UnGet horizontal ModuleMap
        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STATUS)
        if (STATUS /= SUCCESS_)                                                          &  
            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR11")        

        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryFacesU, STAT = STATUS)

        if (STATUS /= SUCCESS_)                                                          &  
            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR12")        

        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryFacesV, STAT = STATUS)

        if (STATUS /= SUCCESS_)                                                          &  
            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR13")        

        call UnGetHorizontalMap(Me%ObjHorizontalMap, ExteriorFacesU, STAT = STATUS)

        if (STATUS /= SUCCESS_)                                                          &  
            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR14")        

        call UnGetHorizontalMap(Me%ObjHorizontalMap, ExteriorFacesV, STAT = STATUS)

        if (STATUS /= SUCCESS_)                                                          &  
            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR15")        



        !UnGets KFloorU, KFloorV
        call UnGetGeometry(Me%ObjGeometry, KFloorU, STAT = STATUS)
        if (STATUS /= SUCCESS_)                                                          &  
            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR16")        


        call UnGetGeometry(Me%ObjGeometry, KFloorV, STAT = STATUS)
        if (STATUS /= SUCCESS_)                                                          &  
            call SetError (FATAL_, INTERNAL_, "UpdateImposedValues - ModuleMap - ERR17")        

    end subroutine UpdateImposedValues     


    !--------------------------------------------------------------------------

    
    subroutine SetComputesFaces3D(Map_ID, ComputeFacesU3D, ComputeFacesV3D,                &
                                  ActualTime, STAT)   
        
        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        integer, dimension(:, :, :), pointer        :: ComputeFacesU3D
        integer, dimension(:, :, :), pointer        :: ComputeFacesV3D
        type (T_Time)                               :: ActualTime
        integer, optional, intent(OUT)              :: STAT      

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            call SetMatrixValue(Me%ComputeFaces3D%U, Me%Size, ComputeFacesU3D)
            call SetMatrixValue(Me%ComputeFaces3D%V, Me%Size, ComputeFacesV3D)

            !Cleans Compute Faces along boundary
            call CleanBoundary()

            !Updates all ComputefacesW
            call UpdateComputeFacesW()

            !Updates all points which have at least one Computeface
            call UpdateOpenPoints3D ()

            Me%ActualTime = ActualTime

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
            

        !----------------------------------------------------------------------

    end subroutine SetComputesFaces3D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetLandPoints3D(Map_ID, LandPoints3D, STAT)  

        !Arguments-------------------------------------------------------------
        integer                                 :: Map_ID
        integer,  pointer, dimension(:,:,:)     :: LandPoints3D
        integer, optional, intent(OUT)          :: STAT

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mMAP_, Me%InstanceID)

            LandPoints3D => Me%LandPoints3D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_
            
        !----------------------------------------------------------------------

    end subroutine GetLandPoints3D

    
    !--------------------------------------------------------------------------

    
    subroutine GetWaterPoints3D(Map_ID, WaterPoints3D, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                 :: Map_ID
        integer,  pointer, dimension(:,:,:)     :: WaterPoints3D
        integer,  optional, intent(OUT)         :: STAT

        !External--------------------------------------------------------------
        integer                                 :: ready_              

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mMAP_, Me%InstanceID)

            WaterPoints3D => Me%WaterPoints3D


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if(present(STAT)) STAT = STAT_
           
        !----------------------------------------------------------------------

    end subroutine GetWaterPoints3D

    !--------------------------------------------------------------------------

    subroutine GetOpenPoints3D(Map_ID, OpenPoints3D, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                 :: Map_ID
        integer,  pointer, dimension(:,:,:)     :: OpenPoints3D
        integer, optional, intent(OUT)          :: STAT

        !External--------------------------------------------------------------
        integer                                 :: ready_              

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mMAP_, Me%InstanceID)

            OpenPoints3D => Me%OpenPoints3D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetOpenPoints3D

    
    !--------------------------------------------------------------------------


    subroutine GetComputeFaces3D(Map_ID, ComputeFacesU3D, ComputeFacesV3D, ComputeFacesW3D, &
                                 ActualTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Map_ID
        integer, optional, intent(OUT)                  :: STAT
        integer, dimension(:,:,:), optional, pointer    :: ComputeFacesU3D
        integer, dimension(:,:,:), optional, pointer    :: ComputeFacesV3D
        integer, dimension(:,:,:), optional, pointer    :: ComputeFacesW3D
        type(T_Time), optional                          :: ActualTime

        !External--------------------------------------------------------------
        integer                                         :: ready_              

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_

        !----------------------------------------------------------------------

        call Ready (Map_ID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(ActualTime)) then
                if (ActualTime .ne. Me%ActualTime) &
                    stop 'GetComputeFaces3D - ModuleMap - ERR01'
            endif

            if (present(ComputeFacesU3D)) then
                call Read_Lock(mMAP_, Me%InstanceID)
                ComputeFacesU3D => Me%ComputeFaces3D%U
            endif

            if (present(ComputeFacesV3D)) then
                call Read_Lock(mMAP_, Me%InstanceID)
                ComputeFacesV3D => Me%ComputeFaces3D%V
            endif

            if (present(ComputeFacesW3D)) then
                call Read_Lock(mMAP_, Me%InstanceID)
                ComputeFacesW3D => Me%ComputeFaces3D%W
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if(present(STAT))STAT = STAT_
            
        !----------------------------------------------------------------------

    end subroutine GetComputeFaces3D

    
    !-------------------------------------------------------------------------

    
    subroutine GetImposedNormalFaces(Map_ID, ImposedNormalFaceU,                         &
                                     ImposedNormalFaceV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        integer, dimension(:,:,:), pointer          :: ImposedNormalFaceU, ImposedNormalFaceV
        integer, optional, intent (OUT)             :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_              

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mMAP_, Me%InstanceID)!ImposedNormalFacesX
            call Read_Lock(mMAP_, Me%InstanceID)!ImposedNormalFacesY


            ImposedNormalFaceU => Me%ImposedNormalFaces%U
            ImposedNormalFaceV => Me%ImposedNormalFaces%V


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if(present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetImposedNormalFaces

    
    !--------------------------------------------------------------------------

    
    subroutine GetImposedTangentialFaces(Map_ID, ImposedTangentialFaceU,                 &
                                         ImposedTangentialFaceV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        integer, dimension(:,:,:), pointer          :: ImposedTangentialFaceU, ImposedTangentialFaceV
        integer, optional, intent (OUT)             :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_              

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)

if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mMAP_, Me%InstanceID)!ImposedTangentialFacesX
            call Read_Lock(mMAP_, Me%InstanceID)!ImposedTangentialFacesY


            ImposedTangentialFaceU => Me%ImposedTangentialFaces%U
            ImposedTangentialFaceV => Me%ImposedTangentialFaces%V


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if(present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetImposedTangentialFaces

    !--------------------------------------------------------------------------
  
    subroutine GetWetFaces(Map_ID, WetFaceU, WetFaceV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        integer, dimension(:,:,:), pointer          :: WetFaceU, WetFaceV
        integer, optional, intent (OUT)             :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_              

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)

if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mMAP_, Me%InstanceID)!WetFacesX
            call Read_Lock(mMAP_, Me%InstanceID)!WetFacesY


            WetFaceU => Me%WetFaces%U
            WetFaceV => Me%WetFaces%V


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if(present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWetFaces

    !--------------------------------------------------------------------------


    subroutine UngetMap3D(Map_ID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        integer, optional, intent(OUT)              :: STAT
        integer, pointer, dimension(:,:,:)          :: Array

        !External--------------------------------------------------------------
        integer                                     :: ready_   

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)

if1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mMAP_, Me%InstanceID, "UngetMap3D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if(present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetMap3D

    !--------------------------------------------------------------------------

    subroutine GetLandBoundaryFaces3D(Map_ID, LandBoundaryFaces3DU, &
                                      LandBoundaryFaces3DV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: Map_ID
        integer, dimension(:, :, :)  , pointer      :: LandBoundaryFaces3DU, LandBoundaryFaces3DV
        integer, optional, intent (OUT)             :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_   

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (Map_ID, ready_)
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mMAP_, Me%InstanceID)!LandBoundaryFaces3DU
            call Read_Lock(mMAP_, Me%InstanceID)!LandBoundaryFaces3DV

            LandBoundaryFaces3DU => Me%LandBoundaryFaces3D%U
            LandBoundaryFaces3DV => Me%LandBoundaryFaces3D%V


            STAT_ = SUCCESS_

        else
         
            STAT_ = ready_

        end if if1

        if(present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetLandBoundaryFaces3D


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillMap(Map_ID, STAT)  

        !Arguments-------------------------------------------------------------
        integer                         :: Map_ID
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: ready_, nUsers            

        !Local-----------------------------------------------------------------
        integer                         :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Map_ID, ready_) 

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mMAP_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%ObjTime /=0) then
                    nUsers = DeassociateInstance (mTIME_,           Me%ObjTime)
                    if (nUsers == 0) stop 'KillMap - ModuleMap - ERR01'
                endif
                
                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillMap - ModuleMap - ERR03'
                
                nUsers = DeassociateInstance (mGEOMETRY_,       Me%ObjGeometry)
                if (nUsers == 0) stop 'KillMap - ModuleMap - ERR04'
                
                
                deallocate(Me%WaterPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR05'

                deallocate(Me%LandPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR06'

                deallocate(Me%OpenPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR07'

                deallocate(Me%ComputeFaces3D%U, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR08'

                deallocate(Me%ComputeFaces3D%V, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR09'

                deallocate(Me%ComputeFaces3D%W, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR10'

                deallocate(Me%LandBoundaryFaces3D%U, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR11'

                deallocate(Me%LandBoundaryFaces3D%V, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR12'

                deallocate(Me%ImposedNormalFaces%U,      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR13'

                deallocate(Me%ImposedNormalFaces%V,      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR14'

                deallocate(Me%ImposedTangentialFaces%U,  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR15'

                deallocate(Me%ImposedTangentialFaces%V,  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR16'

                deallocate(Me%WetFaces%U,  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR17'

                deallocate(Me%WetFaces%V,  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMap - ModuleMap - ERR18'


                nullify (Me%WaterPoints3D   ) 
                nullify (Me%LandPoints3D    )
                nullify (Me%OpenPoints3D    )
                nullify (Me%ComputeFaces3D%U)
                nullify (Me%ComputeFaces3D%V)
                nullify (Me%ComputeFaces3D%W)
                nullify (Me%LandBoundaryFaces3D%U)
                nullify (Me%LandBoundaryFaces3D%V)


                nullify (Me%ImposedNormalFaces%U)
                nullify (Me%ImposedNormalFaces%V)

                nullify (Me%ImposedTangentialFaces%U)
                nullify (Me%ImposedTangentialFaces%V)

                nullify (Me%WetFaces%U)
                nullify (Me%WetFaces%V)


                !Deallocates Instance
                call DeallocateInstance ()

                Map_ID = 0
                STAT_  = SUCCESS_

            end if
        
        else cd1

            STAT_ = ready_

        endif cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillMap
    
    !----------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Map), pointer           :: AuxMap
        type (T_Map), pointer           :: PreviousMap

        !Updates pointers
        if (Me%InstanceID == FirstObjMap%InstanceID) then
            FirstObjMap     => FirstObjMap%Next
        else
            PreviousMap     => FirstObjMap
            AuxMap          => FirstObjMap%Next
            do while (AuxMap%InstanceID /= Me%InstanceID)
                PreviousMap => AuxMap
                AuxMap      => AuxMap%Next
            enddo

            !Now update linked list
            PreviousMap%Next => AuxMap%Next

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

    subroutine Ready (Map_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                             :: Map_ID
        integer                             :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (Map_ID > 0) then
            call LocateObjMap(Map_ID)
            ready_ = VerifyReadLock (mMAP_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjMap (Map_ID)

        !Arguments-------------------------------------------------------------
        integer                             :: Map_ID

        !Local-----------------------------------------------------------------

        Me => FirstObjMap
        do while (associated (Me))
            if (Me%InstanceID == Map_ID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                        &
            stop 'ModuleMap - LocateObjMap - ERR01'

    end subroutine LocateObjMap
    !--------------------------------------------------------------------------

end module ModuleMap

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
