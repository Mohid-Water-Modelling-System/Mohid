!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Triangulation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to perform Delaunay Triangulation 
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

Module ModuleTriangulation

    use ModuleGlobalData

    implicit none 

    private

    !Constructor
    public  :: ConstructTriangulation
    private ::      AllocateInstance
    private ::      VerifyCollinearPoints
    private ::      ConstructDirectedEdges
    private ::          AlreadyInList
    private ::          AddEdge
    private ::      ConstructNodes
    private ::      RemoveDuplicateNodes
    private ::      ConstructTriangles


    !Modifier
    public  :: SetHeightValues
    public  :: InterPolation
    public  :: CalculateReaches
    private ::      AddNewReach
    private ::      InsertReachToList
    private ::      RemoveDepressions
    public  :: CalculateVoronoi

    !Selector
    public  :: GetNumberOfNodes
    public  :: GetNodesList
    public  :: GetNumberOfTriangles
    public  :: GetTriangleList
    public  :: GetNumberOfReaches
    public  :: GetReaches
    public  :: GetNumberOfBoundaryNodes
    public  :: GetBoundaryNodes

    !Destructor
    public  :: KillTriangulation
    private ::      DeallocateInstance

    !Auxilary
    private :: TRMESH
    private ::      ADDNOD
    private ::      trfind
    private ::      left
    private ::      crtri
    private ::      indxcc
    private ::      bdyadd
    private ::      intadd 
    private ::      lstptr
    private ::      swptst
    private ::      swap
    private ::      jrand
    private ::      store
    private ::      insert
    private ::      trlist
    private ::      circum
    private ::      nbcnt
    private ::      bnodes
    private ::      nearnd

    !Management
    private ::      Ready


    interface ConstructTriangulation
        module procedure ConstructTriangulationXY
        module procedure ConstructTriangulationXYZ
    end interface ConstructTriangulation


    !Parameter                              
    integer, parameter                              :: pnRow        = 9  !See Robert Renka Code

    integer, parameter                              :: Pit          = 1
    integer, parameter                              :: Unflooded    = 2
    integer, parameter                              :: Flooded      = 3
    integer, parameter                              :: CurrentLake  = 4

    !Types---------------------------------------------------------------------

    type T_DirectedEdge
        integer                                     :: EdgeID
        integer                                     :: StartNode
        integer                                     :: EndNode
        real                                        :: VoronoiX
        real                                        :: VoronoiY
        real                                        :: Length
        real                                        :: Slope
        logical                                     :: Boundary
        type (T_DirectedEdge), pointer              :: CounterClockEdge
    end type T_DirectedEdge

    type T_Node
        real                                        :: X
        real                                        :: Y
        real                                        :: Z
        integer                                     :: nNeighbor
        logical                                     :: Boundary
        integer                                     :: State
        real                                        :: VoronoiArea
        type (T_DirectedEdge), pointer              :: FirstEdge
        type (T_DirectedEdge), pointer              :: FlowExit
    end type T_Node

    type T_Triangles
        real                                        :: CenterX, CenterY
        real                                        :: Radius, Area
        real                                        :: AspectRatio
    end type T_Triangles

    type T_Reach
        type (T_DirectedEdge), pointer              :: Edge
        integer                                     :: nStrahler            = null_int
        real                                        :: DrainageArea         = null_real
        type (T_Reach), pointer                     :: Next
    end type T_Reach

    type       T_Triangulation
        !Instance ID
        integer                                     :: InstanceID
        real,    dimension(:), pointer              :: XT, YT, ZT
        integer, dimension(:), pointer              :: List, Lptr
        integer, dimension(:), pointer              :: Lend
        integer                                     :: Lnew
        integer, dimension(:), pointer              :: Near, NextTri
        integer, dimension(:, :), pointer           :: Ltri
        integer, dimension(:), pointer              :: BNodes
        real, dimension(:), pointer                 :: Dist
        integer                                     :: NumberOfNodes
        integer                                     :: NumberOfTriangles
        integer                                     :: NumberOfArcs
        integer                                     :: NumberOfBoundaryNodes
        logical                                     :: MustSwap = .false.
        integer                                     :: SwapNode
        real                                        :: MinX  = -null_real
        real                                        :: MaxX  =  null_real
        real                                        :: MinY  = -null_real
        real                                        :: MaxY  =  null_real
        type (T_Node), dimension(:), pointer        :: Nodes
        type (T_DirectedEdge), dimension(:), pointer :: DirectedEdges
        type (T_Triangles), dimension(:), pointer   :: Triangles
        logical                                     :: HaveHeightValues     = .false.
        type (T_Reach), pointer                     :: FirstReach
        integer                                     :: nReaches             = 0
        type (T_Triangulation), pointer             :: Next
    end type T_Triangulation

    !Global Variables
    type (T_Triangulation), pointer                 :: FirstTriangulation
    type (T_Triangulation), pointer                 :: Me

    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !
    !This routine Constructs a Delaunay Triangulation from given X, and Y points
    !Before the routine Interpolation can be called, the routine SetHeigth values 
    !must be called. This routine is designed to serve the Mohid model to 
    !interpolate Gauge levels at the open boundary
    !
    subroutine ConstructTriangulationXY(TriangulationID,                       &
                                        NumberOfNodes, NodeX, NodeY, WriteTriangles, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer                                     :: NumberOfNodes
        real, dimension(:), pointer                 :: NodeX, NodeY
        logical, optional, intent(IN )              :: WriteTriangles
        integer, optional, intent(OUT)              :: STAT     

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_, STAT_         
        real                                        :: AuxX, AuxY
        integer, dimension(1)                       :: lcc, lct

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTriangulation_)) then
            nullify (FirstTriangulation)
            call RegisterModule (mTriangulation_) 
        endif

        call Ready(TriangulationID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance
       
            !Verifies if the Number of nodes is greater then 2
            if (NumberOfNodes < 3) then
                write(*,*)'Insuficient number of triangulation nodes'
                stop 'ConstructTriangulation - ModuleTriangulation - ERR03'
            endif

            !Stores the number of triangulation nodes
            Me%NumberOfNodes = NumberOfNodes

            !Allocates variables of ObjTriangulation
            allocate(Me%XT(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%YT(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%ZT(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%List(1:6*NumberOfNodes-12), STAT = STAT_CALL)
            allocate(Me%Lptr(1:6*NumberOfNodes-12), STAT = STAT_CALL)
            allocate(Me%Lend(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%Near(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%NextTri(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%Dist(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%Ltri(pnRow, 2*NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%BNodes(1:NumberOfNodes), STAT = STAT_CALL)

            !Initializes XT, YT, ZT
            Me%XT(:) = NodeX(:)
            Me%YT(:) = NodeY(:)
            Me%ZT(:) = null_real

            call VerifyCollinearPoints

            if (Me%MustSwap) then
                AuxX               = Me%XT(3)
                AuxY               = Me%YT(3)
                Me%XT(3)           = Me%XT(Me%SwapNode)
                Me%YT(3)           = Me%YT(Me%SwapNode)
                Me%XT(Me%SwapNode) = AuxX
                Me%YT(Me%SwapNode) = AuxY
            endif

            call TRMESH (Me%NumberOfNodes, Me%XT,            &
                         Me%YT, Me%List,                     &
                         Me%Lptr, Me%lend,                   &
                         Me%lnew, Me%near,                   &
                         Me%NextTri, Me%Dist, STAT_CALL)
            if (STAT_CALL /= 0) stop 'ConstructTriangulation - ModuleTriangulation - ERR13'


            if(present(WriteTriangles))then
                if(WriteTriangles)then
                    
                    call TRLIST (0, lcc, Me%NumberOfNodes, Me%List,  &
                                 Me%Lptr, Me%lend, 9,                &
                                 Me%NumberOfTriangles, Me%Ltri,      &
                                 lct, STAT_CALL)
                    if (STAT_CALL /= 0) stop 'ConstructTriangulation - ModuleTriangulation - ERR14'

                end if
            end if


            !Gets information about Boundary node, etc
            call Bnodes (Me%NumberOfNodes, Me%List,          &
                         Me%Lptr, Me%Lend,                   &
                         Me%BNodes, Me%NumberOfBoundaryNodes,&
                         Me%NumberOfArcs, Me%NumberOfTriangles)


            !Returns ID
            TriangulationID    = Me%InstanceID

        else

            stop 'ModuleTriangulation - ConstructTriangulationXY - ERR99' 

        endif

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructTriangulationXY

    !
    !This routine Constructs a Delaunay Triangulation from given X, Y and Z points.
    !At the same time it fills the T_NODES, T_DIRECTEDEDGES and the T_Triangulation
    !structures, which are used by the Basin module

    subroutine ConstructTriangulationXYZ(TriangulationID,                       &
                                         NumberOfNodes, NodeX, NodeY, NodeZ,    &
                                         Tolerance, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer                                     :: NumberOfNodes
        real, dimension(NumberOfNodes)              :: NodeX, NodeY, NodeZ
        real                                        :: Tolerance
        integer, optional, intent(OUT)              :: STAT     

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_, STAT_         
        integer, dimension(1)                       :: lcc, lct
        real                                        :: AuxX, AuxY, AuxZ

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTriangulation_)) then
            nullify (FirstTriangulation)
            call RegisterModule (mTriangulation_) 
        endif

        call Ready(TriangulationID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance
       
            !Verifies if the Number of nodes is greater then 2
            if (NumberOfNodes < 3) then
                write(*,*)'Insuficient number of triangulation nodes'
                stop 'ConstructTriangulation - ModuleTriangulation - ERR03'
            endif

            call RemoveDuplicateNodes (NumberOfNodes, NodeX, NodeY, Tolerance)

            !Stores the number of triangulation nodes
            Me%NumberOfNodes = NumberOfNodes

            !Allocates variables of ObjTriangulation
            allocate(Me%XT(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%YT(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%ZT(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%List(1:6*NumberOfNodes-12), STAT = STAT_CALL)
            allocate(Me%Lptr(1:6*NumberOfNodes-12), STAT = STAT_CALL)
            allocate(Me%Lend(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%Near(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%NextTri(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%Dist(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%BNodes(1:NumberOfNodes), STAT = STAT_CALL)
            allocate(Me%Ltri(pnRow, 2*NumberOfNodes), STAT = STAT_CALL)

            Me%XT(1:NumberOfNodes) = NodeX(1:NumberOfNodes)
            Me%YT(1:NumberOfNodes) = NodeY(1:NumberOfNodes)
            Me%ZT(1:NumberOfNodes) = NodeZ(1:NumberOfNodes)

            call VerifyCollinearPoints

            if (Me%MustSwap) then
                AuxX               = Me%XT(3)
                AuxY               = Me%YT(3)
                AuxZ               = Me%ZT(3)
                Me%XT(3)           = Me%XT(Me%SwapNode)
                Me%YT(3)           = Me%YT(Me%SwapNode)
                Me%ZT(3)           = Me%ZT(Me%SwapNode)
                Me%XT(Me%SwapNode) = AuxX
                Me%YT(Me%SwapNode) = AuxY
                Me%ZT(Me%SwapNode) = AuxZ
            endif

            call TRMESH (Me%NumberOfNodes, Me%XT,            &
                         Me%YT, Me%List,                     &
                         Me%Lptr, Me%lend,                   &
                         Me%lnew, Me%near,                   &
                         Me%NextTri, Me%Dist, STAT_CALL)
            if (STAT_CALL /= 0) stop 'ConstructTriangulation - ModuleTriangulation - ERR13'

            call TRLIST (0, lcc, Me%NumberOfNodes, Me%List,  &
                         Me%Lptr, Me%lend, 9,                &
                         Me%NumberOfTriangles, Me%Ltri,      &
                         lct, STAT_CALL)
            if (STAT_CALL /= 0) stop 'ConstructTriangulation - ModuleTriangulation - ERR14'

            !Gets information about Boundary node, etc
            call Bnodes (Me%NumberOfNodes, Me%List,          &
                         Me%Lptr, Me%Lend,                   &
                         Me%BNodes, Me%NumberOfBoundaryNodes,&
                         Me%NumberOfArcs, Me%NumberOfTriangles)

            !Construct Nodes
            call ConstructNodes         (NodeX, NodeY, NodeZ)

            !Constructs Triangles
            call ConstructTriangles

            !Construct Directed Edges
            call ConstructDirectedEdges

            !Returns ID
            TriangulationID    = Me%InstanceID

        else

            stop 'ModuleTriangulation - ConstructTriangulationXYZ - ERR99' 

        endif

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructTriangulationXYZ

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Triangulation), pointer                :: NewTriangulation
        type (T_Triangulation), pointer                :: PreviousTriangulation


        !Allocates new instance
        allocate (NewTriangulation)
        nullify  (NewTriangulation%Next)
        nullify  (NewTriangulation%FirstReach)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstTriangulation)) then
            FirstTriangulation    => NewTriangulation
            Me                    => NewTriangulation
        else
            PreviousTriangulation => FirstTriangulation
            Me                    => FirstTriangulation%Next
            do while (associated(Me))
                PreviousTriangulation  => Me
                Me                     => Me%Next
            enddo
            Me                         => NewTriangulation
            PreviousTriangulation%Next => NewTriangulation
        endif

        Me%InstanceID = RegisterNewInstance (mTRIANGULATION_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine VerifyCollinearPoints 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer                 :: XT, YT
        integer                                     :: SwapNode
        logical                                     :: NodeFound

        XT => Me%XT
        YT => Me%YT
        
        !Checks if first 3 points are collinear
        if ( .not. left(XT(1),YT(1),XT(2),YT(2),XT(3),YT(3)) ) then
            Me%MustSwap = .false.
        else if ( .not. left(XT(2),YT(2),XT(1),YT(1),XT(3),YT(3)) ) then
            Me%MustSwap = .false.
        else
            Me%MustSwap = .true.
        endif
                    

        !Checks for first point which is not collinear with point 1 and 2
        if (Me%MustSwap) then
            SwapNode = 4
            NodeFound = .false.
            do while (.not. NodeFound)
                if      (.not. left(XT(1),YT(1),XT(2),YT(2),XT(SwapNode), YT(SwapNode)) ) then
                    Me%SwapNode = SwapNode
                    NodeFound = .true.
                elseif  (.not. left(XT(2),YT(2),XT(1),YT(1),XT(SwapNode), YT(SwapNode)) ) then
                    Me%SwapNode = SwapNode
                    NodeFound = .true.
                else
                    SwapNode = SwapNode + 1
                    if (SwapNode > Me%NumberOfNodes) then
                        write(*,*)'All Triangulation Nodes are collinear'
                        stop 'VerifyCollinearPoints - ModuleTriangulation - ERR01'
                    endif
                endif
            enddo

        endif


    end subroutine 
    
    !--------------------------------------------------------------------------

    subroutine ConstructDirectedEdges
    
        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        integer                                     :: v1, v2, v3
        integer                                     :: n1, n2, n3
        integer                                     :: e1, e2, e3
        integer                                     :: EdgeSoFar, iT, iE, iE2, iN
        integer                                     :: lpl, lp, nd, k, i
        integer, dimension(1000)                    :: nabor
        integer                                     :: EndNodeCounterClock
        integer                                     :: StartNode, EndNode

        !Allocates Variable                        
        allocate(Me%DirectedEdges(2*Me%NumberOfArcs))

        !Stores Directed Edge Information
        EdgeSoFar = 0
        do iT = 1, Me%NumberOfTriangles

            !Vertex Pointers
            v1 = Me%LTRI(1, iT)
            v2 = Me%LTRI(2, iT)
            v3 = Me%LTRI(3, iT)
            
            !Neighbor Pointers
            n1 = Me%LTRI(4, iT)
            n2 = Me%LTRI(5, iT)
            n3 = Me%LTRI(6, iT)

            !Edge Pointers    
            e1 = Me%LTRI(7, iT)
            e2 = Me%LTRI(8, iT)
            e3 = Me%LTRI(9, iT)

            !Test all three edges of the triangle if they are already in the Directed 
            !Edge list. If not add them

            !Vertex 1
            if (.not. AlreadyInList(EdgeSoFar, v2, v3)) then
                call AddEdge (EdgeSoFar, e1, v2, v3, n1, iT)
            endif

            !Vertex 2
            if (.not. AlreadyInList(EdgeSoFar, v3, v1)) then
                call AddEdge (EdgeSoFar, e2, v3, v1, n2, iT)
            endif

            !Vertex 3
            if (.not. AlreadyInList(EdgeSoFar, v1, v2)) then
                call AddEdge (EdgeSoFar, e3, v1, v2, n3, iT)
            endif

        enddo

        !Updates information about counter clockwise edge
        do iE = 1, 2 * Me%NumberOfArcs

            !Find the first CounterClockwise Neighbor of the Starting Node
            StartNode = Me%DirectedEdges(iE)%StartNode
            EndNode   = Me%DirectedEdges(iE)%EndNode

            !Robert Renka Stuff
            lpl       = Me%Lend(StartNode)
            lp        = lpl
            k         = 0

            do
                k        = k + 1
                lp       = Me%lptr(lp)
                nd       = Me%list(lp)
                nabor(k) = nd
                if (lp == lpl ) then
                    exit
                end if
            end do
            
            !NODE is a boundary node.
            if (nd <= 0) then
                nabor(k) = -nd
            end if
            
            !Find the node index of counterclock edge the list of Boundary Nodes
            do i = 1, k
                if (nabor(i) == EndNode) then
                    if (i == k) then
                        EndNodeCounterClock = nabor(1)
                    else
                        EndNodeCounterClock = nabor(i+1)
                    endif
                    exit
                endif
            enddo

            !Points to the CounterClock Edge
            do iE2 = 1, 2 * Me%NumberOfArcs
                if (Me%DirectedEdges(iE2)%StartNode   == StartNode) then
                    if (Me%DirectedEdges(iE2)%EndNode == EndNodeCounterClock) then
                        Me%DirectedEdges(iE)%CounterClockEdge =>           &
                            Me%DirectedEdges(iE2)
                        exit
                    endif
                endif
            enddo

        enddo

        !Updates PointerMe%s to first edges
        do iN = 1, Me%NumberOfNodes
doEdge:     do iE = 1, 2 * Me%NumberOfArcs
                if (iN == Me%DirectedEdges(iE)%StartNode) then
                    Me%Nodes(iN)%FirstEdge => Me%DirectedEdges(iE)
                    exit doEdge
                endif
            enddo doEdge
        enddo


    end subroutine ConstructDirectedEdges

    !--------------------------------------------------------------------------

    logical function AlreadyInList (EdgeSoFar, StartNode, EndNode)

        !Arguments-------------------------------------------------------------
        integer                                     :: EdgeSoFar, StartNode, EndNode

        !Local-----------------------------------------------------------------
        integer                                     :: iE

        AlreadyInList = .false.
        do iE = EdgeSoFar, 1, -1 
            if (Me%DirectedEdges(iE)%StartNode   == StartNode) then
                if (Me%DirectedEdges(iE)%EndNode == EndNode  ) then
                    AlreadyInList = .true.
                    return
                endif
            endif
        enddo

    end function AlreadyInList

    !--------------------------------------------------------------------------

    subroutine AddEdge (EdgeSoFar, EdgeID, StartNode, EndNode, Nb, iT)

        !Arguments-------------------------------------------------------------
        integer                                     :: EdgeSoFar, EdgeID
        integer                                     :: StartNode, EndNode
        integer                                     :: Nb, iT

        !Local-----------------------------------------------------------------
        real                                        :: dx, dy, dz, Length, Slope
        logical                                     :: BoundaryEdge

        !Length / Slope between nodes
        dx = Me%Nodes(StartNode)%X - Me%Nodes(EndNode)%X
        dy = Me%Nodes(StartNode)%Y - Me%Nodes(EndNode)%Y
        dz = Me%Nodes(StartNode)%Z - Me%Nodes(EndNode)%Z

        Length = sqrt(dx**2. + dy**2.)
        Slope  = dz / Length

        !Verifies if Edge is a Boundary edge
        if (Me%Nodes(StartNode)%Boundary .and.                             &
            Me%Nodes(EndNode)%Boundary .and. Nb == 0) then
            BoundaryEdge = .true.
        else
            BoundaryEdge = .false.
        endif


        !Increase Number of Edges and stores the Edge
        EdgeSoFar = EdgeSoFar + 1

        Me%DirectedEdges(EdgeSoFar)%EdgeID    = EdgeID
        Me%DirectedEdges(EdgeSoFar)%StartNode = StartNode
        Me%DirectedEdges(EdgeSoFar)%EndNode   = EndNode
        Me%DirectedEdges(EdgeSoFar)%Length    = Length
        Me%DirectedEdges(EdgeSoFar)%Slope     = Slope
        Me%DirectedEdges(EdgeSoFar)%Boundary  = BoundaryEdge
        !Voronoi Center of the right hand side 
        if (Nb /= 0) then
            Me%DirectedEdges(EdgeSoFar)%VoronoiX =                         &
                                                    Me%Triangles(Nb)%CenterX
            Me%DirectedEdges(EdgeSoFar)%VoronoiY =                         &
                                                    Me%Triangles(Nb)%CenterY
        else
            Me%DirectedEdges(EdgeSoFar)%VoronoiX  = null_real
            Me%DirectedEdges(EdgeSoFar)%VoronoiY  = null_real
        endif

        !Stores the complementary Edge
        EdgeSoFar = EdgeSoFar + 1
        Me%DirectedEdges(EdgeSoFar)%EdgeID    = EdgeID
        Me%DirectedEdges(EdgeSoFar)%StartNode = EndNode
        Me%DirectedEdges(EdgeSoFar)%EndNode   = StartNode
        Me%DirectedEdges(EdgeSoFar)%Length    = Length
        Me%DirectedEdges(EdgeSoFar)%Slope     = -1. * Slope
        Me%DirectedEdges(EdgeSoFar)%Boundary  = BoundaryEdge

        !Voronoi Center of the right hand side 
        Me%DirectedEdges(EdgeSoFar)%VoronoiX =                         &
                                                Me%Triangles(iT)%CenterX
        Me%DirectedEdges(EdgeSoFar)%VoronoiY =                         &
                                                Me%Triangles(iT)%CenterY

    
    end subroutine AddEdge

    !--------------------------------------------------------------------------

    subroutine ConstructNodes (NodeX, NodeY, NodeZ)

        !Arguments-------------------------------------------------------------
        real, dimension(Me%NumberOfNodes)               :: NodeX, NodeY, NodeZ
            
        !Local-----------------------------------------------------------------
        integer                                         :: iN, iB

        !Allocates Nodes
        allocate(Me%Nodes(Me%NumberOfNodes))

        do iN = 1, Me%NumberOfNodes
            nullify (Me%Nodes(iN)%FlowExit )
            nullify (Me%Nodes(iN)%FirstEdge)
        enddo
        

        Me%HaveHeightValues = .true.
        do iN = 1, Me%NumberOfNodes
            Me%Nodes(iN)%X         = NodeX(iN)
            Me%Nodes(iN)%Y         = NodeY(iN)
            Me%Nodes(iN)%Z         = NodeZ(iN)
            Me%Nodes(iN)%nNeighbor = nbcnt(Me%lend(iN), Me%lptr)
            Me%Nodes(iN)%Boundary  = .false.
            do iB = 1, Me%NumberOfBoundaryNodes
                if (Me%BNodes(iB) == iN) then
                   Me%Nodes(iN)%Boundary = .true.
                   exit
                endif
            enddo
            
            !Checks Min and Maximum
            if (Me%Nodes(iN)%X < Me%MinX) Me%MinX = Me%Nodes(iN)%X
            if (Me%Nodes(iN)%X > Me%MaxX) Me%MaxX = Me%Nodes(iN)%X
            if (Me%Nodes(iN)%Y < Me%MinY) Me%MinY = Me%Nodes(iN)%Y
            if (Me%Nodes(iN)%Y > Me%MaxY) Me%MaxY = Me%Nodes(iN)%Y

        enddo

    end subroutine ConstructNodes

    !--------------------------------------------------------------------------

    subroutine RemoveDuplicateNodes (NumberOfNodes, NodeX, NodeY, Tolerance)

        !Arguments-------------------------------------------------------------
        integer                                     :: NumberOfNodes
        real, dimension(NumberOfNodes)              :: NodeX, NodeY
        real                                        :: Tolerance

        !Local-----------------------------------------------------------------
        integer                                     :: iN, iT
        real                                        :: dist

        !Removes duplicate nodes
        iN = 1
        do while (iN <= NumberOfNodes)
            iT = iN + 1
            do while (iT <= NumberOfNodes)
                dist = sqrt((NodeX(iN)-NodeX(iT))**2.+(NodeY(iN)-NodeY(iT))**2.)
                if (dist <= Tolerance) then
                    NodeX (iT) = NodeX(NumberOfNodes)
                    NodeY (iT) = NodeY(NumberOfNodes)
                    NumberOfNodes = NumberOfNodes - 1
                else
                    iT = iT + 1
                endif
            enddo
            iN = iN + 1
        enddo

    end subroutine RemoveDuplicateNodes
    
    !--------------------------------------------------------------------------

    subroutine ConstructTriangles

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: iT

        !Allocates Triangles
        allocate(Me%Triangles(Me%NumberOfTriangles))

        !Calculates the Area, Aspect Ratio, CircumCenter of each Triangle
        do iT = 1, Me%NumberOfTriangles
            call Circum (x1 = Me%Nodes(Me%Ltri(1, iT))%X,     &
                         y1 = Me%Nodes(Me%Ltri(1, iT))%Y,     &
                         x2 = Me%Nodes(Me%Ltri(2, iT))%X,     &
                         y2 = Me%Nodes(Me%Ltri(2, iT))%Y,     &
                         x3 = Me%Nodes(Me%Ltri(3, iT))%X,     &
                         y3 = Me%Nodes(Me%Ltri(3, iT))%Y,     &
                         ratio = .true.,                                                  &
                         xc = Me%Triangles(iT)%CenterX,                     &
                         yc = Me%Triangles(iT)%CenterY,                     &
                         cr = Me%Triangles(iT)%Radius,                      &
                         sa = Me%Triangles(iT)%Area,                        &
                         ar = Me%Triangles(iT)%AspectRatio)
        enddo
        


    end subroutine ConstructTriangles

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine SetHeightValues (TriangulationID, NodeZ, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        real, dimension(:), pointer                 :: NodeZ
        integer, optional, intent(OUT)              :: STAT     
    
        !Local-----------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_            
        real                                        :: AuxZ

        !----------------------------------------------------------------------                         

        call Ready(TriangulationID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Verifies the number of triangulation nodes
            if (size(NodeZ) /= Me%NumberOfNodes) then 
                write(*,*)'Number of triangulation nodes diferent from the initial number'
                stop 'SetHeightValues - ModuleTriangulation - ERR01'
            endif

            Me%ZT(:)            = NodeZ(:)
            Me%HaveHeightValues = .true.

            if (Me%MustSwap) then
                AuxZ = Me%ZT(3)
                Me%ZT(3) = Me%ZT(Me%SwapNode)
                Me%ZT(Me%SwapNode) = AuxZ
            endif

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine SetHeightValues
    
    !--------------------------------------------------------------------------

    real function InterPolation(TriangulationID, PX, PY, FillOutsidePoints, Default, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        real                                        :: PX, PY
        logical                                     :: FillOutsidePoints
        real, optional                              :: Default
        integer, optional, intent(OUT)              :: STAT     
    
        !Local-----------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_            
        real                                        :: a,b,c,d
        real                                        :: x1, y1, z1
        real                                        :: x2, y2, z2
        real                                        :: x3, y3, z3
        integer                                     :: i1, i2, i3
        real                                        :: dist, distNearest, distNextPoint, distPrevPoint
        integer                                     :: NearestNode, NextBNodeID, PrevBNodeID
        integer                                     :: iN, iNode

        !----------------------------------------------------------------------                         

        call Ready(TriangulationID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (.not. Me%HaveHeightValues) then
                write(*,*)'Height Values not set or not up to date'
                stop 'ModuleTriangulation - InterPolation - ERR01'
            endif


            !Locates the triangle in which the point is located
            call TRFIND (1, PX, PY, Me%NumberOfNodes, Me%XT,    &
                         Me%YT, Me%List, Me%lptr, &
                         Me%lend, i1, i2, i3)
            
            if (i1 == 0 .or. i2 == 0 .or. i3 == 0) then

                !See minimum distance to closed points
                if (i1 /= 0 .and. i2 /= 0 .and. FillOutsidePoints) then

                    !Finds nearest boundary Node
                    distNearest       = -null_real
                    do iN = 1, Me%NumberOfBoundaryNodes
                        iNode = Me%BNodes(iN)
                        dist  = sqrt((Me%XT(iNode)-PX)**2. + (Me%YT(iNode)-PY)**2.)
                        if (dist < distNearest) then
                            distNearest       = dist
                            NearestNode       = iNode
                        endif                        
                    enddo

                    !Finds Prev and next node
                    do iN = 1, Me%NumberOfBoundaryNodes
                        if (Me%BNodes(iN) == NearestNode) then
                            if (iN == Me%NumberOfBoundaryNodes) then
                                NextBNodeID = Me%BNodes(1)
                            else
                                NextBNodeID = Me%BNodes(iN + 1)
                            endif
                            if (iN == 1) then
                                PrevBNodeID = Me%BNodes(Me%NumberOfBoundaryNodes)
                            else
                                PrevBNodeID = Me%BNodes(iN-1)
                            endif
                            exit 
                        endif
                    enddo

                    distPrevPoint = sqrt((Me%XT(PrevBNodeID)-PX)**2.0 + (Me%YT(PrevBNodeID)-PY)**2.0)
                    distNextPoint = sqrt((Me%XT(NextBNodeID)-PX)**2.0 + (Me%YT(NextBNodeID)-PY)**2.0)


                    !Point on the left of Nearest -> Next?
                    if (left(Me%XT(NearestNode), Me%YT(NearestNode), Me%XT(NextBNodeID), Me%YT(NextBNodeID), PX, PY) .and. &
                        left(Me%XT(NearestNode), Me%YT(NearestNode), Me%XT(PrevBNodeID), Me%YT(PrevBNodeID), PX, PY) ) then
                       Interpolation = (Me%ZT(NearestNode) * distPrevPoint +  &
                                         Me%ZT(PrevBNodeID) * distNearest) / &
                                         (distPrevPoint  + distNearest)
                    elseif (.not. left(Me%XT(NearestNode), Me%YT(NearestNode), Me%XT(NextBNodeID),  &
                                       Me%YT(NextBNodeID), PX, PY) .and.                            &
                            .not. left(Me%XT(NearestNode), Me%YT(NearestNode), Me%XT(PrevBNodeID),  &
                                       Me%YT(PrevBNodeID), PX, PY) ) then
                        Interpolation = (Me%ZT(NearestNode) * distNextPoint +     &
                                         Me%ZT(NextBNodeID) * distNearest) / &
                                         (distNextPoint + distNearest)
                    else
                        Interpolation = Me%ZT(NearestNode)
                    endif


                    STAT_ = SUCCESS_

                else

                    if (present(Default)) then
                        InterPolation = Default
                        STAT_ = SUCCESS_
                    else
                        STAT_ = UNKNOWN_
                    endif

                endif
            else

                x1 = Me%XT(i1)
                y1 = Me%YT(i1)
                z1 = Me%ZT(i1)

                x2 = Me%XT(i2)
                y2 = Me%YT(i2)
                z2 = Me%ZT(i2)

                x3 = Me%XT(i3)
                y3 = Me%YT(i3)
                z3 = Me%ZT(i3)


                a =  y1 * (z2-z3) + y2 * (z3-z1) + y3 * (z1-z2)
                b =  z1 * (x2-x3) + z2 * (x3-x1) + z3 * (x1-x2)
                c =  x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)
                d = -x1 * (y2*z3-y3*z2) - x2*(y3*z1-y1*z3) - x3*(y1*z2-y2*z1)

                if (c.ne.0) then
                    InterPolation = - (a * PX + b * PY + d) / c
                    STAT_ = SUCCESS_
                else
                    if (present(Default)) then
                        InterPolation = Default
                        STAT_ = SUCCESS_
                    else
                        STAT_ = UNKNOWN_
                    endif
                endif
            endif

          
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end function InterPolation

    !--------------------------------------------------------------------------

    subroutine CalculateVoronoi (TriangulationID, OutputFile, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        character(len=*), optional                  :: OutputFile
        integer, optional, intent(OUT)              :: STAT     

        !Local-----------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_, STAT_CALL
        integer                                     :: iUnit
        integer                                     :: iN
        real                                        :: Area
        logical                                     :: OutputOn
        integer                                     :: nPoints, iP, iE
        type (T_DirectedEdge), pointer              :: CurrentEdge
        real, dimension(:), allocatable             :: TmpX, TmpY
        logical                                     :: Inside

        call Ready(TriangulationID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Opens the results file
            if (present(OutputFile)) then
                call UnitsManager (iUnit, OPEN_FILE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CalculateVoronoi - ModuleTriangulation - ERR01'

                open (unit=iUnit, status='unknown', file=trim(adjustl(OutputFile)))

                OutputOn = .true.
            else
                OutputOn = .false.
            endif


            do iN = 1, Me%NumberOfNodes

                Area = 0.
                if (.not. Me%Nodes(iN)%Boundary) then
                    
                    allocate (TmpX(Me%Nodes(iN)%nNeighbor+1))
                    allocate (TmpY(Me%Nodes(iN)%nNeighbor+1))
                    nPoints = Me%Nodes(iN)%nNeighbor
                    CurrentEdge => Me%Nodes(iN)%FirstEdge
                    iP = 0
                    do iE = 1, Me%Nodes(iN)%nNeighbor
                        iP = iP + 1
                        TmpX (iP) = CurrentEdge%VoronoiX
                        TmpY (iP) = CurrentEdge%VoronoiY
                        CurrentEdge => CurrentEdge%CounterClockEdge
                    enddo

                    !Close Polygon
                    TmpX(nPoints + 1) = TmpX(1)
                    TmpY(nPoints + 1) = TmpY(1)

                    !Verifies if polygon is inside the domain
                    Inside = .true.
                    if (maxval(TmpX) > Me%MaxX) then
                        Inside = .false.
                    endif
                    if (minval(TmpX) < Me%MinX) then
                        Inside = .false.
                    endif
                    if (maxval(TmpY) > Me%MaxY) then
                        Inside = .false.
                    endif
                    if (minval(TmpY) < Me%MinY) then
                        Inside = .false.
                    endif

                    if (Inside) then

                        !Calculates Area
                        do iP = 1, nPoints
                            Area = Area + TmpX(iP) * TmpY(iP+1) - TmpX(iP+1) * TmpY(iP)
                        enddo    

                        Me%Nodes(iN)%VoronoiArea = 0.5 * Area

                        if (OutputOn) then
                            write(iUnit, *)'<begin_polygon>'
                            do iP = 1, nPoints + 1
                                write(iUnit,*)TmpX(iP), TmpY(iP)
                            enddo
                            write(iUnit, *)'<end_polygon>'
                        endif

                    else

                        Me%Nodes(iN)%VoronoiArea = null_real

                    endif


                    deallocate (TmpX)
                    deallocate (TmpY)

                endif

            enddo

            if (OutputOn) then
                call UnitsManager (iUnit, CLOSE_FILE, STAT = STAT_CALL)
            endif

            STAT_ = SUCCESS_
        else               

            STAT_ = ready_

        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_


    end subroutine CalculateVoronoi

    !--------------------------------------------------------------------------

    subroutine CalculateReaches (TriangulationID, RemovePits, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        logical                                     :: RemovePits
        integer, optional, intent(OUT)              :: STAT     
    
        !Local-----------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_            
        integer                                     :: iN, iE
        type (T_Reach), pointer                     :: CurrentReach, UpStreamReach
        integer                                     :: MinOrder, MaxOrder
        type (T_DirectedEdge), pointer              :: CurrentEdge
        real                                        :: MaiorSlope

        !----------------------------------------------------------------------                         

        call Ready(TriangulationID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (.not. Me%HaveHeightValues) then
                write(*,*)'Height Values not set or not up to date'
                stop 'ModuleTriangulation - CalculateReaches - ERR01'
            endif


            !For each triangle calculates the edge where the reach is located. Adds the 
            !reach to the reach list
            do iN = 1, Me%NumberOfNodes

                if (.not. Me%Nodes(iN)%Boundary) then

                    !Looks for the Flow Exit Edge
                    CurrentEdge => Me%Nodes(iN)%FirstEdge
                    MaiorSlope  = 0.0
                    do iE = 1, Me%Nodes(iN)%nNeighbor
                        if (CurrentEdge%Slope > MaiorSlope .and. .not. CurrentEdge%Boundary) then
                            Me%Nodes(iN)%FlowExit => CurrentEdge
                            MaiorSlope = CurrentEdge%Slope
                        endif
                        CurrentEdge => CurrentEdge%CounterClockEdge
                    enddo

                    !If the current node does not have Flow exit set is as pit
                    if (associated(Me%Nodes(iN)%FlowExit)) then
                        Me%Nodes(iN)%State = UnFlooded
                    else
                        Me%Nodes(iN)%State = Pit
                    endif
                
                endif

            enddo


            !Adds all directed edge which are flow exits to the reaches list
            do iN = 1, Me%NumberOfNodes
                if (associated(Me%Nodes(iN)%FlowExit)) then
                    call AddNewReach(Me%Nodes(iN)%FlowExit)
                endif
            enddo

            if (RemovePits) then
                call RemoveDepressions
            endif


            !Stream Classification after Strahler
            CurrentReach => Me%FirstReach
            do while (associated(CurrentReach))
                
                MinOrder = - null_int
                MaxOrder =   null_int

                !Looks for upstream reaches
                UpStreamReach => Me%FirstReach
                do while (associated(UpStreamReach)) 
                    if (UpStreamReach%Edge%EndNode == CurrentReach%Edge%StartNode) then
                        if (UpStreamReach%nStrahler < MinOrder) then
                            MinOrder = UpStreamReach%nStrahler
                        endif
                        if (UpStreamReach%nStrahler > MaxOrder) then
                            MaxOrder = UpStreamReach%nStrahler
                        endif
                    endif
                    UpStreamReach => UpStreamReach%Next
                enddo
                
                !If Min Order = -null_int and MaxOrder == null_int there is no upstream reach
                !and the classification will be 1
                if (MinOrder == - null_int .and. MaxOrder == null_int) then
                    CurrentReach%nStrahler = 1
                
                !If Min Order <> Max Order the order will be equal MaxOrder
                elseif (MinOrder /= MaxOrder) then
                    CurrentReach%nStrahler = MaxOrder
                
                !MinOrder == MaxOrder -> nStrahler = MaxOrder + 1
                else
                    CurrentReach%nStrahler = MaxOrder + 1
                endif

                CurrentReach => CurrentReach%Next
            enddo
            

            !Upstream Drainage Area (1). Set all Reaches with nStrahler = 1 to zero
            CurrentReach => Me%FirstReach
            do while (associated(CurrentReach))
                if (CurrentReach%nStrahler == 1) then
                    CurrentReach%DrainageArea = 0.0
                endif
                CurrentReach => CurrentReach%Next
            enddo
                
            !Upstream Drainage Area (2). Calcualtes Upstream Areas for other reaches
            CurrentReach => Me%FirstReach
            do while (associated(CurrentReach))
                if (CurrentReach%nStrahler /= 1) then

                    !Looks for upstream reaches
                    UpStreamReach => Me%FirstReach
                    do while (associated(UpStreamReach)) 
                        if (UpStreamReach%Edge%EndNode == CurrentReach%Edge%StartNode) then
                        endif
                        UpStreamReach => UpStreamReach%Next
                    enddo
                endif
                CurrentReach => CurrentReach%Next
            enddo


            STAT_ = SUCCESS_

        else               

            STAT_ = ready_

        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

    end subroutine CalculateReaches

    !--------------------------------------------------------------------------

    subroutine AddNewReach (DirectedEdge)

        !Arguments-------------------------------------------------------------
        type (T_DirectedEdge), pointer              :: DirectedEdge

        !Local-----------------------------------------------------------------
        type (T_Reach), pointer                     :: NewReach
        
        !If no reach with the given characteristic was found, add a new Reach to the list
        allocate (NewReach)
        nullify  (NewReach%Next)

        NewReach%Edge => DirectedEdge

        call InsertReachToList (Me%FirstReach, NewReach,                   &
                                Me%nReaches)

    end subroutine AddNewReach 

    !--------------------------------------------------------------------------

    subroutine InsertReachToList(FirstReach, NewReach, nReaches)

        !Arguments-------------------------------------------------------------
        type (T_Reach), pointer                     :: FirstReach
        type (T_Reach), pointer                     :: NewReach
        integer                                     :: nReaches

        !Local-----------------------------------------------------------------
        type (T_Reach), pointer                     :: CurrentReach
        type (T_Reach), pointer                     :: PreviousReach


        !Inserts a new Reach to the list if Reachs
        if (.not. associated(FirstReach)) then
            FirstReach => NewReach
        else
            PreviousReach => FirstReach
            CurrentReach  => PreviousReach%Next
            do while (associated(CurrentReach))
                PreviousReach => CurrentReach
                CurrentReach  => PreviousReach%Next
            enddo
            PreviousReach%Next => NewReach
        endif

        !Increments Reach number by one and sets ID
        nReaches          = nReaches + 1 

    end subroutine InsertReachToList

    !--------------------------------------------------------------------------

    subroutine RemoveDepressions

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: iN, nN
        logical                                     :: FoundOutlet
        integer                                     :: LowestNode
        integer                                     :: iNeig
        real                                        :: LowestElev
        type (T_DirectedEdge), pointer              :: CurrEdge
        type (T_DirectedEdge), pointer              :: FlowExitEdge
        logical                                     :: Impossible
        logical                                     :: ThereArePits
        integer                                     :: nPits, nPrevPits

        !Applies a lake fill algorithm to find outlets of pit nodes
        !See Appendix A - An object oriented framework for distributed hydrologic
        !and geomorphic modeling using triangulated irregular networks
        !Tucker, Lancaster
        nPrevPits = -1
        ThereArePits = .true.
        do while (nPrevPits /= nPits) 
        write(*,*)'nPits',nPits
        ThereArePits = .false.
        nPrevPits = nPits
        nPits = 0
        do iN = 1, Me%NumberOfNodes
            
            if (Me%Nodes(iN)%State == Pit) then

                ThereArePits = .true.
                nPits = nPits + 1


                FoundOutlet = .false.
                nN = iN
                do while (.not. FoundOutlet)

                    !Start by finding the lowest neighbor that does not drain to to node 
                    !itself
                    LowestNode = iN
                    LowestElev = -null_real
                    CurrEdge => Me%Nodes(nN)%FirstEdge
                    do iNeig = 1, Me%Nodes(nN)%nNeighbor
                        nN = CurrEdge%EndNode 
                        if (Me%Nodes(nN)%Z <= LowestElev) then
                            if (associated(Me%Nodes(nN)%FlowExit)) then
                                if (Me%Nodes(nN)%FlowExit%EndNode /= nN) then
                                    LowestElev = Me%Nodes(nN)%Z
                                    LowestNode = nN
                                endif
                            endif
                        endif
                        CurrEdge => CurrEdge%CounterClockEdge
                    enddo

                    !Find the first node in the drainage network of this node which is lower 
                    !then the current
                    Impossible = .false.
                    if (associated(Me%Nodes(LowestNode)%FlowExit)) then
                        nN = Me%Nodes(LowestNode)%FlowExit%EndNode
                        do while (Me%Nodes(nN)%Z >= Me%Nodes(iN)%Z)
                            if (associated(Me%Nodes(nN)%FlowExit)) then
                                nN = Me%Nodes(nN)%FlowExit%EndNode
                            else
                                Impossible = .true.
                                exit
                            endif
                        enddo
                        if (.not. Impossible) then
                            allocate(FlowExitEdge)
                            FlowExitEdge%StartNode = iN
                            FlowExitEdge%EndNode   = LowestNode
                            call AddNewReach (FlowExitEdge)
                            Me%Nodes(iN)%State = UnFlooded
                            Me%Nodes(iN)%FlowExit => FlowExitEdge
                        else
!                            write (*,*)'Unable to solve Drainage for Node:', iN
!                            write (*,*)'X', Me%Nodes(iN)%X
!                            write (*,*)'Y', Me%Nodes(iN)%Y
!                            write (*,*)'Z', Me%Nodes(iN)%Z
                        endif
                    else
!                            write (*,*)'Unable to solve Drainage for Node:', iN
!                            write (*,*)'X', Me%Nodes(iN)%X
!                            write (*,*)'Y', Me%Nodes(iN)%Y
!                            write (*,*)'Z', Me%Nodes(iN)%Z
                    endif
                    FoundOutlet = .true.

!                    if (Me%Nodes(nN)%Z < Me%Nodes(iN)%Z) then
!                        FoundOutlet = .true.
!                        allocate(FlowExitEdge)
!                        FlowExitEdge%StartNode = iN
!                        FlowExitEdge%EndNode   = nN
!                        call AddNewReach (ObjTriangulation, FlowExitEdge)
!                    endif
                enddo
            endif
        enddo
        enddo             
       
    end subroutine RemoveDepressions

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetNumberOfNodes (TriangulationID, nNodes, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer                                     :: nNodes
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_            

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nNodes = Me%NumberOfNodes

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetNumberOfNodes

    !--------------------------------------------------------------------------

    subroutine GetNodesList (TriangulationID, xt, yt, zt, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        real, dimension(:), pointer                 :: xt, yt, zt
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_
        integer                                     :: iN            

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            do iN = 1, Me%NumberOfNodes
                XT(iN) = Me%XT(iN)
                YT(iN) = Me%YT(iN)
                ZT(iN) = Me%ZT(iN)
            enddo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetNodesList

    !--------------------------------------------------------------------------

    subroutine GetNumberOfTriangles (TriangulationID, nTriangles, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer                                     :: nTriangles
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_            

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nTriangles = Me%NumberOfTriangles

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

    end subroutine GetNumberOfTriangles

    !--------------------------------------------------------------------------

    subroutine GetTriangleList (TriangulationID, v1, v2, v3, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer, dimension(:), pointer              :: v1, v2, v3
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_
        integer                                     :: iT            

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            do iT = 1, Me%NumberOfTriangles
                v1(iT) = Me%LTRI(1, iT)
                v2(iT) = Me%LTRI(2, iT)
                v3(iT) = Me%LTRI(3, iT)
            enddo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

        !------------------

    end subroutine GetTriangleList

    !--------------------------------------------------------------------------

    subroutine GetNumberOfReaches (TriangulationID, nReaches, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer                                     :: nReaches
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_            

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nReaches = Me%nReaches


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

    end subroutine GetNumberOfReaches

    !--------------------------------------------------------------------------

    subroutine GetReaches (TriangulationID, StartNode, EndNode, Strahler, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer, dimension(:), pointer              :: StartNode, EndNode, Strahler
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_
        integer                                     :: iR
        type (T_Reach), pointer                     :: CurrentReach

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            iR = 0
            CurrentReach => Me%FirstReach
CurReach:   do while (associated(CurrentReach))
                iR = iR + 1
                StartNode(iR) = CurrentReach%Edge%StartNode
                EndNode(iR)   = CurrentReach%Edge%EndNode
                Strahler(iR)  = CurrentReach%nStrahler
                CurrentReach => CurrentReach%Next
            enddo CurReach


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

        !------------------

    end subroutine GetReaches

        !--------------------------------------------------------------------------

    subroutine GetNumberOfBoundaryNodes (TriangulationID, NumberOfBoundaryNodes, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer                                     :: NumberOfBoundaryNodes
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_            

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            NumberOfBoundaryNodes = Me%NumberOfBoundaryNodes


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

    end subroutine GetNumberOfBoundaryNodes

    !--------------------------------------------------------------------------

    subroutine GetBoundaryNodes (TriangulationID, Bnodes, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer, dimension(:), pointer              :: Bnodes
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, iN             
        integer                                     :: STAT_


        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            do iN = 1, Me%NumberOfBoundaryNodes
                Bnodes(iN) = Me%Bnodes(iN)
            enddo


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

        !------------------

    end subroutine GetBoundaryNodes


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillTriangulation(TriangulationID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_            
        integer                                     :: nUsers

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TriangulationID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTRIANGULATION_,  Me%InstanceID)

            if (nUsers == 0) then

                !Dealloctes variables of ObjTriangulation
                if (associated  (Me%XT))       deallocate(Me%XT)
                if (associated  (Me%YT))       deallocate(Me%YT)
                if (associated  (Me%ZT))       deallocate(Me%ZT)

                if (associated  (Me%List))     deallocate(Me%List)
                if (associated  (Me%Lptr))     deallocate(Me%Lptr)
                if (associated  (Me%Lend))     deallocate(Me%Lend)
                if (associated  (Me%Near))     deallocate(Me%Near)
                if (associated  (Me%NextTri))  deallocate(Me%NextTri)
                if (associated  (Me%Ltri))     deallocate(Me%Ltri)
                if (associated  (Me%BNodes))   deallocate(Me%BNodes)
                if (associated  (Me%Dist))     deallocate(Me%Dist)

                if (associated  (Me%Nodes))                                 &
                     deallocate(Me%Nodes)
                if (associated  (Me%DirectedEdges))                         &
                     deallocate(Me%DirectedEdges)
                if (associated  (Me%Triangles))                             &
                     deallocate(Me%Triangles)

                !Deallocates Instance
                call DeallocateInstance

                TriangulationID = 0

                STAT_ = SUCCESS_

            end if

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillTriangulation

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Triangulation), pointer             :: AuxTriangulation
        type (T_Triangulation), pointer             :: PreviousTriangulation

        !Updates pointers
        if (Me%InstanceID == FirstTriangulation%InstanceID) then
            FirstTriangulation => FirstTriangulation%Next
        else
            PreviousTriangulation => FirstTriangulation
            AuxTriangulation      => FirstTriangulation%Next
            do while (AuxTriangulation%InstanceID /= Me%InstanceID)
                PreviousTriangulation => AuxTriangulation
                AuxTriangulation      => AuxTriangulation%Next
            enddo

            !Now update linked list
            PreviousTriangulation%Next => AuxTriangulation%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 
            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine trmesh ( n, x, y, list, lptr, lend, lnew, near, next, dist, ier )

        !*******************************************************************************
        !
        !! TRMESH triangulates a set of points in the plane.
        !
        !
        !  Discussion:
        !
        !    This subroutine creates a Delaunay triangulation of a
        !    set of N arbitrarily distributed points in the plane re-
        !    ferred to as nodes.  The Delaunay triangulation is defined
        !    as a set of triangles with the following five properties:
        !
        !    1)  The triangle vertices are nodes.
        !    2)  No triangle contains a node other than its vertices.
        !    3)  The interiors of the triangles are pairwise disjoint.
        !    4)  The union of triangles is the convex hull of the set
        !        of nodes (the smallest convex set which contains
        !        the nodes).
        !    5)  The interior of the circumcircle of each triangle
        !        contains no node.
        !
        !    The first four properties define a triangulation, and the
        !    last property results in a triangulation which is as close
        !    as possible to equiangular in a certain sense and which is
        !    uniquely defined unless four or more nodes lie on a common
        !    circle.  This property makes the triangulation well-suited
        !    for solving closest point problems and for triangle-based
        !    interpolation.
        !
        !    The triangulation can be generalized to a constrained
        !    Delaunay triangulation by a call to Subroutine ADDCST.
        !    This allows for user-specified boundaries defining a non-
        !    convex and/or multiply connected region.
        !
        !    The algorithm for constructing the triangulation has
        !    expected time complexity O(N*log(N)) for most nodal dis-
        !    tributions.  Also, since the algorithm proceeds by adding
        !    nodes incrementally, the triangulation may be updated with
        !    the addition (or deletion) of a node very efficiently.
        !    The adjacency information representing the triangulation
        !    is stored as a linked list requiring approximately 13N
        !    storage locations.
        !
        !
        !    The following is a list of the software package modules
        !    which a user may wish to call directly:
        !
        !    ADDCST - Generalizes the Delaunay triangulation to allow
        !             for user-specified constraints.
        !
        !    ADDNOD - Updates the triangulation by appending or
        !             inserting a new node.
        !
        !    AREAP  - Computes the area bounded by a closed polygonal
        !             curve such as the boundary of the triangula-
        !             tion or of a constraint region.
        !
        !    BNODES - Returns an array containing the indexes of the
        !             boundary nodes in counterclockwise order.
        !             Counts of boundary nodes, triangles, and arcs
        !             are also returned.
        !
        !    CIRCUM - Computes the area, circumcenter, circumradius,
        !             and, optionally, the aspect ratio of a trian-
        !             gle defined by user-specified vertices.
        !
        !    DELARC - Deletes a boundary arc from the triangulation.
        !
        !    DELNOD - Updates the triangulation with the deletion of a
        !             node.
        !
        !    EDGE   - Forces a pair of nodes to be connected by an arc
        !             in the triangulation.
        !
        !    GETNP  - Determines the ordered sequence of L closest
        !             nodes to a given node, along with the associ-
        !             ated distances.  The distance between nodes is
        !             taken to be the length of the shortest connec-
        !             ting path which intersects no constraint
        !             region.
        !
        !    INTSEC - Determines whether or not an arbitrary pair of
        !             line segments share a common point.
        !
        !    JRAND  - Generates a uniformly distributed pseudo-random
        !             integer.
        !
        !    LEFT   - Locates a point relative to a line.
        !
        !    NEARND - Returns the index of the nearest node to an
        !             arbitrary point, along with its squared
        !             distance.
        !
        !    STORE  - Forces a value to be stored in main memory so
        !             that the precision of floating point numbers
        !             in memory locations rather than registers is
        !             computed.
        !
        !    TRLIST - Converts the triangulation data structure to a
        !             triangle list more suitable for use in a fin-
        !             ite element code.
        !
        !    TRLPRT - Prints the triangle list created by TRLIST.
        !
        !    TRMESH - Creates a Delaunay triangulation of a set of nodes.
        !
        !    TRMSHR - Creates a Delaunay triangulation (more effici-
        !             ently than TRMESH) of a set of nodes lying at
        !             the vertices of a (possibly skewed) rectangu-
        !             lar grid.
        !
        !    TRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
        !             file containing a triangulation plot.
        !
        !    TRPRNT - Prints the triangulation data structure and,
        !             optionally, the nodal coordinates.
        !
        !  Author:
        !
        !    Robert Renka,
        !    Department of Computer Science,
        !    University of North Texas,
        !    renka@cs.unt.edu
        !
        !  Parameters:
        !
        !    Input, integer N, the number of nodes in the triangulation.  N >= 3.
        !
        !    Input, real X(N), Y(N), the coordinates of the nodes.  (X(K),Y(K)) is 
        !    referred to as node K, and K is referred to as a nodal index.  The 
        !    first three nodes must not be collinear.
        !
        !    Output, integer LIST(6*N-12), nodal indexes which, along with LPTR,
        !    LEND, and LNEW, define the triangulation as a set of N adjacency 
        !    lists; counterclockwise-ordered sequences of neighboring nodes such
        !    that the first and last neighbors of a boundary node are boundary 
        !    nodes (the first neighbor of an interior node is arbitrary).  In
        !    order to distinguish between interior and boundary nodes, the last 
        !    neighbor of each boundary node is represented by the negative
        !    of its index.
        !
        !    Output, integer LPTR(6*N-12), pointers (LIST indexes) in one-to-one
        !    correspondence with the elements of LIST.  LIST(LPTR(I)) indexes the 
        !    node which follows LIST(I) in cyclical counterclockwise order
        !    (the first neighbor follows the last neighbor).
        !
        !    Output, integer LEND(N), pointers to adjacency lists.  LEND(K)
        !    points to the last neighbor of node K for K = 1,...,N.  Thus, 
        !    LIST(LEND(K)) < 0 if and only if K is a boundary node.
        !
        !    Output, integer LNEW, pointer to the first empty location in LIST
        !    and LPTR (list length plus one).  LIST, LPTR, LEND, and LNEW are 
        !    not altered if IER < 0, and are incomplete if IER > 0.
        !
        !    Workspace NEAR(N), NEXT(N), DIST(N).  The space is used to efficiently
        !    determine the nearest triangulation node to each unprocessed node for 
        !    use by ADDNOD.
        !
        !    Output, integer IER = Error indicator:
        !     0 if no errors were encountered.
        !    -1 if N < 3 on input.
        !    -2 if the first three nodes are collinear.
        !    -4 if an error flag was returned by a call to SWAP in ADDNOD.  This is 
        !      an internal error and should be reported to the programmer.
        !     L if nodes L and M coincide for some M > L.  The linked list represents
        !      a triangulation of nodes 1 to M-1 in this case.
        !
        !  Local parameters:
        !
        !    D =        Squared distance from node K to node I
        !    D1,D2,D3 = Squared distances from node K to nodes 1, 2,
        !               and 3, respectively
        !    EPS =      Half the machine precision
        !    I,J =      Nodal indexes
        !    I0 =       Index of the node preceding I in a sequence of
        !               unprocessed nodes:  I = NEXT(I0)
        !    K =        Index of node to be added and DO-loop index: K > 3
        !    KM1 =      K-1
        !    LCC(1) =   Dummy array
        !    LP =       LIST index (pointer) of a neighbor of K
        !    LPL =      Pointer to the last neighbor of K
        !    NCC =      Number of constraint curves
        !    NEXTI =    NEXT(I)
        !    NN =       Local copy of N
        !    SWTOL =    Tolerance for function SWPTST
        !
        integer n

        real d
        real d1
        real d2
        real d3
        real dist(n)
        real eps
        integer i
        integer i0
        integer ier
        integer j
        integer k
        integer km1
        integer lcc(1)
        integer lend(n)
        integer list(*)
        integer lnew
        integer lp
        integer lpl
        integer lptr(*)
        integer ncc
        integer near(n)
        integer next(n)
        integer nexti
        integer nn
        real swtol
        real x(n)
        real y(n)

        nn = n

        if (nn < 3) then
            ier = -1
            return
        end if
!
!  Compute a tolerance for function SWPTST:  SWTOL = 10*
!  (machine precision)
!
        eps = epsilon ( eps )

        swtol = eps * 20.0
!
!  Store the first triangle in the linked list.
!
        if ( .not. left(x(1),y(1),x(2),y(2),x(3),y(3)) ) then
!
!  The initial triangle is (3,2,1) = (2,1,3) = (1,3,2).
!
            list(1) = 3
            lptr(1) = 2
            list(2) = -2
            lptr(2) = 1
            lend(1) = 2

            list(3) = 1
            lptr(3) = 4
            list(4) = -3
            lptr(4) = 3
            lend(2) = 4

            list(5) = 2
            lptr(5) = 6
            list(6) = -1
            lptr(6) = 5
            lend(3) = 6

        else if ( .not. left(x(2),y(2),x(1),y(1),x(3),y(3)) ) then
!
!  The initial triangle is (1,2,3).
!
            list(1) = 2
            lptr(1) = 2
            list(2) = -3
            lptr(2) = 1
            lend(1) = 2

            list(3) = 3
            lptr(3) = 4
            list(4) = -1
            lptr(4) = 3
            lend(2) = 4

            list(5) = 1
            lptr(5) = 6
            list(6) = -2
            lptr(6) = 5
            lend(3) = 6

        else
!
!  The first three nodes are collinear.
!
            ier = -2
            return
        end if
!
!  Initialize LNEW and test for N = 3.
!
        lnew = 7
        if (nn == 3) then
            ier = 0
            return
        end if
!
!  A nearest-node data structure (NEAR, NEXT, and DIST) is
!  used to obtain an expected-time (N*log(N)) incremental
!  algorithm by enabling constant search time for locating
!  each new node in the triangulation.
!
!  For each unprocessed node K, NEAR(K) is the index of the
!  triangulation node closest to K (used as the starting
!  point for the search in Subroutine TRFIND) and DIST(K)
!  is an increasing function of the distance between nodes
!  K and NEAR(K).
!
!  Since it is necessary to efficiently find the subset of
!  unprocessed nodes associated with each triangulation
!  node J (those that have J as their NEAR entries), the
!  subsets are stored in NEAR and NEXT as follows:  for
!  each node J in the triangulation, I = NEAR(J) is the
!  first unprocessed node in J's set (with I = 0 if the
!  set is empty), L = NEXT(I) (if I > 0) is the second,
!  NEXT(L) (if L > 0) is the third, etc.  The nodes in each
!  set are initially ordered by increasing indexes (which
!  maximizes efficiency) but that ordering is not main-
!  tained as the data structure is updated.
!
!  Initialize the data structure for the single triangle.
!
        near(1) = 0
        near(2) = 0
        near(3) = 0

        do k = nn, 4, -1

            d1 = (x(k)-x(1))**2 + (y(k)-y(1))**2
            d2 = (x(k)-x(2))**2 + (y(k)-y(2))**2
            d3 = (x(k)-x(3))**2 + (y(k)-y(3))**2

            if (d1 <= d2  .and.  d1 <= d3) then
                near(k) = 1
                dist(k) = d1
                next(k) = near(1)
                near(1) = k
            else if (d2 <= d1  .and.  d2 <= d3) then
                near(k) = 2
                dist(k) = d2
                next(k) = near(2)
                near(2) = k
            else
                near(k) = 3
                dist(k) = d3
                next(k) = near(3)
                near(3) = k
            end if

        end do
!
!  Add the remaining nodes.  Parameters for ADDNOD are as follows:
!
!  K = Index of the node to be added.
!  NEAR(K) = Index of the starting node for the search in TRFIND.
!  NCC = Number of constraint curves.
!  LCC = Dummy array (since NCC = 0).
!  KM1 = Number of nodes in the triangulation.
!
        ncc = 0

        do 7 k = 4,nn

            km1 = k-1
            call addnod (k,x(k),y(k),near(k),ncc, lcc,km1,x,y, &
                         list,lptr,lend,lnew, ier)

            if (ier /= 0) return
!
!  Remove K from the set of unprocessed nodes associated with NEAR(K).
!
            i = near(k)

            if (near(i) == k) then

                near(i) = next(k)

            else

                i = near(i)

                do

                    i0 = i
                    i = next(i0)
                    if (i == k) then
                        exit
                    end if

                end do

                next(i0) = next(k)

             end if

            near(k) = 0
!
!  Loop on neighbors J of node K.
!
            lpl = lend(k)
            lp = lpl

    4       continue

            lp = lptr(lp)
            j  = abs(list(lp))
!
!  Loop on elements I in the sequence of unprocessed nodes
!  associated with J:  K is a candidate for replacing J
!  as the nearest triangulation node to I.  The next value
!  of I in the sequence, NEXT(I), must be saved before I
!  is moved because it is altered by adding I to K's set.
!
            i = near(j)

    5       continue

            if (i == 0) go to 6
            nexti = next(i)
!
!  Test for the distance from I to K less than the distance
!  from I to J.
!
            d = (x(k)-x(i))**2. + (y(k)-y(i))**2.
            if (d < dist(i)) then
!
!  Replace J by K as the nearest triangulation node to I:
!  update NEAR(I) and DIST(I), and remove I from J's set
!  of unprocessed nodes and add it to K's set.
!
                    near(i) = k
                    dist(i) = d
                    if (i == near(j)) then
                        near(j) = nexti
                    else
                        next(i0) = nexti
                    end if
                    next(i) = near(k)
                    near(k) = i
            else
                i0 = i
            end if
            !
!  Bottom of loop on I.
!
            i = nexti
            go to 5
!
!  Bottom of loop on neighbors J.
!
    6       continue

            if ( lp /= lpl ) then
                go to 4
            end if

    7   continue

        return

    end subroutine TRMESH


    !--------------------------------------------------------------------------

    subroutine addnod ( k, xk, yk, ist, ncc, lcc, n, x, y, list, lptr, lend, lnew, &
                        ier )
        !
        !
        !! ADDNOD adds a node to a triangulation.
        !
        !
        !  Discussion:
        !
        !    Given a triangulation of N nodes in the plane created by
        !    subroutine TRMESH or TRMSHR, this subroutine updates the
        !    data structure with the addition of a new node in position
        !    K.  If node K is inserted into X and Y (K <= N) rather
        !    than appended (K = N+1), then a corresponding insertion
        !    must be performed in any additional arrays associated
        !    with the nodes.  For example, an array of data values Z
        !    must be shifted down to open up position K for the new
        !    value:  set Z(I+1) to Z(I) for I = N,N-1,...,K.  For
        !    optimal efficiency, new nodes should be appended whenever
        !    possible.  Insertion is necessary, however, to add a non-
        !    constraint node when constraints are present (refer to
        !    subroutine ADDCST).
        !
        !    Note that a constraint node cannot be added by this
        !    routine.  In order to insert a constraint node, it is
        !    necessary to add the node with no constraints present
        !    (call this routine with NCC = 0), update LCC by increment-
        !    ing the appropriate entries, and then create (or restore)
        !    the constraints by a call to ADDCST.
        !
        !    The algorithm consists of the following steps:  node K
        !    is located relative to the triangulation (TRFIND), its
        !    index is added to the data structure (INTADD or BDYADD),
        !    and a sequence of swaps (SWPTST and SWAP) are applied to
        !    the arcs opposite K so that all arcs incident on node K
        !    and opposite node K (excluding constraint arcs) are local-
        !    ly optimal (satisfy the circumcircle test).  Thus, if a
        !    (constrained) Delaunay triangulation is input, a (con-
        !    strained) Delaunay triangulation will result.  All indexes
        !    are incremented as necessary for an insertion.
        !
        !  Modified:
        !
        !    29 March 2002
        !
        !  Author:
        !
        !    Robert Renka,
        !    Department of Computer Science,
        !    University of North Texas,
        !    renka@cs.unt.edu
        !
        !  Parameters:
        !
        !    Input, integer K, the nodal index (index for X, Y, and LEND) of the
        !    new node to be added.  1 <= K <= LCC(1).  (K <= N+1 if NCC=0).
        !
        !    Input, real XK, YK, the coordinates of the new node (to be stored in 
        !    X(K) and Y(K)).  The node must not lie in a constraint region.
        !
        !    Input, integer IST, the index of a node at which TRFIND begins the
        !    search.  Search time depends on the proximity
        !    of this node to node K.  1 <= IST <= N.
        !
        !    Input, integer NCC, the number of constraint curves.  NCC >= 0.
        !
        !    Input/output, integer LCC(*), list of constraint curve starting indexes 
        !    (or dummy array of length 1 if NCC = 0).  Refer to subroutine ADDCST.
        !    On output, starting indexes incremented by 1 to reflect the insertion 
        !    of K unless NCC = 0 or (IER /= 0 and IER /= -4).
        !
        !    Input/output, integer N, the number of nodes in the triangulation.  
        !    N >= 3.  Note that N will be incremented following the addition of node K.
        !
        !    Input, real X(N+1), real Y(N+1), containing the coordinates of the 
        !    nodes in the first N positions with non-constraint nodes
        !    in the first LCC(1)-1 locations if NCC > 0.  On output, updated with
        !    the insertion of XK and YK in the K-th positions (node I+1 was node 
        !    I before the insertion for I = K to N if K <= N)
        !    unless IER /= 0 and IER /= -4.
        !
        !    Input/output, integer LIST(*), LPTR(*), LEND(N), LNEW, the data 
        !    structure associated with the triangulation of nodes 1 to N.  The 
        !    arrays must have sufficient length for N+1 nodes.  Refer to TRMESH.
        !    On output, updated with the addition of node K unless
        !    IER /= 0 and IER /= -4.
        !
        !    Output, integer IER = Error indicator:
        !     0 if no errors were encountered.
        !    -1 if K, IST, NCC, N, or an LCC entry is outside its valid range on input.
        !    -2 if all nodes (including K) are collinear.
        !     L if nodes L and K coincide for some L.
        !    -3 if K lies in a constraint region.
        !    -4 if an error flag is returned by SWAP implying that the triangulation
        !      (geometry) was bad on input.
        !
        integer i
        integer i1
        integer i2
        integer i3
        integer ibk
        integer ier
        integer in1
        integer io1
        integer io2
        integer ist
        integer k
        integer kk
        integer l
        integer lcc(*)
        integer lccip1
        integer lend(*)
        integer list(*)
        integer lnew
        integer lp
        integer lpf
        integer lpo1
        integer lptr(*)
        integer n
        integer ncc
        integer nm1
        real x(*)
        real xk
        real y(*)
        real yk
!
        kk = k
!
!  Test for an invalid input parameter.
!
        if ( kk < 1  .or.  ist < 1  .or.  ist > n &
            .or.  ncc < 0  .or.  n < 3 ) then
            ier = -1
            return
        end if

        lccip1 = n + 1

        do i = ncc, 1, -1
            if ( lccip1-lcc(i) < 3 ) then
                ier = -1
                return
            end if
            lccip1 = lcc(i)
        end do

        if ( kk > lccip1 ) then
            ier = -1
            return
        end if
!
!  Find a triangle (I1,I2,I3) containing K or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed from node K.
!
        call trfind ( ist, xk, yk, n, x, y, list, lptr, lend, i1, i2, i3 )
!
!  Test for collinear nodes, duplicate nodes, and K lying in
!  a constraint region.
!
        if ( i1 == 0 ) then
            ier = -2
            return
        end if

        if ( i3 /= 0 ) then

            l = i1
            if ( xk == x(l)  .and.  yk == y(l) ) then
                ier = l
                return
            end if

            l = i2
            if ( xk == x(l)  .and.  yk == y(l) ) then
                ier = l
                return
            end if

            l = i3
            if ( xk == x(l)  .and.  yk == y(l) ) then
                ier = l
                return
            end if

            if ( ncc > 0  .and.  crtri(ncc,lcc,i1,i2,i3) ) then
                ier = -3
                return
            end if

        else
!
!  K is outside the convex hull of the nodes and lies in a
!  constraint region iff an exterior constraint curve is present.
!
            if ( ncc > 0  .and.  indxcc(ncc,lcc,n,list,lend) /= 0 ) then
                ier = -3
                return
            end if

        end if
!
!  No errors encountered.
!
        ier = 0
        nm1 = n
        n = n + 1

        if (kk < n) then
!
!  Open a slot for K in X, Y, and LEND, and increment all
!  nodal indexes which are greater than or equal to K.
!
!  Note that LIST, LPTR, and LNEW are not yet updated with
!  either the neighbors of K or the edges terminating on K.
!
            do ibk = nm1, kk, -1
                x(ibk+1) = x(ibk)
                y(ibk+1) = y(ibk)
                lend(ibk+1) = lend(ibk)
            end do

            do i = 1, ncc
                lcc(i) = lcc(i) + 1
            end do

            l = lnew - 1

            do i = 1, l

                if ( list(i) >= kk ) then
                    list(i) = list(i) + 1
                end if

                if ( list(i) <= -kk ) then
                    list(i) = list(i) - 1
                end if

            end do

            if ( i1 >= kk ) then
                i1 = i1 + 1
            end if

            if ( i2 >= kk ) then
                i2 = i2 + 1
            end if

            if ( i3 >= kk ) then
                i3 = i3 + 1
            end if

        end if
!
!  Insert K into X and Y, and update LIST, LPTR, LEND, and
!  LNEW with the arcs containing node K.
!
        x(kk) = xk
        y(kk) = yk

        if ( i3 == 0 ) then
            call bdyadd ( kk, i1, i2, list, lptr, lend, lnew )
        else
            call intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )
        end if
!
!  Initialize variables for optimization of the triangulation.
!
        lp = lend(kk)
        lpf = lptr(lp)
        io2 = list(lpf)
        lpo1 = lptr(lpf)
        io1 = abs ( list(lpo1) )
!
!  Begin loop:  find the node opposite K.
!
        do

            lp = lstptr ( lend(io1), io2, list, lptr )

            if ( list(lp) < 0 ) then
                go to 6
            end if

            lp = lptr(lp)
            in1 = abs ( list(lp) )

            if ( crtri ( ncc, lcc, io1, io2, in1 ) ) then
                go to 6
            end if
!
!  Swap test:  if a swap occurs, two new arcs are
!  opposite K and must be tested.
!
            if ( .not. swptst(in1,kk,io1,io2,x,y) ) then
                go to 6
            end if

            call swap ( in1, kk, io1, io2, list, lptr, lend, lpo1 )

            if ( lpo1 == 0 ) then
                ier = -4
                exit
            end if

            io1 = in1

            cycle
!
!  No swap occurred.  Test for termination and reset IO2 and IO1.
!
        6   continue

            if ( lpo1 == lpf  .or.  list(lpo1) < 0 ) then
                exit
            end if

            io2 = io1
            lpo1 = lptr(lpo1)
            io1 = abs(list(lpo1))

        end do

        return
    
    end subroutine addnod

    logical function frwrd(xa,ya,xb,yb,xc,yc) 
    
        !Arguments-------------------------------------------------------------
        real xa
        real xb
        real xc
        real ya
        real yb
        real yc

        if ((xb-xa)*(xc-xa) + (yb-ya)*(yc-ya) >= 0.0) then
            frwrd = .true.
        else
            frwrd = .false.
        endif


    end function frwrd


    subroutine trfind ( nst, px, py, n, x, y, list, lptr, lend, i1, i2, i3 )
        !
        !*******************************************************************************
        !
        !! TRFIND locates a point relative to a triangulation.
        !
        !
        !  Discussion:
        !
        !    This subroutine locates a point P relative to a triangu-
        !    lation created by subroutine TRMESH or TRMSHR.  If P is
        !    contained in a triangle, the three vertex indexes are
        !    returned.  Otherwise, the indexes of the rightmost and
        !    leftmost visible boundary nodes are returned.
        !
        !  Author:
        !
        !    Robert Renka,
        !    Department of Computer Science,
        !    University of North Texas,
        !    renka@cs.unt.edu
        !
        !  Parameters:
        !
        !    Input, integer NST, the index of a node at which TRFIND begins the
        !    search.  Search time depends on the proximity of this node to P.
        !
        !    Input, real PX, PY, the coordinates of the point P to be located.
        !
        !    Input, integer N, the number of nodes in the triangulation.  N >= 3.
        !
        !    Input, real X(N), Y(N), the coordinates of the nodes in the triangulation.
        !
        !    Input, integer LIST(*), LPTR(*), LEND(N), the data structure defining 
        !    the triangulation.  Refer to subroutine TRMESH.
        !
        !    Output, integer I1, I2, I3, nodal indexes, in counterclockwise order,
        !    of the vertices of a triangle containing P if P is contained in a 
        !    triangle.  If P is not in the convex hull of the nodes, I1 indexes 
        !    the rightmost visible boundary node, I2 indexes the leftmost visible
        !    boundary node, and I3 = 0.  Rightmost and leftmost are defined from 
        !    the perspective of P, and a pair of points are visible from each 
        !    other if and only if the line segment joining them intersects no 
        !    triangulation arc.  If P and all of the nodes lie on a common line, 
        !    then I1 = I2 = I3 = 0 on output.
        !
        !  Local parameters:
        !
        !    B1,B2 =    Unnormalized barycentric coordinates of P with respect 
        !               to (N1,N2,N3)
        !    IX,IY,IZ = Integer seeds for JRAND
        !    LP =       LIST pointer
        !    N0,N1,N2 = Nodes in counterclockwise order defining a
        !               cone (with vertex N0) containing P
        !    N1S,N2S =  Saved values of N1 and N2
        !    N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
        !    NB =       Index of a boundary node -- first neighbor of
        !               NF or last neighbor of NL in the boundary traversal loops
        !    NF,NL =    First and last neighbors of N0, or first
        !               (rightmost) and last (leftmost) nodes
        !               visible from P when P is exterior to the triangulation
        !    NP,NPP =   Indexes of boundary nodes used in the boundary traversal loops
        !    XA,XB,XC = Dummy arguments for FRWRD
        !    YA,YB,YC = Dummy arguments for FRWRD
        !    XP,YP =    Local variables containing the components of P

!
        integer n
!
  real b1
  real b2
!  logical frwrd
  integer i1
  integer i2
  integer i3
  integer, save :: ix = 1
  integer, save :: iy = 2
  integer, save :: iz = 3
  integer lend(n)
  integer list(*)
  integer lp
  integer lptr(*)
  integer n0
  integer n1
  integer n1s
  integer n2
  integer n2s
  integer n3
  integer n4
  integer nb
  integer nf
  integer nl
  integer np
  integer npp
  integer nst
  real px
  real py
  real x(n)
!  real xa
!  real xb
!  real xc
  real xp
  real y(n)
!  real ya
!  real yb
!  real yc
  real yp
!
!  Statement function:
!
!  FRWRD = TRUE iff C is forward of A->B iff <A->B,A->C> >= 0.
!
!  frwrd(xa,ya,xb,yb,xc,yc) = (xb-xa)*(xc-xa) + (yb-ya)*(yc-ya) >= 0.0
!
!  Initialize variables.
!
  xp = px
  yp = py
  n0 = nst

  if ( n0 < 1  .or.  n0 > n ) then
    n0 = jrand ( n, ix, iy, iz )
  end if
!
!  Set NF and NL to the first and last neighbors of N0, and
!  initialize N1 = NF.
!
1 continue

  lp = lend(n0)
  nl = list(lp)
  lp = lptr(lp)
  nf = list(lp)
  n1 = nf
!
!  Find a pair of adjacent neighbors N1,N2 of N0 that define
!  a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
!
  if (nl > 0) go to 2
!
!   N0 is a boundary node.  Test for P exterior.
!
  nl = -nl

  if ( .not. left(x(n0),y(n0),x(nf),y(nf),xp,yp) ) then
    nl = n0
    go to 9
  end if

  if ( .not. left(x(nl),y(nl),x(n0),y(n0),xp,yp) ) then
    nb = nf
    nf = n0
    np = nl
    npp = n0
    go to 11
  end if

  go to 3
!
!  N0 is an interior node.  Find N1.
!
2 continue

    do

      if ( left(x(n0),y(n0),x(n1),y(n1),xp,yp) ) then
        exit
      end if

      lp = lptr(lp)
      n1 = list(lp)

      if ( n1 == nl ) then
        go to 6
      end if

    end do
!
!  P is to the left of edge N0->N1.  Initialize N2 to the
!  next neighbor of N0.
!
3 continue

    lp = lptr(lp)
    n2 = abs(list(lp))

    if ( .not. left(x(n0),y(n0),x(n2),y(n2),xp,yp) ) then
      go to 7
    end if

    n1 = n2
    if (n1 /= nl) go to 3

  if ( .not. left(x(n0),y(n0),x(nf),y(nf),xp,yp) ) then
    go to 6
  end if

  if (xp == x(n0) .and. yp == y(n0)) go to 5
!
!  P is left of or on edges N0->NB for all neighbors NB of N0.
!  All points are collinear iff P is left of NB->N0 for
!  all neighbors NB of N0.  Search the neighbors of N0.
!  NOTE: N1 = NL and LP points to NL.
!
4   continue

    if ( .not. left(x(n1),y(n1),x(n0),y(n0),xp,yp) ) then
      go to 5
    end if

    lp = lptr(lp)
    n1 = abs(list(lp))

    if ( n1 == nl ) then
      i1 = 0
      i2 = 0
      i3 = 0
      return
    end if

    go to 4
!
!  P is to the right of N1->N0, or P=N0.  Set N0 to N1 and start over.
!
5 continue

  n0 = n1
  go to 1
!
!  P is between edges N0->N1 and N0->NF.
!
6 continue

  n2 = nf
!
!  P is contained in the wedge defined by line segments
!  N0->N1 and N0->N2, where N1 is adjacent to N2.  Set
!  N3 to the node opposite N1->N2, and save N1 and N2 to
!  test for cycling.
!
7 continue

  n3 = n0
  n1s = n1
  n2s = n2
!
!  Top of edge hopping loop.  Test for termination.
!
8 continue

  if ( left(x(n1),y(n1),x(n2),y(n2),xp,yp) ) then
!
!  P LEFT N1->N2 and hence P is in (N1,N2,N3) unless an
!  error resulted from floating point inaccuracy and
!  collinearity.  Compute the unnormalized barycentric
!  coordinates of P with respect to (N1,N2,N3).
!
    b1 = (x(n3)-x(n2))*(yp-y(n2)) - (xp-x(n2))*(y(n3)-y(n2))
    b2 = (x(n1)-x(n3))*(yp-y(n3)) - (xp-x(n3))*(y(n1)-y(n3))

    if ( store ( b1 + 1.0 ) >= 1.0  .and. store ( b2 + 1.0 ) >= 1.0 ) then
      go to 16
    end if
!
!   Restart with N0 randomly selected.
!
    n0 = jrand(n, ix,iy,iz )
    go to 1

  end if
!
!  Set N4 to the neighbor of N2 which follows N1 (node
!  opposite N2->N1) unless N1->N2 is a boundary edge.
!
  lp = lstptr(lend(n2),n1,list,lptr)

  if (list(lp) < 0) then
    nf = n2
    nl = n1
    go to 9
  end if

  lp = lptr(lp)
  n4 = abs(list(lp))
!
!  Select the new edge N1->N2 which intersects the line
!  segment N0-P, and set N3 to the node opposite N1->N2.
!
  if ( left(x(n0),y(n0),x(n4),y(n4),xp,yp) ) then
    n3 = n1
    n1 = n4
    n2s = n2
    if (n1 /= n1s  .and.  n1 /= n0) go to 8
  else
    n3 = n2
    n2 = n4
    n1s = n1
    if (n2 /= n2s  .and.  n2 /= n0) go to 8
  end if
!
!  The starting node N0 or edge N1-N2 was encountered
!  again, implying a cycle (infinite loop).  Restart
!  with N0 randomly selected.
!
  n0 = jrand(n, ix,iy,iz )
  go to 1
!
!  Boundary traversal loops.  NL->NF is a boundary edge and
!  P RIGHT NL->NF.  Save NL and NF.

9 continue

  np = nl
  npp = nf
!
!  Find the first (rightmost) visible boundary node NF.  NB
!  is set to the first neighbor of NF, and NP is the last neighbor.
!
10 continue

  lp = lend(nf)
  lp = lptr(lp)
  nb = list(lp)

  if ( .not. left(x(nf),y(nf),x(nb),y(nb),xp,yp) ) then
    go to 12
  end if
!
!  P LEFT NF->NB and thus NB is not visible unless an error
!  resulted from floating point inaccuracy and collinear-
!  ity of the 4 points NP, NF, NB, and P.
!
11 continue

  if ( frwrd(x(nf),y(nf),x(np),y(np),xp,yp)  .or. &
       frwrd(x(nf),y(nf),x(np),y(np),x(nb),y(nb)) ) then
    i1 = nf
    go to 13
  end if
!
!  Bottom of loop.
!
12 continue

  np = nf
  nf = nb
  go to 10
!
!  Find the last (leftmost) visible boundary node NL.  NB
!  is set to the last neighbor of NL, and NPP is the first
!  neighbor.
!
13 continue

  lp = lend(nl)
  nb = -list(lp)

  if ( .not. left(x(nb),y(nb),x(nl),y(nl),xp,yp) ) then
    go to 14
  end if
!
!  P LEFT NB->NL and thus NB is not visible unless an error
!  resulted from floating point inaccuracy and collinear-
!  ity of the 4 points P, NB, NL, and NPP.
!
  if ( frwrd(x(nl),y(nl),x(npp),y(npp),xp,yp)  .or. &
       frwrd(x(nl),y(nl),x(npp),y(npp),x(nb),y(nb)) ) &
    go to 15
!
!  Bottom of loop.
!
14 continue

  npp = nl
  nl = nb
  go to 13
!
!  NL is the leftmost visible boundary node.
!
15 continue

  i2 = nl
  i3 = 0
  return
!
!  P is in the triangle (N1,N2,N3).
!
   16 i1 = n1
  i2 = n2
  i3 = n3

  return

    end subroutine trfind


logical function left ( x1, y1, x2, y2, x0, y0 )
!
!*******************************************************************************
!
!! LEFT determines whether a node is to the left of a line.
!
!
!  Discussion:
!
!    This function determines whether node N0 is to the left
!    or to the right of the line through N1-N2 as viewed by an
!    observer at N1 facing N2.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, real X1, Y1, coordinates of N1.
!
!    Input, real X2, Y2, coordinates of N2.
!
!    Input, real X0, Y0, coordinates of N0.
!
!    Output, logical LEFT, is .TRUE. if and only if (X0,Y0) is on or 
!    to the left of the directed line N1->N2.
!
!  Local parameters:
!
!    DX1,DY1 = X,Y components of the vector N1->N2
!    DX2,DY2 = X,Y components of the vector N1->N0
!
!
  real dx1
  real dx2
  real dy1
  real dy2
  real x0
  real x1
  real x2
  real y0
  real y1
  real y2
!
  dx1 = x2-x1
  dy1 = y2-y1
  dx2 = x0-x1
  dy2 = y0-y1
!
!  If the sign of the vector cross product of N1->N2 and
!  N1->N0 is positive, then sin(A) > 0, where A is the
!  angle between the vectors, and thus A is in the range
!  (0,180) degrees.
!
  left = dx1*dy2 >= dx2*dy1

  return
end function left

logical function crtri ( ncc, lcc, i1, i2, i3 )
!
!*******************************************************************************
!
!! CRTRI determines if a triangle lies in a constraint region.
!
!
!  Discussion:
!
!    This function returns TRUE if and only if triangle (I1,
!    I2,I3) lies in a constraint region.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, NCC, LCC = Constraint data structure.  Refer to Sub-
!    routine ADDCST.
!
!    Input, I1, I2, I3 = Nodal indexes of the counterclockwise-
!    ordered vertices of a triangle.
!
!    Output, logical CRTRI, is TRUE if and only if (I1,I2,I3) is a 
!    constraint region triangle.
!
  integer i
  integer i1
  integer i2
  integer i3
  integer imax
  integer imin
  integer lcc(*)
  integer ncc
!
  imax = max ( i1, i2, i3 )
!
!  Find the index I of the constraint containing IMAX.
!
  i = ncc + 1

  do

    i = i - 1

    if ( i <= 0 ) then
      crtri = .false.
      return
    end if

    if ( imax >= lcc(i) ) then
      exit
    end if

  end do

  imin = min ( i1, i2, i3 )
!
!  P lies in a constraint region iff I1, I2, and I3 are nodes
!  of the same constraint (IMIN >= LCC(I)), and (IMIN,IMAX)
!  is (I1,I3), (I2,I1), or (I3,I2).
!
  crtri = imin >= lcc(i)  .and.  ((imin == i1 .and. &
          imax == i3)  .or.  (imin == i2  .and. &
          imax == i1)  .or.  (imin == i3  .and. &
          imax == i2))

  return

end function crtri

integer function indxcc ( ncc, lcc, n, list, lend )
!
!*******************************************************************************
!
!! INDXCC returns the index of an exterior constraint curve.
!
!
!  Discussion:
!
!    Given a constrained Delaunay triangulation, this func-
!    tion returns the index, if any, of an exterior constraint
!    curve (an unbounded constraint region).  An exterior con-
!    straint curve is assumed to be present if and only if the
!    clockwise-ordered sequence of boundary nodes is a subse-
!    quence of a constraint node sequence.  The triangulation
!    adjacencies corresponding to constraint edges may or may
!    not have been forced by a call to ADDCST, and the con-
!    straint region may or may not be valid (contain no nodes).
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer NCC, the number of constraints.  NCC >= 0.
!
!    Input, integer LCC(*), list of constraint curve starting indexes (or
!    dummy array of length 1 if NCC = 0).  Refer to subroutine ADDCST.
!
!    Input, integer N, the number of nodes in the triangulation.  N >= 3.
!
!    Input, integer LIST(*), LEND(N), the data structure defining the 
!    triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer INDXCC, index of the exterior constraint curve, if
!    present, or 0 otherwise.
!
  implicit none
!
  integer n
!
  integer i
  integer ifrst
  integer ilast
  integer lcc(*)
  integer lend(n)
  integer list(*)
  integer lp
  integer n0
  integer ncc
  integer nst
  integer nxt
!
  indxcc = 0

  if ( ncc < 1 ) then
    return
  end if
!
!  Set N0 to the boundary node with smallest index.
!
  n0 = 0

  do

    n0 = n0 + 1
    lp = lend(n0)

    if ( list(lp) <= 0 ) then
      exit
    end if

  end do
!
!  Search in reverse order for the constraint I, if any, that
!  contains N0.  IFRST and ILAST index the first and last
!  nodes in constraint I.
!
  i = ncc
  ilast = n

  do

    ifrst = lcc(i)

    if ( n0 >= ifrst ) then
      exit
    end if

    if ( i == 1 ) then
      return
    end if

    i = i - 1
    ilast = ifrst - 1

  end do
!
!  N0 is in constraint I which indexes an exterior constraint
!  curve iff the clockwise-ordered sequence of boundary
!  node indexes beginning with N0 is increasing and bounded
!  above by ILAST.
!
  nst = n0

  do

    nxt = -list(lp)

    if ( nxt == nst ) then
      exit
    end if

    if ( nxt <= n0  .or.  nxt > ilast ) then
      return
    end if

    n0 = nxt
    lp = lend(n0)

  end do
!
!  Constraint I contains the boundary node sequence as a subset.
!
  indxcc = i

  return

end function INDXCC


subroutine bdyadd ( kk, i1, i2, list, lptr, lend, lnew )
!
!*******************************************************************************
!
!! BDYADD adds a boundary node to a triangulation.
!
!
!  Discussion:
!
!    This subroutine adds a boundary node to a triangulation
!    of a set of points in the plane.  The data structure is
!    updated with the insertion of node KK, but no optimization
!    is performed.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer KK, the index of a node to be connected to the sequence
!    of all visible boundary nodes.  KK >= 1 and
!    KK must not be equal to I1 or I2.
!
!    Input, integer I1, the first (rightmost as viewed from KK) boundary
!    node in the triangulation which is visible from
!    node KK (the line segment KK-I1 intersects no arcs.
!
!    Input, integer I2, the last (leftmost) boundary node which is visible
!    from node KK.  I1 and I2 may be determined by subroutine TRFIND.
!
!    Input/output, integer LIST(*), LPTR(*), LEND(N), LNEW.  The 
!    triangulation data structure created by TRMESH or TRMSHR.
!    On input, nodes I1 and I2 must be included in the triangulation.
!    On output, the data structure has been updated with the addition 
!    of node KK.  Node KK is connected to I1, I2, and all boundary 
!    nodes in between.
!
!
  integer i1
  integer i2
  integer k
  integer kk
  integer lend(*)
  integer list(*)
  integer lnew
  integer lp
  integer lptr(*)
  integer lsav
  integer n1
  integer n2
  integer next
  integer nsav
!
  k = kk
  n1 = i1
  n2 = i2
!
!  Add K as the last neighbor of N1.
!
  lp = lend(n1)
  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = -k
  lptr(lnew) = lsav
  lend(n1) = lnew
  lnew = lnew + 1
  next = -list(lp)
  list(lp) = next
  nsav = next
!
!  Loop on the remaining boundary nodes between N1 and N2,
!  adding K as the first neighbor.
!
  do

    lp = lend(next)

    call insert ( k, lp, list, lptr, lnew )

    if ( next == n2 ) then
      exit
    end if

    next = -list(lp)
    list(lp) = next

  end do
!
!  Add the boundary nodes between N1 and N2 as neighbors
!  of node K.
!
  lsav = lnew
  list(lnew) = n1
  lptr(lnew) = lnew + 1
  lnew = lnew + 1
  next = nsav

  do

    if ( next == n2 ) then
      exit
    end if

    list(lnew) = next
    lptr(lnew) = lnew + 1
    lnew = lnew + 1
    lp = lend(next)
    next = list(lp)

  end do

  list(lnew) = -n2
  lptr(lnew) = lsav
  lend(k) = lnew
  lnew = lnew + 1

  return
end subroutine bdyadd

subroutine intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )
!
!*******************************************************************************
!
!! INTADD adds an interior point to a triangulation.
!
!
!  Discussion:
!
!    This subroutine adds an interior node to a triangulation
!    of a set of points in the plane.  The data structure is
!    updated with the insertion of node KK into the triangle
!    whose vertices are I1, I2, and I3.  No optimization of the
!    triangulation is performed.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer KK, the index of the node to be inserted.  KK >= 1
!    and KK must not be equal to I1, I2, or I3.
!
!    Input, integer I1, I2, I3, indexes of the counterclockwise-ordered
!    sequence of vertices of a triangle which contains node KK.
!
!    Input/output, integer LIST(*), LPTR(*), LEND(N), LNEW, the data 
!    structure defining the triangulation.  Refer to subroutine TRMESH. 
!    Triangle (I1,I2,I3) must be included in the triangulation.
!    On output, updated with the addition of node KK.  KK
!    will be connected to nodes I1, I2, and I3.
!
!
  integer i1
  integer i2
  integer i3
  integer k
  integer kk
  integer lend(*)
  integer list(*)
  integer lnew
  integer lp
  integer lptr(*)
  integer n1
  integer n2
  integer n3
!
  k = kk
!
!  Initialization.
!
  n1 = i1
  n2 = i2
  n3 = i3
!
!  Add K as a neighbor of I1, I2, and I3.
!
  lp = lstptr(lend(n1),n2,list,lptr)
  call insert (k,lp,list,lptr,lnew)
  lp = lstptr(lend(n2),n3,list,lptr)
  call insert (k,lp,list,lptr,lnew)
  lp = lstptr(lend(n3),n1,list,lptr)
  call insert (k,lp,list,lptr,lnew)
!
!  Add I1, I2, and I3 as neighbors of K.
!
  list(lnew) = n1
  list(lnew+1) = n2
  list(lnew+2) = n3
  lptr(lnew) = lnew + 1
  lptr(lnew+1) = lnew + 2
  lptr(lnew+2) = lnew
  lend(k) = lnew + 2
  lnew = lnew + 3

  return
end subroutine intadd

integer function lstptr ( lpl, nb, list, lptr )
!
!*******************************************************************************
!
!! LSTPTR returns the index of NB in the adjacency list for N0.
!
!
!  Discussion:
!
!    This function returns the index (LIST pointer) of NB in
!    the adjacency list for N0, where LPL = LEND(N0).
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer LPL = LEND(N0).
!
!    Input, integer NB, the index of the node whose pointer is to be 
!    returned.  NB must be connected to N0.
!
!    Input, integer LIST(*), LPTR(*), the data structure defining the 
!    triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer LSTPTR, pointer such that LIST(LSTPTR) = NB or
!    LIST(LSTPTR) = -NB, unless NB is not a neighbor of N0, in which 
!    case LSTPTR = LPL.
!
!
  integer list(*)
  integer lp
  integer lpl
  integer lptr(*)
  integer nb
  integer nd
!
  lp = lptr(lpl)

  do

    nd = list(lp)

    if ( nd == nb ) then
      exit
    end if

    lp = lptr(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do

  lstptr = lp

  return

end function lstptr 

logical function swptst ( in1, in2, io1, io2, x, y )
!
!*******************************************************************************
!
!! SWPTST applies the circumcircle test to a quadrilateral.
!
!
!  Discussion:
!
!    This function applies the circumcircle test to a quadri-
!    lateral defined by a pair of adjacent triangles.  The
!    diagonal arc (shared triangle side) should be swapped for
!    the other diagonl if and only if the fourth vertex is
!    strictly interior to the circumcircle of one of the
!    triangles (the decision is independent of the choice of
!    triangle).  Equivalently, the diagonal is chosen to maxi-
!    mize the smallest of the six interior angles over the two
!    pairs of possible triangles (the decision is for no swap
!    if the quadrilateral is not strictly convex).
!
!    When the four vertices are nearly cocircular (the
!    neutral case), the preferred decision is no swap -- in
!    order to avoid unnecessary swaps and, more important, to
!    avoid cycling in subroutine OPTIM which is called by
!    DELNOD and EDGE.  Thus, a tolerance SWTOL (stored in
!    SWPCOM by TRMESH or TRMSHR) is used to define 'nearness'
!    to the neutral case.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer IN1, IN2, IO1, IO2, the nodal indexes of the vertices of
!    the quadrilateral.  IO1-IO2 is the triangulation arc (shared triangle
!    side) to be replaced by IN1-IN2 if the decision is to swap.  The
!    triples (IO1,IO2,IN1) and (IO2,IO1,IN2) must define triangles (be
!    in counterclockwise order) on input.
!
!    Input, real X(*), Y(*), the nodal coordinates.
!
!    Output, logical SWPTST, .TRUE. if and only if the arc connecting
!    IO1 and IO2 is to be replaced.
!
!  Local parameters:
!
!    DX11,DY11 = X,Y components of the vector IN1->IO1
!    DX12,DY12 = X,Y components of the vector IN1->IO2
!    DX22,DY22 = X,Y components of the vector IN2->IO2
!    DX21,DY21 = X,Y components of the vector IN2->IO1
!    SIN1 =      Cross product of the vectors IN1->IO1 and
!                IN1->IO2 -- proportional to sin(T1), where
!                T1 is the angle at IN1 formed by the vectors
!    COS1 =      Inner product of the vectors IN1->IO1 and
!                IN1->IO2 -- proportional to cos(T1)
!    SIN2 =      Cross product of the vectors IN2->IO2 and
!                IN2->IO1 -- proportional to sin(T2), where
!                T2 is the angle at IN2 formed by the vectors
!    COS2 =      Inner product of the vectors IN2->IO2 and
!                IN2->IO1 -- proportional to cos(T2)
!    SIN12 =     SIN1*COS2 + COS1*SIN2 -- proportional to sin(T1+T2)
!
!
  real cos1
  real cos2
  real dx11
  real dx12
  real dx21
  real dx22
  real dy11
  real dy12
  real dy21
  real dy22
  integer in1
  integer in2
  integer io1
  integer io2
  real sin1
  real sin12
  real sin2
  real swtol
  real x(*)
  real y(*)
  real eps
!
!  Tolerance stored by TRMESH or TRMSHR.
!
  eps = epsilon ( eps )
  swtol = eps * 20.0
!
!  Compute the vectors containing the angles T1 and T2.
!
  dx11 = x(io1) - x(in1)
  dx12 = x(io2) - x(in1)
  dx22 = x(io2) - x(in2)
  dx21 = x(io1) - x(in2)

  dy11 = y(io1) - y(in1)
  dy12 = y(io2) - y(in1)
  dy22 = y(io2) - y(in2)
  dy21 = y(io1) - y(in2)
!
!  Compute inner products.
!
  cos1 = dx11*dx12 + dy11*dy12
  cos2 = dx22*dx21 + dy22*dy21
!
!  The diagonals should be swapped iff (T1+T2) > 180
!  degrees.  The following two tests ensure numerical
!  stability:  the decision must be FALSE when both
!  angles are close to 0, and TRUE when both angles
!  are close to 180 degrees.
!
  if ( cos1 >= 0.0  .and.  cos2 >= 0.0 ) then
    swptst = .false.
    return
  end if

  if ( cos1 < 0.0  .and.  cos2 < 0.0 ) then
    swptst = .true.
    return
  end if
!
!  Compute vector cross products (Z-components).
!
  sin1 = dx11*dy12 - dx12*dy11
  sin2 = dx22*dy21 - dx21*dy22
  sin12 = sin1*cos2 + cos1*sin2

  if ( sin12 >= -swtol ) then
    swptst = .false.
  else
    swptst = .true.
  end if

  return
end function swptst


subroutine swap ( in1, in2, io1, io2, list, lptr, lend, lp21 )
!
!*******************************************************************************
!
!! SWAP adjusts a triangulation by swapping a diagonal arc.
!
!
!  Discussion:
!
!    Given a triangulation of a set of points on the unit
!    sphere, this subroutine replaces a diagonal arc in a
!    strictly convex quadrilateral (defined by a pair of adja-
!    cent triangles) with the other diagonal.  Equivalently, a
!    pair of adjacent triangles is replaced by another pair
!    having the same union.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer IN1, IN2, IO1, IO2, the nodal indexes of the vertices of
!    the quadrilateral.  IO1-IO2 is replaced by IN1-IN2.  (IO1,IO2,IN1)
!    and (IO2,IO1,IN2) must be triangles on input.
!
!    Input/output, integer LIST(*), LPTR(*), LEND(N), the data structure
!    defining the triangulation.  Refer to subroutine TRMESH.  On output,
!    updated with the swap; triangles (IO1,IO2,IN1) and (IO2,IO1,IN2) are
!    replaced by (IN1,IN2,IO2) and (IN2,IN1,IO1) unless LP21 = 0.
!
!    Output, integer LP21, the index of IN1 as a neighbor of IN2 after the
!    swap is performed unless IN1 and IN2 are adjacent on input, in which 
!    case LP21 = 0.
!
!  Local parameters:
!
!    LP, LPH, LPSAV = LIST pointers
!
!
  integer in1
  integer in2
  integer io1
  integer io2
  integer lend(*)
  integer list(*)
  integer lp
  integer lp21
  integer lph
  integer lpsav
  integer lptr(*)
!
!  Test for IN1 and IN2 adjacent.
!
  lp = lstptr(lend(in1),in2,list,lptr)

  if (abs(list(lp)) == in2) then
    lp21 = 0
    return
  end if
!
!  Delete IO2 as a neighbor of IO1.
!
  lp = lstptr(lend(io1),in2,list,lptr)
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO2 is the last neighbor of IO1, make IN2 the last neighbor.
!
  if (lend(io1) == lph) then
    lend(io1) = lp
  end if
!
!  Insert IN2 as a neighbor of IN1 following IO1
!  using the hole created above.
!
  lp = lstptr(lend(in1),io1,list,lptr)
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in2
  lptr(lph) = lpsav
!
!  Delete IO1 as a neighbor of IO2.
!
  lp = lstptr(lend(io2),in1,list,lptr)
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO1 is the last neighbor of IO2, make IN1 the last neighbor.
!
  if (lend(io2) == lph) then
    lend(io2) = lp
  end if
!
!  Insert IN1 as a neighbor of IN2 following IO2.
!
  lp = lstptr(lend(in2),io2,list,lptr)
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in1
  lptr(lph) = lpsav
  lp21 = lph

  return
end subroutine swap

integer function jrand ( n, ix, iy, iz )
!
!*******************************************************************************
!
!! JRAND returns a uniformly distributed random integer between 1 and N.
!
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    B. A. Wichmann and I. D. Hill, 
!    An Efficient and Portable Pseudo-random Number Generator,
!    Applied Statistics, 
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input, integer N, the maximum value to be returned.
!
!    Input/output, integer IX, IY, IZ, seeds initialized to values in
!    the range 1 to 30,000 before the first call to JRAND, and not altered 
!    by the user between subsequent calls (unless a sequence of random 
!    numbers is to be repeated by reinitializing the seeds).
!
!    Output, integer JRAND, random integer in the range 1 to N.
!
!  Local parameters:
!
!    U = Pseudo-random number uniformly distributed in the interval (0,1).
!    X = Pseudo-random number in the range 0 to 3 whose fractional part is U.
!
!
  integer ix
  integer iy
  integer iz
  integer n
  real u
  real x
!
  ix = mod ( 171 * ix, 30269 )
  iy = mod ( 172 * iy, 30307 )
  iz = mod ( 170 * iz, 30323 )

  x = ( real ( ix ) / 30269.0E+00 ) &
    + ( real ( iy ) / 30307.0E+00 ) &
    + ( real ( iz ) / 30323.0E+00 )
 
  u = x - int ( x )
  jrand = real ( n ) * u + 1.0E+00

  return
end function jrand

real function store ( x )
!
!*******************************************************************************
!
!! STORE forces its argument to be stored.
!
!
!  Discussion:
!
!    This function forces its argument X to be stored in a
!    memory location, thus providing a means of determining
!    floating point number characteristics (such as the machine
!    precision) when it is necessary to avoid computation in
!    high precision registers.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, real X, the value to be stored.
!
!    Output, real STORE, the value of X after it has been stored and
!    possibly truncated or rounded to the single precision word length.
!
!
  real x
  real y
!
  common /stcom/ y
!
  y = x
  store = y

  return
end function store

subroutine insert ( k, lp, list, lptr, lnew )
!
!*******************************************************************************
!
!! INSERT inserts K as a neighbor of N1.
!
!
!  Discussion:
!
!    This subroutine inserts K as a neighbor of N1 following
!    N2, where LP is the LIST pointer of N2 as a neighbor of
!    N1.  Note that, if N2 is the last neighbor of N1, K will
!    become the first neighbor (even if N1 is a boundary node).
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer K, the index of the node to be inserted.
!
!    Input, integer LP, the LIST pointer of N2 as a neighbor of N1.
!
!    Input/output, integer LIST(*), LPTR(*), LNEW, the data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.  On output,
!    the data structure has been updated to include node K.
!
!
  integer k
  integer list(*)
  integer lnew
  integer lp
  integer lptr(*)
  integer lsav
!
  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = k
  lptr(lnew) = lsav
  lnew = lnew + 1

  return
end subroutine insert 

subroutine trlist ( ncc, lcc, n, list, lptr, lend, nrow, nt, ltri, lct, ier )
!
!*******************************************************************************
!
!! TRLIST converts a triangulation to triangle list form.
!
!
!  Discussion:
!
!    This subroutine converts a triangulation data structure
!    from the linked list created by subroutine TRMESH or
!    TRMSHR to a triangle list.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer NCC, the number of constraints.  NCC >= 0.
!
!    Input, integer LCC(*), list of constraint curve starting indexes (or
!    dummy array of length 1 if NCC = 0).  Refer to subroutine ADDCST.
!
!    Input, integer N, the number of nodes in the triangulation.  N >= 3.
!
!    Input, integer LIST(*), LPTR(*), LEND(N), linked list data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.
!
!    Input, integer NROW, the number of rows (entries per triangle) 
!    reserved for the triangle list LTRI.  The value must be 6 if only 
!    the vertex indexes and neighboring triangle indexes are to be
!    stored, or 9 if arc indexes are also to be assigned and stored.  
!    Refer to LTRI.
!
!    Input, integer LTRI(NROW*NT), where NT is at most 2N-5.  (A sufficient
!    length is 12*N if NROW=6 or 18*N if NROW=9.)
!
!    Output, integer NT, the number of triangles in the triangulation unless
!    IER /= 0, in which case NT = 0.  NT = 2N - NB- 2, where NB is the number 
!    of boundary nodes.
!
!    Output, integer LTRI(NROW,NT), whose J-th column contains the vertex nodal
!    indexes (first three rows), neighboring triangle indexes (second three
!    rows), and, if NROW = 9, arc indexes (last three rows) associated with
!    triangle J for J = 1,...,NT.  The vertices are ordered counterclockwise
!    with the first vertex taken to be the one with smallest index.  Thus,
!    LTRI(2,J) and LTRI(3,J) are larger than LTRI(1,J) and index adjacent
!    neighbors of node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J) and LTRI(I+6,J)
!    index the triangle and arc, respectively, which are opposite (not shared
!    by) node LTRI(I,J), with LTRI(I+3,J) = 0 if LTRI(I+6,J) indexes a boundary
!    arc.  Vertex indexes range from 1 to N, triangle indexes from 0 to NT,
!    and, if included, arc indexes from 1 to NA = NT+N-1.  The triangles are 
!    ordered on first (smallest) vertex indexes, except that the sets of
!    constraint triangle (triangles contained in the closure of a constraint
!    region) follow the non-constraint triangles.
!
!    Output, integer LCT(NCC), containing the triangle index of the first
!    triangle of constraint J in LCT(J).  Thus, the number of non-constraint
!    triangles is LCT(1)-1, and constraint J contains LCT(J+1)-LCT(J) 
!    triangles, where LCT(NCC+1) = NT+1.
!
!    Output, integer IER = Error indicator.
!    0, if no errors were encountered.
!    1, if NCC, N, NROW, or an LCC entry is outside its valid range on input.
!    2, if the triangulation data structure (LIST,LPTR,LEND) is invalid.  
!
!  Local Parameters:
!
!    ARCS = TRUE iff arc indexes are to be stored.
!    KA,KT = Numbers of currently stored arcs and triangles.
!    N1ST = Starting index for the loop on nodes (N1ST = 1 on
!           pass 1, and N1ST = LCC1 on pass 2).
!    NM2 = Upper bound on candidates for N1.
!    PASS2 = TRUE iff constraint triangles are to be stored.
!
  implicit none
!
  integer n
  integer nrow
!
  logical arcs
  logical cstri
  integer i
  integer i1
  integer i2
  integer i3
  integer ier
  integer isv
  integer j
  integer jlast
  integer ka
  integer kn
  integer kt
  integer l
  integer lcc(*)
  integer lcc1
  integer lct(*)
  integer lend(n)
  integer list(*)
  integer lp
  integer lp2
  integer lpl
  integer lpln1
  integer lptr(*)
  integer ltri(nrow,*)
  integer n1
  integer n1st
  integer n2
  integer n3
  integer ncc
  integer nm2
  integer nn
  integer nt
  logical pass2
!
!  Test for invalid input parameters and store the index
!  LCC1 of the first constraint node (if any).
!
  nn = n

  if (ncc < 0  .or.  (nrow /= 6  .and. nrow /= 9)) then
    nt = 0
    ier = 1
    return
  end if

  lcc1 = nn+1

  if (ncc == 0) then

    if ( nn < 3 ) then
      nt = 0
      ier = 1
      return
    end if

  else

    do i = ncc,1,-1
      if (lcc1-lcc(i) < 3) then
        nt = 0
        ier = 1
        return
      end if
      lcc1 = lcc(i)
    end do

    if (lcc1 < 1) then
      nt = 0
      ier = 1
      return
    end if

  end if
!
!  Initialize parameters for loop on triangles KT = (N1,N2,
!  N3), where N1 < N2 and N1 < N3.  This requires two
!  passes through the nodes with all non-constraint
!  triangles stored on the first pass, and the constraint
!  triangles stored on the second.
!
  arcs = nrow == 9
  ka = 0
  kt = 0
  n1st = 1
  nm2 = nn-2
  pass2 = .false.
!
!  Loop on nodes N1:  
!  J = constraint containing N1,
!  JLAST = last node in constraint J.
!
2 continue

  j = 0
  jlast = lcc1 - 1

  do 11 n1 = n1st,nm2

    if (n1 > jlast) then
!
!  N1 is the first node in constraint J+1.  Update J and
!  JLAST, and store the first constraint triangle index
!  if in pass 2.
!
      j = j + 1

      if (j < ncc) then
        jlast = lcc(j+1) - 1
      else
        jlast = nn
      end if

      if (pass2) then
        lct(j) = kt + 1
      end if

    end if
!
!  Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
!  to the last neighbor of N1, and LP2 points to N2.
!
    lpln1 = lend(n1)
    lp2 = lpln1

    3 continue

      lp2 = lptr(lp2)
      n2 = list(lp2)
      lp = lptr(lp2)
      n3 = abs(list(lp))

      if (n2 < n1  .or.  n3 < n1) go to 10
!
!  (N1,N2,N3) is a constraint triangle iff the three nodes
!  are in the same constraint and N2 < N3.  Bypass con-
!  straint triangles on pass 1 and non-constraint triangles
!  on pass 2.
!
      cstri = n1 >= lcc1  .and.  n2 < n3  .and. n3 <= jlast

      if ((cstri  .and.  .not. pass2)  .or. &
          (.not. cstri  .and.  pass2)) go to 10
!
!  Add a new triangle KT = (N1,N2,N3).
!
      kt = kt + 1
      ltri(1,kt) = n1
      ltri(2,kt) = n2
      ltri(3,kt) = n3
!
!  Loop on triangle sides (I1,I2) with neighboring triangles
!  KN = (I1,I2,I3).
!
      do i = 1,3

        if (i == 1) then
          i1 = n3
          i2 = n2
        else if (i == 2) then
          i1 = n1
          i2 = n3
        else
          i1 = n2
          i2 = n1
        end if
!
!  Set I3 to the neighbor of I1 which follows I2 unless
!  I2->I1 is a boundary arc.
!
        lpl = lend(i1)
        lp = lptr(lpl)

4       continue

          if (list(lp) == i2) go to 5
          lp = lptr(lp)
          if (lp /= lpl) go to 4
!
!  I2 is the last neighbor of I1 unless the data structure
!  is invalid.  Bypass the search for a neighboring
!  triangle if I2->I1 is a boundary arc.
!
        if (abs(list(lp)) /= i2) go to 13
        kn = 0
        if (list(lp) < 0) go to 8
!
!  I2->I1 is not a boundary arc, and LP points to I2 as
!  a neighbor of I1.
!
5   continue

        lp = lptr(lp)
        i3 = abs(list(lp))
!
!  Find L such that LTRI(L,KN) = I3 (not used if KN > KT),
!  and permute the vertex indexes of KN so that I1 is
!  smallest.
!
        if (i1 < i2  .and.  i1 < i3) then
          l = 3
        else if (i2 < i3) then
          l = 2
          isv = i1
          i1 = i2
          i2 = i3
          i3 = isv
        else
          l = 1
          isv = i1
          i1 = i3
          i3 = i2
          i2 = isv
        end if
!
!  Test for KN > KT (triangle index not yet assigned).
!
        if (i1 > n1  .and.  .not. pass2) go to 9
! 
!  Find KN, if it exists, by searching the triangle list in
!  reverse order.
!
        do kn = kt-1,1,-1
          if (ltri(1,kn) == i1  .and.  ltri(2,kn) == &
              i2  .and.  ltri(3,kn) == i3) then
            go to 7
          end if
        end do

        go to 9
!
!  Store KT as a neighbor of KN.
!
7       continue

        ltri(l+3,kn) = kt
!
!  Store KN as a neighbor of KT, and add a new arc KA.
!
8       continue

        ltri(i+3,kt) = kn

        if (arcs) then
          ka = ka + 1
          ltri(i+6,kt) = ka
          if (kn /= 0) ltri(l+6,kn) = ka
        end if

9       continue

    end do
! 
!  Bottom of loop on triangles.
!
   10     continue

          if (lp2 /= lpln1) go to 3

   11     continue
!
!  Bottom of loop on nodes.
!
  if (.not. pass2  .and.  ncc > 0) then
    pass2 = .true.
    n1st = lcc1
    go to 2
  end if
!
!  No errors encountered.
!
  nt = kt
  ier = 0
  return
!
!  Invalid triangulation data structure:  I1 is a neighbor of
!  I2, but I2 is not a neighbor of I1.
!
   13 continue

  nt = 0
  ier = 2

  return

end subroutine trlist

subroutine circum ( x1, y1, x2, y2, x3, y3, ratio, xc, yc, cr, sa, ar )
!
!*******************************************************************************
!
!! CIRCUM determines the circumcenter (and more) of a triangle.
!
!
!  Discussion:
!
!    Given three vertices defining a triangle, this subrou-
!    tine returns the circumcenter, circumradius, signed
!    triangle area, and, optionally, the aspect ratio of the
!    triangle.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of the vertices.
!
!    Input, logical RATIO, is TRUE if and only if the aspect ratio is 
!    to be computed.
!
!    Output, real XC, YC, coordinates of the circumcenter (center of the 
!    circle defined by the three points) unless SA = 0, in which XC and YC
!    are not altered.
!
!    Output, real CR, the circumradius (radius of the circle defined by
!    the three points) unless SA = 0 (infinite radius), in which case 
!    CR is not altered.
!
!    Output, real SA, the signed triangle area with positive value if
!    and only if the vertices are specified in counterclockwise order:  
!    (X3,Y3) is strictly to the left of the directed line from (X1,Y1)
!    toward (X2,Y2).
!
!    Output, real AR, the aspect ratio r/CR, where r is the radius of the
!    inscribed circle, unless RATIO = FALSE, in which case AR is not 
!    altered.  AR is in the range 0 to 0.5, with value 0 iff SA = 0 and
!    value 0.5 iff the vertices define an equilateral triangle.
!
  implicit none
!
  real ar
  real cr
  real ds(3)
  real fx
  real fy
  logical ratio
  real sa
  real u(3)
  real v(3)
  real x1
  real x2
  real x3
  real xc
  real y1
  real y2
  real y3
  real yc
!
!  Set U(K) and V(K) to the x and y components, respectively,
!  of the directed edge opposite vertex K.
!
  u(1) = x3 - x2
  u(2) = x1 - x3
  u(3) = x2 - x1
  v(1) = y3 - y2
  v(2) = y1 - y3
  v(3) = y2 - y1
!
!  Set SA to the signed triangle area.
!
  sa = ( u(1) * v(2) - u(2) * v(1) ) / 2.0E+00

  if ( sa == 0.0E+00 ) then
    if ( ratio ) then
      ar = 0.0E+00
    end if
    return
  end if
!
!  Set DS(K) to the squared distance from the origin to vertex K.
!
  ds(1) = x1 * x1 + y1 * y1
  ds(2) = x2 * x2 + y2 * y2
  ds(3) = x3 * x3 + y3 * y3
!
!  Compute factors of XC and YC.
!
  fx = - dot_product ( ds(1:3), v(1:3) )
  fy =   dot_product ( ds(1:3), u(1:3) )

  xc = fx / ( 4.0E+00 * sa )
  yc = fy / ( 4.0E+00 * sa )
  cr = sqrt ( ( xc - x1 )**2 + ( yc - y1 )**2 )

  if ( .not. ratio ) then
    return
  end if
!
!  Compute the squared edge lengths and aspect ratio.
!
  ds(1:3) = u(1:3)**2 + v(1:3)**2

  ar = 2.0E+00 * abs ( sa ) / &
       ( ( sqrt ( ds(1) ) + sqrt ( ds(2) ) + sqrt ( ds(3) ) ) * cr )

  return
end subroutine circum

integer function nbcnt ( lpl, lptr )
!
!*******************************************************************************
!
!! NBCNT returns the number of neighbors of a node.
!
!
!  Discussion:
!
!    This function returns the number of neighbors of a node
!    N0 in a triangulation created by Subroutine TRMESH (or
!    TRMSHR).
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer LPL, the LIST pointer to the last neighbor of N0.
!    LPL = LEND(N0).
!
!    Input, integer LPTR(*), pointers associated with LIST.
!
!    Output, integer NBCNT, the  number of neighbors of N0.
!
  implicit none
!
  integer k
  integer lp
  integer lpl
  integer lptr(*)
!
  lp = lpl
  k  = 1

  do
    lp = lptr(lp)
    if ( lp == lpl ) then
      exit
    end if
    k = k + 1

  end do

  nbcnt = k

  return
end function nbcnt


subroutine bnodes ( n, list, lptr, lend, nodes, nb, na, nt )
!
!*******************************************************************************
!
!! BNODES returns a list of the boundary nodes.
!
!
!  Discussion:
!
!    Given a triangulation of N points in the plane, this
!    subroutine returns an array containing the indexes, in
!    counterclockwise order, of the nodes on the boundary of
!    the convex hull of the set of points.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, integer N, the number of nodes in the triangulation.  N >= 3.
!
!    Input, integer LIST(*), LPTR(*), LEND(N), the data structure defining 
!    the triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer NODES(NB), ordered sequence of boundary node indexes
!    in the range 1 to N.
!
!    Output, integer NB, the number of boundary nodes.
!
!    Output, integer NA, NT, the number of arcs and triangles, respectively,
!    in the triangulation.
!
  implicit none
!
  integer n
!
  integer k
  integer lend(n)
  integer list(*)
  integer lp
  integer lptr(*)
  integer n0
  integer na
  integer nb
  integer nodes(*)
  integer nst
  integer nt
!
!  Set NST to the first boundary node encountered.
!
  nst = 1

  do

    lp = lend(nst)

    if ( list(lp) < 0 ) then
      exit
    end if

    nst = nst + 1

  end do
!
!  Initialization.
!
  nodes(1) = nst
  k = 1
  n0 = nst
!
!  Traverse the boundary in counterclockwise order.
!
  do

    lp = lend(n0)
    lp = lptr(lp)
    n0 = list(lp)

    if ( n0 == nst ) then
      exit
    end if

    k = k + 1
    nodes(k) = n0

  end do
!
!  Termination.
!
  nb = k
  nt = 2*n - nb - 2
  na = nt + n - 1

  return
end subroutine bnodes

function nearnd ( xp, yp, ist, n, x, y, list, lptr, lend, dsq )
!
!*******************************************************************************
!
!! NEARND finds the nearest triangulation node to a point.
!
!
!  Discussion:
!
!    Given a point P in the plane and a Delaunay triangula-
!    tion created by Subroutine TRMESH or TRMSHR, this function
!    returns the index of the nearest triangulation node to P.
!
!    The algorithm consists of implicitly adding P to the
!    triangulation, finding the nearest neighbor to P, and
!    implicitly deleting P from the triangulation.  Thus, it
!    is based on the fact that, if P is a node in a Delaunay
!    triangulation, the nearest node to P is a neighbor of P.
!
!    Note that the number of candidates for NEARND
!    (neighbors of P) is limited to LMAX defined in
!    the PARAMETER statement below.
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Parameters:
!
!    Input, real XP, YP, the coordinates of the point P.
!
!    Input, integer IST, the index of a node at which TRFIND begins 
!    the search.  Search time depends on the proximity
!    of this node to P.
!
!    Input, integer N, the number of nodes in the triangulation.
!
!    Input, real X(N), Y(N), the coordinates of the nodes.
!
!    Input, integer LIST(*), LPTR(*), LEND(N), a data structure 
!    defining the triangulation.  Refer to TRMESH.
!
!    Output, real DSQ, the square of the distance between P and
!    node NEARND.
!
!    Output, integer NEARND, the index of the nearest node to P, 
!    or 0 if N < 3 or the triangulation data structure is invalid.
!
  implicit none
!
  integer, parameter :: lmax = 25
  integer n
!
  real cos1
  real cos2
  real ds1
  real dsq
  real dsr
  real dx11
  real dx12
  real dx21
  real dx22
  real dy11
  real dy12
  real dy21
  real dy22
  integer i1
  integer i2
  integer i3
  integer ist
  integer l
  integer lend(n)
  integer list(*)
  integer listp(lmax)
  integer lp
  integer lp1
  integer lp2
  integer lpl
  integer lptr(*)
  integer lptrp(lmax)
  integer n1
  integer n2
  integer n3
  integer nearnd
  integer nr
  integer nst
  real sin1
  real sin2
  real x(n)
  real xp
  real y(n)
  real yp
!
!  Store local parameters and test for N invalid.
!
  nearnd = 0
  dsq = -1.0E+00

  if ( n < 3 ) then
    return
  end if

  nst = ist

  if ( nst < 1  .or.  nst > n ) then
    nst = 1
  end if
!
!  Find a triangle (I1,I2,I3) containing P, or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed from P.
!
  call trfind (nst,xp,yp,n,x,y,list,lptr,lend, i1,i2,i3)
!
!  Test for collinear nodes.
!
  if ( i1 == 0 ) then
    return
  end if
!
!  Store the linked list of 'neighbors' of P in LISTP and
!  LPTRP.  I1 is the first neighbor, and 0 is stored as
!  the last neighbor if P is not contained in a triangle.
!  L is the length of LISTP and LPTRP, and is limited to LMAX.
!
  if ( i3 /= 0 ) then

    listp(1) = i1
    lptrp(1) = 2
    listp(2) = i2
    lptrp(2) = 3
    listp(3) = i3
    lptrp(3) = 1
    l = 3

  else

    n1 = i1
    l = 1
    lp1 = 2
    listp(l) = n1
    lptrp(l) = lp1
!
!  Loop on the ordered sequence of visible boundary nodes
!  N1 from I1 to I2.
!
    do

      lpl = lend(n1)
      n1 = -list(lpl)
      l = lp1
      lp1 = l+1
      listp(l) = n1
      lptrp(l) = lp1

      if ( n1 == i2  .or.  lp1 >= lmax ) then
        exit
      end if

    end do

    l = lp1
    listp(l) = 0
    lptrp(l) = 1

  end if
!
!  Initialize variables for a loop on arcs N1-N2 opposite P
!  in which new 'neighbors' are 'swapped' in.  N1 follows
!  N2 as a neighbor of P, and LP1 and LP2 are the LISTP
!  indexes of N1 and N2.
!
  lp2 = 1
  n2 = i1
  lp1 = lptrp(1)
  n1 = listp(lp1)
!
!  Begin loop:  find the node N3 opposite N1->N2.
!
  do

    lp = lstptr(lend(n1),n2,list,lptr)

    if (list(lp) < 0) go to 4

    lp = lptr(lp)
    n3 = abs(list(lp))
!
!  Swap test:  Exit the loop if L = LMAX.
!
    if ( l >= lmax ) then
      exit
    end if

    dx11 = x(n1) - x(n3)
    dx12 = x(n2) - x(n3)
    dx22 = x(n2) - xp
    dx21 = x(n1) - xp

    dy11 = y(n1) - y(n3)
    dy12 = y(n2) - y(n3)
    dy22 = y(n2) - yp
    dy21 = y(n1) - yp

    cos1 = dx11*dx12 + dy11*dy12
    cos2 = dx22*dx21 + dy22*dy21

    if ( cos1 >= 0.0 .and. cos2 >= 0.0 ) go to 4

    if ( cos1 < 0.0 .and. cos2 < 0.0 ) go to 3

    sin1 = dx11*dy12 - dx12*dy11
    sin2 = dx22*dy21 - dx21*dy22
    if (sin1*cos2 + cos1*sin2 >= 0.0 ) go to 4
!
!  Swap:  Insert N3 following N2 in the adjacency list for P.
!  The two new arcs opposite P must be tested.
!
3   continue

    l = l+1
    lptrp(lp2) = l
    listp(l) = n3
    lptrp(l) = lp1
    lp1 = l
    n1 = n3
    cycle
!
!  No swap:  Advance to the next arc and test for termination
!  on N1 = I1 (LP1 = 1) or N1 followed by 0.
!
4   continue

    if ( lp1 == 1 ) then
      exit
    end if

    lp2 = lp1
    n2 = n1
    lp1 = lptrp(lp1)
    n1 = listp(lp1)

    if ( n1 == 0 ) then
      exit
    end if

  end do
!
!  Set NR and DSR to the index of the nearest node to P and
!  its squared distance from P, respectively.
!
  nr = i1
  dsr = (x(nr)-xp)**2 + (y(nr)-yp)**2

  do lp = 2, l

    n1 = listp(lp)

    if ( n1 == 0 ) then
      cycle
    end if

    ds1 = ( x(n1) - xp )**2 + ( y(n1) - yp )**2

    if ( ds1 < dsr ) then
      nr = n1
      dsr = ds1
    end if

  end do

  dsq = dsr
  nearnd = nr

  return

end function nearnd



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine Ready (TriangulationID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (TriangulationID > 0) then
            call LocateObjTriangulation (TriangulationID)
            ready_ = VerifyReadLock (mTRIANGULATION_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjTriangulation (TriangulationID)

        !Arguments-------------------------------------------------------------
        integer                                     :: TriangulationID

        !Local-----------------------------------------------------------------

        Me => FirstTriangulation
        do while (associated (Me))
            if (Me%InstanceID == TriangulationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleTriangulation - LocateObjTriangulation - ERR01'

    end subroutine LocateObjTriangulation

    !--------------------------------------------------------------------------


end Module ModuleTriangulation

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
