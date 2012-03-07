!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : Horizontal Map 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module which defines grid points into several categories (see below)
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
! WaterPoints2D         - All points where:
!                           - Bathymetry > -55.
!
! OpenPoints2D          - All points where 
!                           - ComputeFaces2D%U(i, j, k) + ComputeFaces2D%U(i, j+1, k) + 
!                             ComputeFaces2D%V(i, j, k) + ComputeFaces2D%V(i+1, j, k) > 0 
!
! LandPoints2D          - All points where:
!                           - Bathymetry < -90.
!
! BoundaryPoints2D      - All points which are:
!                           - Waterpoints at the limit of the domain
!                           - Waterpoints close to Bathymetry points of the value -80.
!
! ExteriorPoints2D      - All points where:
!                           - -90 < Bathymetry < -55
!                           - I = Size%ILB or I = Size%IUB
!                           - J = Size%JLB or J = Size%JUB
!
! ComputeFaces2D        - All faces which, simultaneous:
!                           - on both sides have Waterpoints
!                           - at least one of the Waterpoints have a waterlevel above H_Min
!                           - are NO face between boundary points
!
! ExteriorBoundaryFaces - All faces which:
!                           - Exterior faces of open boundary points
!                           - Faces between boundary points and points of the bathymetry with 
!                             the value of -80.
!
! BoundaryFaces         - All faces which:
!                           - have in one side a interior point and in another a boundary point
!
!
! WaterFaces            - All faces which:
!                           - can have water in both sides
                         
! IMin(:), IMax(:)      - Using these arrays it is possible to make:
!                           - DO I = IMin(J), IMax(J)
!                           - IMin is the first grid point with water
!                           - IMax is the last  grid point with water
!
!
module ModuleHorizontalMap

    use ModuleGlobalData
    use ModuleTime                 
    use ModuleGridData,          only : GetGridData, UngetGridData
    use ModuleHorizontalGrid,    only : GetHorizontalGridSize
    use ModuleStopWatch,         only : StartWatch, StopWatch         
    use ModuleFunctions,         only : Chunk_J, Chunk_K

    implicit none
    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructHorizontalMap
    private ::      AllocateInstance
    private ::      AllocateVariables
    private ::      ConstructBoundary
    private ::      ComputeIMaxIMin


    !Selector
    public  :: GetLandPoints2D
    public  :: GetWaterPoints2D
    public  :: GetOpenPoints2D
    public  :: GetExteriorPoints2D
    public  :: GetComputeFaces2D
    public  :: GetBoundaries
    public  :: GetExteriorBoundaryFaces
    public  :: GetWaterFaces2D    
    public  :: GetBoundaryFaces

    public  :: UngetHorizontalMap


    !Modifier
    public  :: UpdateComputeFaces2D
    private ::      UpdateOpenPoints2D
    public  :: UpdateWaterPoints2D

    !Destructor
    public  :: KillHorizontalMap


    !Management
    private ::      Ready
    private ::          LocateMe

    !Interfaces----------------------------------------------------------------
    private :: ConstructHorizontalMapWater
    private :: ConstructHorizontalMapLand
!    interface  ConstructHorizontalMap
!        module procedure ConstructHorizontalMapWater
!        module procedure ConstructHorizontalMapLand
!    end interface  ConstructHorizontalMap
   
    interface  UpdateComputeFaces2D
        module procedure UpdateComputeFaces2D_R4
        module procedure UpdateComputeFaces2D_R8
    end interface  UpdateComputeFaces2D


    
    
    private :: UngetHorizontalMap2D
    interface  UngetHorizontalMap
        module procedure UngetHorizontalMap2D
    end interface  UngetHorizontalMap


    !Parameter-----------------------------------------------------------------
    integer, parameter                          :: MapMohidWater = 1
    integer, parameter                          :: MapMohidLand  = 2        

    !Types---------------------------------------------------------------------
    type       T_2D_INT
        integer, pointer, dimension(:, :)       :: U, V
    end type T_2D_INT


    type      T_HorizontalMap
        integer                             :: InstanceID
        type(T_Size2D  )                    :: Size 
        type(T_Size2D  )                    :: WorkSize
        type(T_Time    )                    :: ActualTime
        integer                             :: MapType              = null_int

        !Instance of other modules
        integer                             :: ObjTopography        = 0
        integer                             :: ObjHorizontalGrid    = 0
                                            
        integer, pointer, dimension(:,:)    :: WaterPoints2D                !Mohid Water and Mohid Land
        integer, pointer, dimension(:,:)    :: OpenPoints2D                 !Mohid Water and Mohid Land
        integer, pointer, dimension(:,:)    :: LandPoints2D                 !Mohid Water and Mohid Land
        integer, pointer, dimension(:,:)    :: BoundaryPoints2D             !Mohid Water
        integer, pointer, dimension(:,:)    :: ExteriorPoints2D             !Mohid Water
        integer, pointer, dimension(:  )    :: IMin                         !Mohid Water
        integer, pointer, dimension(:  )    :: IMax                         !Mohid Water

        type(T_2D_INT  )                    :: ComputeFaces2D               !Mohid Water and Mohid Land
        type(T_2D_INT  )                    :: ExteriorBoundaryFaces        !Mohid Water
        type(T_2D_INT  )                    :: BoundaryFaces                !Mohid Water
        type(T_2D_INT  )                    :: WaterFaces                   !Mohid Water

        type (T_HorizontalMap), pointer     :: Next

    end type T_HorizontalMap

    !Global Module Variables
    type (T_HorizontalMap), pointer         :: Me, FirstHorizontalMap

    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine ConstructHorizontalMap(HorizontalMapID, GridDataID, HorizontalGridID,        &
                                      ActualTime, Points, Model, STAT) 
                                              
        !Arguments-------------------------------------------------------------
        integer, intent(INOUT)                     :: HorizontalMapID
        integer, intent(INOUT)                     :: GridDataID
        integer, intent(INOUT)                     :: HorizontalGridID
        type (T_Time), optional, intent(IN)        :: ActualTime
        integer, optional, dimension(:,:), pointer :: Points
        integer, optional, intent(OUT)             :: STAT
        integer, optional                          :: Model !1 - Water, 2 - Land

        !Local----------------------------------------------------------------- 
        integer                                    :: model_
        
        !----------------------------------------------------------------------
        if (present(Model)) then
            model_ = Model
        else
            model_ = 1 !Mohid Water
        endif
        
        if (model_ == 1) then
            call ConstructHorizontalMapWater (HorizontalMapID,  &
                                              GridDataID,       &
                                              HorizontalGridID, &
                                              ActualTime,       &
                                              Points,           &
                                              STAT)
        else
            call ConstructHorizontalMapLand  (HorizontalMapID,  &
                                              GridDataID,       &
                                              HorizontalGridID, &
                                              ActualTime,       &
                                              Points,           &
                                              STAT)        
        endif                                                                                                         
        !----------------------------------------------------------------------

    end subroutine ConstructHorizontalMap

    !--------------------------------------------------------------------------

    subroutine ConstructHorizontalMapWater(HorizontalMapID, GridDataID, HorizontalGridID,   &
                                           ActualTime, WaterPoints, STAT)  

        !Arguments-------------------------------------------------------------
        integer, intent(INOUT)                :: HorizontalMapID
        integer, intent(INOUT)                :: GridDataID
        integer, intent(INOUT)                :: HorizontalGridID
        type (T_Time), optional,   intent(IN) :: ActualTime
        integer, optional, dimension(:,:), pointer :: WaterPoints
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: ready_         
        real, pointer, dimension(: , :)     :: Bathymetry    
        integer                             :: i, j
        integer                             :: STAT_ 
        integer                             :: ILB, IUB, JLB, JUB
                                            
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mHorizontalMap_)) then
            nullify (FirstHorizontalMap)
            call RegisterModule (mHorizontalMap_) 
        endif

        call Ready(HorizontalMapID, ready_)    

cd2 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance           ()

            !Associated external Instances 
            Me%ObjTopography      = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalGrid  = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)

            !Actualize the time
            if (present(ActualTime)) Me%ActualTime = ActualTime

            !Sets Mapping to Water
            Me%MapType    = MapMohidWater

            !Gets the horizontal size from the Bathymetry
            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                             &
                                       Size        = Me%Size,                            &
                                       WorkSize    = Me%WorkSize,                        &
                                       STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHorizontalMapWater - ModuleHorizontalMap - ERR02'

            !Auxiliar variables for the do loops
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB
            
            !Allocates Variables
            call AllocateVariables()
            

            if (present(WaterPoints)) then

do5:            do j = JLB, JUB
do6:            do i = ILB, IUB 

                    Me%WaterPoints2D(i, j)    =     WaterPoints(i, j)
                    Me%LandPoints2D (i, j)    = 1 - WaterPoints(i, j)

                enddo do6
                enddo do5
            
            else 

                !Recieves the Bathymetry
                call GetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHorizontalMapWater - ModuleHorizontalMap - ERR03'

                
do1:            do j = JLB, JUB
do2:            do i = ILB, IUB 
cd3 :               if      ( Bathymetry(i, j) >-55.0) then

                        Me%WaterPoints2D(i, j)    = 1

                    else if ( Bathymetry(i, j) <-90.0) then cd3

                        Me%LandPoints2D(i, j)     = 1

                    else if ((Bathymetry(i, j) <-55.0) .AND. (Bathymetry(i, j) >-90.0)) then cd3

                        Me%ExteriorPoints2D(i, j) = 1

                    end if cd3
                end do do2
                end do do1


do3 :           do I = Me%Size%ILB, Me%Size%IUB
                    Me%ExteriorPoints2D(I,Me%Size%JLB) = 1
                    Me%ExteriorPoints2D(I,Me%Size%JUB) = 1
                end do do3

do4 :           do J = Me%Size%JLB, Me%Size%JUB
                    Me%ExteriorPoints2D(Me%Size%ILB,J) = 1
                    Me%ExteriorPoints2D(Me%Size%IUB,J) = 1
                end do do4


                !Computes Imax & Imin. These arrays are used to speed up the model avoiding the computation of land points. 
                ! Important for applications with small Water/Domain ratios.
                call ComputeIMaxIMin()
                
                
                !This subroutine search the matrix Bathymetry for boundary points
                call ConstructBoundary(Bathymetry)

                !This subroutine defines the water faces
                call ConstructWaterFaces


                !Ungets bathymetry
                call UngetGridData  (Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHorizontalMapWater - ModuleHorizontalMap - ERR03'
            
            endif


            !Returns ID
            HorizontalMapID = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            
            stop 'ModuleHorizontalMap - ConstructHorizontalMapWater - ERR99' 

        end if cd2


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructHorizontalMapWater

    !--------------------------------------------------------------------------

    subroutine ConstructHorizontalMapLand(HorizontalMapID, GridDataID, HorizontalGridID,   &
                                          ActualTime, BasinPoints, STAT)  

        !Arguments-------------------------------------------------------------
        integer                             :: HorizontalMapID
        integer                             :: GridDataID
        integer                             :: HorizontalGridID
        type (T_Time)                       :: ActualTime
        integer, dimension(:,:), pointer    :: BasinPoints
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        real, pointer, dimension(: , :)     :: Bathymetry    
        integer                             :: STAT_CALL
        integer                             :: ready_         
        integer                             :: STAT_ 
        integer                             :: i, j
                                            
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mHorizontalMap_)) then
            nullify (FirstHorizontalMap)
            call RegisterModule (mHorizontalMap_) 
        endif

        call Ready(HorizontalMapID, ready_)    

cd2 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance           ()

            !Associated external Instances 
            Me%ObjTopography      = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalGrid  = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)

            !Actualize the time
            Me%ActualTime = ActualTime

            !Sets Mapping to Water
            Me%MapType    = MapMohidLand

            !Gets the horizontal size from the Bathymetry
            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                             &
                                       Size        = Me%Size,                            &
                                       WorkSize    = Me%WorkSize,                        &
                                       STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHorizontalMapLand - ModuleHorizontalMap - ERR02'

            !Allocates Variables
            call AllocateVariables()

            !Sets WaterPoints equal to BasinsPoints
            Me%WaterPoints2D = BasinPoints

do3 :       do I = Me%Size%ILB, Me%Size%IUB
                Me%ExteriorPoints2D(I,Me%Size%JLB) = 1
                Me%ExteriorPoints2D(I,Me%Size%JUB) = 1
            end do do3

do4 :       do J = Me%Size%JLB, Me%Size%JUB
                Me%ExteriorPoints2D(Me%Size%ILB,J) = 1
                Me%ExteriorPoints2D(Me%Size%IUB,J) = 1
            end do do4

            !Recieves the Bathymetry
            call GetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHorizontalMapWater - ModuleHorizontalMap - ERR03'

            !This subroutine search the matrix Bathymetry for boundary points
            call ConstructBoundary(Bathymetry)

            !Ungets bathymetry
            call UngetGridData  (Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHorizontalMapWater - ModuleHorizontalMap - ERR03'


            !Returns ID
            HorizontalMapID = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            
            stop 'ModuleHorizontalMap - ConstructHorizontalMapLand - ERR99' 

        end if cd2


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructHorizontalMapLand

    !--------------------------------------------------------------------------

    subroutine AllocateInstance ()

        !Arguments-------------------------------------------------------------

    
        !Local-----------------------------------------------------------------
        type (T_HorizontalMap), pointer            :: NewMe
        type (T_HorizontalMap), pointer            :: PreviousMe


        !Allocates new instance
        allocate (NewMe)
        nullify  (NewMe%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstHorizontalMap)) then
            FirstHorizontalMap            => NewMe
            Me              => NewMe
        else
            PreviousMe      => FirstHorizontalMap
            Me              => FirstHorizontalMap%Next
            do while (associated(Me))
                PreviousMe  => Me
                Me          => Me%Next
            enddo
            Me              => NewMe
            PreviousMe%Next => NewMe
        endif

        Me%InstanceID = RegisterNewInstance (mHORIZONTALMAP_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ConstructBoundary (Bathymetry)   

        !Arguments-------------------------------------------------------------
        real, dimension(: ,:), pointer              :: Bathymetry

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        Me%BoundaryPoints2D       (:, :)   = 0
        Me%ExteriorBoundaryFaces%U(:, :)   = 0        
        Me%ExteriorBoundaryFaces%V(:, :)   = 0        

        !See if the matrix boundaries are also point boundaries
        i    = ILB
do1 :   do j = JLB , JUB
cd1 :       if (Bathymetry(i, j) > -55.) then
                Me%BoundaryPoints2D(i, j)          = 1
                Me%ExteriorBoundaryFaces%V(i, j)   = 1
            end if cd1
        end do do1

        i = IUB
do2 :   do j = JLB , JUB
cd2 :       if (Bathymetry(i, j) > -55.) then
                Me%BoundaryPoints2D(i, j)          = 1
                Me%ExteriorBoundaryFaces%V(i+1, j) = 1
            end if cd2
        end do do2

        j = JLB
do3 :   do i = ILB , IUB
cd3 :       if (Bathymetry(i, j) > -55.) then
                Me%BoundaryPoints2D(i, j)          = 1
                Me%ExteriorBoundaryFaces%U(i, j)   = 1
            end if cd3
        end do do3

        j = JUB
do4 :   do i = ILB , IUB
cd4 :       if (Bathymetry(i, j) > -55.) then
                Me%BoundaryPoints2D(i, j)          = 1
                Me%ExteriorBoundaryFaces%U(i, j+1) = 1
            end if cd4
        end do do4

        !Search for more boundaries points
do5 :   do j = JLB+1 , JUB-1
do6 :   do i = ILB+1 , IUB-1

            !lower cell of a waterpoint with -80?
            if (Bathymetry(i, j)    > -55.0     .and.  &
                Bathymetry(i-1, j)  < -55.0     .and.  &
                Bathymetry(i-1, j)  > -90.0)    then
                    Me%BoundaryPoints2D(i, j)          = 1
                    Me%ExteriorBoundaryFaces%V(i, j)   = 1
            endif
                    
            !left call of a waterpoint with -80?
            if (Bathymetry(i, j)    > -55.0     .and.  &
                Bathymetry(i, j-1)  < -55.0     .and.  & 
                Bathymetry(i, j-1)  > -90.0)    then
                    Me%BoundaryPoints2D(i, j)          = 1
                    Me%ExteriorBoundaryFaces%U(i, j)   = 1
            endif

            !upper cell of a waterpoint with -80?
            if (Bathymetry(i, j)    > -55.0     .and.  &
                Bathymetry(i+1, j)  < -55.0     .and.  &
                Bathymetry(i+1, j)  > -90.0)    then
                    Me%BoundaryPoints2D(i, j)          = 1
                    Me%ExteriorBoundaryFaces%V(i+1, j) = 1
            endif
            

            !right cell of a waterpoint with -80?
            if (Bathymetry(i, j)    > -55.0     .and.  &
                Bathymetry(i, j+1)  < -55.0     .and.  &
                Bathymetry(i, j+1)  > -90.0)    then
                    Me%BoundaryPoints2D(i, j)          = 1
                    Me%ExteriorBoundaryFaces%U(i, j+1) = 1
            endif
            


        end do do6
        end do do5


        Me%BoundaryFaces%U(:, :) = Not_Boundary
        Me%BoundaryFaces%V(:, :) = Not_Boundary

do7 :   do j = JLB+1 , JUB
do8 :   do i = ILB+1 , IUB
       
            if (Me%BoundaryPoints2D(i, j)     == Boundary    .and.  & 
                Me%WaterPoints2D   (i, j - 1) == WaterPoint  .and.  &
                Me%BoundaryPoints2D(i, j - 1) == Not_Boundary)  then 

                Me%BoundaryFaces%U(i, j) = Boundary

            endif


            if (Me%BoundaryPoints2D(i, j - 1) == Boundary    .and.    & 
                Me%WaterPoints2D   (i, j)     == WaterPoint  .and.    &
                Me%BoundaryPoints2D(i, j)     == Not_Boundary)  then 

                Me%BoundaryFaces%U(i, j) = Boundary

            endif

            if (Me%BoundaryPoints2D(i, j)     == Boundary     .and.  & 
                Me%WaterPoints2D   (i - 1, j) == WaterPoint   .and.  &
                Me%BoundaryPoints2D(i - 1, j) == Not_Boundary)  then 
                Me%BoundaryFaces%V(i, j) = Boundary
            endif


            if (Me%BoundaryPoints2D(i - 1, j) == Boundary     .and.    & 
                Me%WaterPoints2D   (i, j)     == WaterPoint   .and.    &
                Me%BoundaryPoints2D(i, j)     == Not_Boundary)  then 
                Me%BoundaryFaces%V(i, j) = Boundary
            endif

        end do do8
        end do do7

        !----------------------------------------------------------------------

    end subroutine ConstructBoundary

    !--------------------------------------------------------------------------
    subroutine ConstructWaterFaces

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !U direction
        do j = JLB , JUB + 1
        do i = ILB , IUB
                    
            if ((Me%WaterPoints2D(i,j-1) == WaterPoint .and.  Me%WaterPoints2D(i,j) == WaterPoint) .or. &
                Me%ExteriorBoundaryFaces%U(i, j)  == 1) then
             
                Me%WaterFaces%U(i, j) = 1
            else                
                Me%WaterFaces%U(i, j) = 0
            end if 
        end do
        enddo


        !V direction
        do j = JLB , JUB
        do i = ILB , IUB + 1
                    
            if ((Me%WaterPoints2D(i-1,j) == WaterPoint .and.  Me%WaterPoints2D(i,j) == WaterPoint) .or. &
                Me%ExteriorBoundaryFaces%V(i, j)  == 1) then
             
                Me%WaterFaces%V(i, j) = 1
            else
                Me%WaterFaces%V(i, j) = 0            
            end if 
        end do
        enddo


        !----------------------------------------------------------------------

    end subroutine ConstructWaterFaces

    !--------------------------------------------------------------------------

    subroutine AllocateVariables() 

        !External--------------------------------------------------------------

        !Local----------------------------------------------------------------

        integer :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        !Allocates Point variables common to MOHID Land and MOHID Water
        allocate(Me%WaterPoints2D               (ILB:IUB, JLB:JUB))
        allocate(Me%OpenPoints2D                (ILB:IUB, JLB:JUB))
        allocate(Me%LandPoints2D                (ILB:IUB, JLB:JUB))
        allocate(Me%BoundaryPoints2D            (ILB:IUB, JLB:JUB))
        allocate(Me%ExteriorPoints2D            (ILB:IUB, JLB:JUB))

        !Allocates face variables common to MOHID Land and MOHID Water
        allocate(Me%ComputeFaces2D%U            (ILB:IUB, JLB:JUB))
        allocate(Me%ComputeFaces2D%V            (ILB:IUB, JLB:JUB))

        !Allocates face variables of MOHID Water
        allocate(Me%ExteriorBoundaryFaces%U     (ILB:IUB, JLB:JUB))
        allocate(Me%ExteriorBoundaryFaces%V     (ILB:IUB, JLB:JUB))
        allocate(Me%BoundaryFaces%U             (ILB:IUB, JLB:JUB))
        allocate(Me%BoundaryFaces%V             (ILB:IUB, JLB:JUB))
        allocate(Me%WaterFaces%U                (ILB:IUB, JLB:JUB))
        allocate(Me%WaterFaces%V                (ILB:IUB, JLB:JUB))


        !Allocates Imin & Imax
        allocate(Me%Imin                        (         JLB:JUB))
        allocate(Me%Imax                        (         JLB:JUB))

        Me%WaterPoints2D(:,:)           = 0     !By Default there are no WaterPoints
        Me%OpenPoints2D (:,:)           = 0     
        Me%LandPoints2D(:,:)            = 0     !By Default there are no LandPoints
        
        Me%BoundaryPoints2D(:,:)        = 0         !By Default there are no BoundaryPoints
        Me%ExteriorPoints2D(:,:)        = 0 

        Me%ComputeFaces2D%U             = 0
        Me%ComputeFaces2D%V             = 0
        
        Me%ExteriorBoundaryFaces%U(:,:) = 0   
        Me%ExteriorBoundaryFaces%V(:,:) = 0

        Me%BoundaryFaces%U(:,:)         = 0    
        Me%BoundaryFaces%V(:,:)         = 0

        Me%WaterFaces%U(:,:)            = 0    
        Me%WaterFaces%V(:,:)            = 0


        Me%Imin(:)                      = 0    
        Me%Imax(:)                      = 0

        !----------------------------------------------------------------------

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine ComputeIMaxIMin()   
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, i, j
                                
        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !If a row has no waterpoints it will compute only the first point
        Me%IMin(:) = ILB
        Me%IMax(:) = ILB

        !Fills IMin
do1 :   do J = JLB, JUB
do2 :   do I = ILB, IUB
if1 :       if (Me%WaterPoints2D(I, J) == 1) then
                Me%IMin(J) = I
                exit do2
            endif   if1
        enddo       do2
        enddo       do1

        !Fills IMax
do3 :   do J = JLB, JUB
do4 :   do I = IUB, ILB, -1
if3 :       if (Me%WaterPoints2D(I, J) == 1) then
                Me%IMax(J) = I+1
                exit do4
            endif   if3
        enddo       do4
        enddo       do3

        !----------------------------------------------------------------------

    end subroutine ComputeIMaxIMin

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine UpdateComputeFaces2D_R4(HorizontalMapID, WaterElevation, H_MIN, ActualTime, STAT)   
        
        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        real(4), dimension(:, :), pointer           :: WaterElevation      
        real                                        :: H_MIN      
        type (T_Time)                               :: ActualTime
        integer, optional, intent(OUT)              :: STAT      

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        real                                        :: Min_Level_East, Min_Level_West
        real                                        :: Min_Level_South, Min_Level_North
        real                                        :: Level_East, Level_West
        real                                        :: Level_North, Level_South
        integer                                     :: ComputeFaceEastWest, ComputeFaceWestEast
        integer                                     :: ComputeFaceNorthSouth, ComputeFaceSouthNorth
        integer                                     :: ILB, IUB, JLB, JUB, I, J
        integer                                     :: STAT_
        real, dimension(:,:), pointer               :: GridData
        integer                                     :: STAT_CALL
        real                                        :: Factor
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            Me%ComputeFaces2D%U = 0
            Me%ComputeFaces2D%V = 0

            !Actualize the time
            Me%ActualTime = ActualTime

            !Gets Bathymetry/Topography
            call GetGridData              (Me%ObjTopography, GridData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateComputeFaces2D_R4 - ModuleHorizontalMap - ERR01'

            if      (Me%MapType == MapMohidWater) then
                Factor = -1.0
            elseif  (Me%MapType == MapMohidLand ) then
                Factor =  1.0
            else
                stop 'UpdateComputeFaces2D_R4 - ModuleHorizontalMap - ERR02'
            endif

            if (MonitorPerformance) call StartWatch ("ModuleHorizontalMap", "UpdateComputeFaces2D_R4")

            CHUNK = CHUNK_J(JLB,JUB)
            !$OMP PARALLEL PRIVATE(I,J,Min_Level_East,Min_Level_West,Level_East) &
            !$OMP PRIVATE(Level_West,ComputeFaceEastWest,ComputeFaceWestEast) &
            !$OMP PRIVATE(Min_Level_North,Min_Level_South,Level_North,Level_South) &
            !$OMP PRIVATE(ComputeFaceNorthSouth,ComputeFaceSouthNorth)

            !Checks for ComputeFaces2D in X direction
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do1:        do J = JLB+1, JUB
do2:        do I = ILB, IUB

                !By default can discharge from East to West, neither from West to East
                ComputeFaceEastWest = 0
                ComputeFaceWestEast = 0

                Min_Level_East = Factor * GridData(I, J)   + H_MIN
                Min_Level_West = Factor * GridData(I, J-1) + H_MIN
                Level_East     = WaterElevation(I, J)
                Level_West     = WaterElevation(I, J-1)

                !If on the east side the waterlevel is higher then the eastern minimum level AND
                !   on the east side the waterlevel is higher then the western minimum level THEN
                !its possible to transfer water from EAST to WEST
                if (Level_East > Min_Level_East .and. &
                    Level_East > Min_Level_West)      &
                    ComputeFaceEastWest = 1

                !If on the west side the waterlevel is higher then the western minimum level AND
                !   on the west side the waterlevel is higher then the eastern minimum level THEN
                !its possible to transfer water from WEST to EAST
                if (Level_West > Min_Level_West .and. &
                    Level_West > Min_Level_East)      &
                    ComputeFaceWestEast = 1

                !When its possible to transfer water from WEST to EAST or from EAST to WEST then
                !its a compute face
                if (ComputeFaceEastWest + ComputeFaceWestEast > 0) then
                    Me%ComputeFaces2D%U(I, J)   = 1
                endif


                !ComputeFaces2D just exists, if both sides are waterpoints
                Me%ComputeFaces2D%U(I, J) = Me%ComputeFaces2D%U(I, J) * Me%WaterPoints2D(I, J)

                Me%ComputeFaces2D%U(I, J) = Me%ComputeFaces2D%U(I, J) * Me%WaterPoints2D(I, J-1)

                !Faces between BoundaryPoints aren't considered
                if ((Me%BoundaryPoints2D(i, j) == 1) .and. (Me%BoundaryPoints2D(i, j-1) == 1)) then
                    Me%ComputeFaces2D%U(i, j) = 0
                endif
                                
            enddo do2
            enddo do1
            !$OMP END DO

            !Checks for ComputeFaces2D in Y direction
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3:        do J = JLB, JUB
do4:        do I = ILB+1, IUB


                !By default can discharge from North to South, neither from South to North
                ComputeFaceNorthSouth = 0
                ComputeFaceSouthNorth = 0

                Min_Level_North = Factor * GridData(I, J)   + H_MIN
                Min_Level_South = Factor * GridData(I-1, J) + H_MIN
                Level_North     = WaterElevation(I, J)
                Level_South     = WaterElevation(I-1, J)

                !If on the north side the waterlevel is higher then the northern minimum level AND
                !   on the north side the waterlevel is higher then the southern minimum level THEN
                !its possible to transfer water from NORTH to SOUTH
                if (Level_North > Min_Level_North .and. &
                    Level_North > Min_Level_South)      &
                    ComputeFaceNorthSouth = 1

                !If on the north side the waterlevel is higher then the northern minimum level AND
                !   on the north side the waterlevel is higher then the southern minimum level THEN
                !its possible to transfer water from NORTH to SOUTH
                if (Level_South > Min_Level_South .and. &
                    Level_South > Min_Level_North)      &
                    ComputeFaceSouthNorth = 1

                !When its possible to transfer water from NORTH to SOUTH or from SOUTH to NORTH then
                !its a compute face
                if (ComputeFaceNorthSouth + ComputeFaceSouthNorth > 0) then
                    Me%ComputeFaces2D%V(I, J)   = 1
                endif


                !ComputeFaces2D just exists, if both sides are waterpoints
                Me%ComputeFaces2D%V(I, J) = Me%ComputeFaces2D%V(I, J) * Me%WaterPoints2D(I, J)

                Me%ComputeFaces2D%V(I, J) = Me%ComputeFaces2D%V(I, J) * Me%WaterPoints2D(I-1, J)

                !Faces between BoundaryPoints aren't considered
                if ((Me%BoundaryPoints2D(i, j) == 1) .and. (Me%BoundaryPoints2D(i-1, j) == 1)) then
                    Me%ComputeFaces2D%V(I, J) = 0
                endif

            enddo do4
            enddo do3
            !$OMP END DO NOWAIT

            !$OMP END PARALLEL

            if (MonitorPerformance) call StopWatch ("ModuleHorizontalMap", "UpdateComputeFaces2D_R4")

            call UpdateOpenPoints2D()

            !Ungets Bathymetry/Topography
            call UnGetGridData     (Me%ObjTopography, GridData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateComputeFaces2D_R4 - ModuleHorizontalMap - ERR02'


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UpdateComputeFaces2D_R4

    !----------------------------------------------------------------------

    subroutine UpdateComputeFaces2D_R8(HorizontalMapID, WaterElevation, H_MIN, ActualTime, STAT)   
        
        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        real(8), dimension(:, :), pointer           :: WaterElevation      
        real                                        :: H_MIN      
        type (T_Time)                               :: ActualTime
        integer, optional, intent(OUT)              :: STAT      

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        real                                        :: Min_Level_East, Min_Level_West
        real                                        :: Min_Level_South, Min_Level_North
        real                                        :: Level_East, Level_West
        real                                        :: Level_North, Level_South
        integer                                     :: ComputeFaceEastWest, ComputeFaceWestEast
        integer                                     :: ComputeFaceNorthSouth, ComputeFaceSouthNorth
        integer                                     :: ILB, IUB, JLB, JUB, I, J
        integer                                     :: STAT_
        real, dimension(:,:), pointer               :: GridData
        integer                                     :: STAT_CALL
        real                                        :: Factor
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            Me%ComputeFaces2D%U = 0
            Me%ComputeFaces2D%V = 0

            !Actualize the time
            Me%ActualTime = ActualTime

            !Gets Bathymetry/Topography
            call GetGridData              (Me%ObjTopography, GridData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateComputeFaces2D_R8 - ModuleHorizontalMap - ERR00'

            if      (Me%MapType == MapMohidWater) then
                Factor = -1.0
            elseif  (Me%MapType == MapMohidLand ) then
                Factor =  1.0
            else
                stop 'UpdateComputeFaces2D_R8 - ModuleHorizontalMap - ERR01'
            endif

            if (MonitorPerformance) call StartWatch ("ModuleHorizontalMap", "UpdateComputeFaces2D_R8")

            CHUNK = CHUNK_J(JLB,JUB)
            !$OMP PARALLEL PRIVATE(I,J,Min_Level_East,Min_Level_West,Level_East) &
            !$OMP PRIVATE(Level_West,ComputeFaceEastWest,ComputeFaceWestEast) &
            !$OMP PRIVATE(Min_Level_North,Min_Level_South,Level_North,Level_South) &
            !$OMP PRIVATE(ComputeFaceNorthSouth,ComputeFaceSouthNorth)

            !Checks for ComputeFaces2D in X direction
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do1:        do J = JLB+1, JUB
do2:        do I = ILB, IUB

                !By default can discharge from East to West, neither from West to East
                ComputeFaceEastWest = 0
                ComputeFaceWestEast = 0

                Min_Level_East = Factor * GridData(I, J)   + H_MIN
                Min_Level_West = Factor * GridData(I, J-1) + H_MIN
                Level_East     = WaterElevation(I, J)
                Level_West     = WaterElevation(I, J-1)

                !If on the east side the waterlevel is higher then the eastern minimum level AND
                !   on the east side the waterlevel is higher then the western minimum level THEN
                !its possible to transfer water from EAST to WEST
                if (Level_East > Min_Level_East .and. &
                    Level_East > Min_Level_West)      &
                    ComputeFaceEastWest = 1

                !If on the west side the waterlevel is higher then the western minimum level AND
                !   on the west side the waterlevel is higher then the eastern minimum level THEN
                !its possible to transfer water from WEST to EAST
                if (Level_West > Min_Level_West .and. &
                    Level_West > Min_Level_East)      &
                    ComputeFaceWestEast = 1

                !When its possible to transfer water from WEST to EAST or from EAST to WEST then
                !its a compute face
                if (ComputeFaceEastWest + ComputeFaceWestEast > 0) then
                    Me%ComputeFaces2D%U(I, J)   = 1
                endif


                !ComputeFaces2D just exists, if both sides are waterpoints
                Me%ComputeFaces2D%U(I, J) = Me%ComputeFaces2D%U(I, J) * Me%WaterPoints2D(I, J)

                Me%ComputeFaces2D%U(I, J) = Me%ComputeFaces2D%U(I, J) * Me%WaterPoints2D(I, J-1)

                !Faces between BoundaryPoints aren't considered
                if ((Me%BoundaryPoints2D(i, j) == 1) .and. (Me%BoundaryPoints2D(i, j-1) == 1)) then
                    Me%ComputeFaces2D%U(i, j) = 0
                endif
                                
            enddo do2
            enddo do1
            !$OMP END DO

            !Checks for ComputeFaces2D in Y direction
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3:        do J = JLB, JUB
do4:        do I = ILB+1, IUB


                !By default can discharge from North to South, neither from South to North
                ComputeFaceNorthSouth = 0
                ComputeFaceSouthNorth = 0

                Min_Level_North = Factor * GridData(I, J)   + H_MIN
                Min_Level_South = Factor * GridData(I-1, J) + H_MIN
                Level_North     = WaterElevation(I, J)
                Level_South     = WaterElevation(I-1, J)

                !If on the north side the waterlevel is higher then the northern minimum level AND
                !   on the north side the waterlevel is higher then the southern minimum level THEN
                !its possible to transfer water from NORTH to SOUTH
                if (Level_North > Min_Level_North .and. &
                    Level_North > Min_Level_South)      &
                    ComputeFaceNorthSouth = 1

                !If on the north side the waterlevel is higher then the northern minimum level AND
                !   on the north side the waterlevel is higher then the southern minimum level THEN
                !its possible to transfer water from NORTH to SOUTH
                if (Level_South > Min_Level_South .and. &
                    Level_South > Min_Level_North)      &
                    ComputeFaceSouthNorth = 1

                !When its possible to transfer water from NORTH to SOUTH or from SOUTH to NORTH then
                !its a compute face
                if (ComputeFaceNorthSouth + ComputeFaceSouthNorth > 0) then
                    Me%ComputeFaces2D%V(I, J)   = 1
                endif


                !ComputeFaces2D just exists, if both sides are waterpoints
                Me%ComputeFaces2D%V(I, J) = Me%ComputeFaces2D%V(I, J) * Me%WaterPoints2D(I, J)

                Me%ComputeFaces2D%V(I, J) = Me%ComputeFaces2D%V(I, J) * Me%WaterPoints2D(I-1, J)

                !Faces between BoundaryPoints aren't considered
                if ((Me%BoundaryPoints2D(i, j) == 1) .and. (Me%BoundaryPoints2D(i-1, j) == 1)) then
                    Me%ComputeFaces2D%V(I, J) = 0
                endif

            enddo do4
            enddo do3
            !$OMP END DO NOWAIT

            !$OMP END PARALLEL

            if (MonitorPerformance) call StopWatch ("ModuleHorizontalMap", "UpdateComputeFaces2D_R8")

            call UpdateOpenPoints2D()

            !Ungets Bathymetry/Topography
            call UnGetGridData     (Me%ObjTopography, GridData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateComputeFaces2D_R8 - ModuleHorizontalMap - ERR02'


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UpdateComputeFaces2D_R8


    !--------------------------------------------------------------------------

    subroutine UpdateOpenPoints2D()   
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, i, j
        integer                                     :: CHUNK
                                
        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        if (MonitorPerformance) call StartWatch ("ModuleHorizontalMap", "UpdateOpenPoints2D")

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(I,J)

        !Searches for open points
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ComputeFaces2D%U(i, j) + Me%ComputeFaces2D%U(i  , j+1) +  & 
                Me%ComputeFaces2D%V(i, j) + Me%ComputeFaces2D%V(i+1, j  ) > 0 ) then
                Me%OpenPoints2D(i, j) = 1
            else
                Me%OpenPoints2D(i, j) = 0
            endif
        enddo
        enddo

        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleHorizontalMap", "UpdateOpenPoints2D")

        !----------------------------------------------------------------------

    end subroutine UpdateOpenPoints2D


    !--------------------------------------------------------------------------

    subroutine UpdateWaterPoints2D(HorizontalMapID, WaterPoints, STAT)  

        !Arguments-------------------------------------------------------------
        integer, intent(INOUT)                :: HorizontalMapID
        integer, optional, dimension(:,:), pointer :: WaterPoints
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: ready_         
        real, pointer, dimension(: , :)     :: Bathymetry    
        integer                             :: i, j
        integer                             :: STAT_ 
        integer                             :: ILB, IUB, JLB, JUB
                                            
        !----------------------------------------------------------------------


        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            if (present(WaterPoints)) then

do5:            do j = JLB, JUB
do6:            do i = ILB, IUB 

                    Me%WaterPoints2D(i, j)    =     WaterPoints(i, j)
                    Me%LandPoints2D (i, j)    = 1 - WaterPoints(i, j)

                enddo do6
                enddo do5
            
            else 

                !Recieves the Bathymetry
                call GetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'UpdateWaterPoints2D - ModuleHorizontalMap - ERR10'

                
do1:            do j = JLB, JUB
do2:            do i = ILB, IUB 

                    Me%WaterPoints2D   (i, j) = 0
                    Me%LandPoints2D    (i, j) = 0
                    Me%ExteriorPoints2D(i, j) = 0
                    
cd3 :               if      ( Bathymetry(i, j) >-55.0) then

                        Me%WaterPoints2D(i, j)    = 1

                    else if ( Bathymetry(i, j) <-90.0) then cd3

                        Me%LandPoints2D(i, j)     = 1

                    else if ((Bathymetry(i, j) <-55.0) .AND. (Bathymetry(i, j) >-90.0)) then cd3

                        Me%ExteriorPoints2D(i, j) = 1

                    end if cd3
                end do do2
                end do do1


do3 :           do I = Me%Size%ILB, Me%Size%IUB
                    Me%ExteriorPoints2D(I,Me%Size%JLB) = 1
                    Me%ExteriorPoints2D(I,Me%Size%JUB) = 1
                end do do3

do4 :           do J = Me%Size%JLB, Me%Size%JUB
                    Me%ExteriorPoints2D(Me%Size%ILB,J) = 1
                    Me%ExteriorPoints2D(Me%Size%IUB,J) = 1
                end do do4


                !Computes Imax & Imin. These arrays are used to speed up the model avoiding the computation of land points. 
                ! Important for applications with small Water/Domain ratios.
                call ComputeIMaxIMin()
                
                
                !This subroutine search the matrix Bathymetry for boundary points
                call ConstructBoundary(Bathymetry)

                !This subroutine defines the water faces
                call ConstructWaterFaces


                !Ungets bathymetry
                call UngetGridData  (Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'UpdateWaterPoints2D - ModuleHorizontalMap - ERR20'

            endif

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine UpdateWaterPoints2D

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetLandPoints2D(HorizontalMapID, LandPoints2D, STAT)  

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer,  pointer, dimension(: , :)         :: LandPoints2D
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)
            LandPoints2D => Me%LandPoints2D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetLandPoints2D

    !--------------------------------------------------------------------------

    subroutine GetWaterPoints2D(HorizontalMapID, WaterPoints2D, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer,  pointer, dimension(: , :)         :: WaterPoints2D
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)

            WaterPoints2D => Me%WaterPoints2D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWaterPoints2D

    !----------------------------------------------------------------------------

    subroutine GetOpenPoints2D(HorizontalMapID, OpenPoints2D, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer,  pointer, dimension(: , :)         :: OpenPoints2D
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)

            OpenPoints2D => Me%OpenPoints2D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetOpenPoints2D

    !----------------------------------------------------------------------------

    subroutine GetExteriorPoints2D(HorizontalMapID, ExteriorPoints2D, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer,  pointer, dimension(: , :)         :: ExteriorPoints2D
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
                
cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)

            ExteriorPoints2D => Me%ExteriorPoints2D

            STAT_ = SUCCESS_
        else  cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetExteriorPoints2D

    !----------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine GetComputeFaces2D(HorizontalMapID, ComputeFaces2DU, &
                                 ComputeFaces2DV, ActualTime, STAT)   

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer, dimension(:,:), optional, pointer  :: ComputeFaces2DU
        integer, dimension(:,:), optional, pointer  :: ComputeFaces2DV
        type (T_Time)          , optional           :: ActualTime
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)
            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)


            !Verifies if the distances are up to date
            if (present(ActualTime)) then
                if (ActualTime .ne. Me%ActualTime) &
                    stop 'GetComputeFaces2D - ModuleHorizontalMap - ERR01'
            endif

            ComputeFaces2DU => Me%ComputeFaces2D%U
            ComputeFaces2DV => Me%ComputeFaces2D%V


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetComputeFaces2D

    !----------------------------------------------------------------------

    subroutine GetBoundaries(HorizontalMapID, BoundaryPoints2D, STAT)  

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer, dimension(:, :) , pointer          :: BoundaryPoints2D
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)

            BoundaryPoints2d => Me%BoundaryPoints2D

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetBoundaries

    !--------------------------------------------------------------------------

    subroutine GetExteriorBoundaryFaces(HorizontalMapID, BoundaryPointsFaceU, &
                                        BoundaryPointsFaceV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer, dimension(:, :), pointer, optional :: BoundaryPointsFaceU
        integer, dimension(:, :), pointer, optional :: BoundaryPointsFaceV
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

cd2:        if (present(BoundaryPointsFaceU)) then

                call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)
                BoundaryPointsFaceU => Me%ExteriorBoundaryFaces%U

            endif cd2

cd3:        if (present(BoundaryPointsFaceV)) then

                call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)
                BoundaryPointsFaceV => Me%ExteriorBoundaryFaces%V

            endif cd3


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetExteriorBoundaryFaces

    !--------------------------------------------------------------------------

    subroutine GetBoundaryFaces(HorizontalMapID, BoundaryFacesU, &
                                BoundaryFacesV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer, dimension(:, :)  , pointer         :: BoundaryFacesU, BoundaryFacesV
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)
            BoundaryFacesU => Me%BoundaryFaces%U

            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)
            BoundaryFacesV => Me%BoundaryFaces%V


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetBoundaryFaces


    !--------------------------------------------------------------------------

    subroutine GetWaterFaces2D(HorizontalMapID, WaterFaces2DU, &
                                WaterFaces2DV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer, dimension(:, :)  , pointer         :: WaterFaces2DU, WaterFaces2DV
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)
            WaterFaces2DU => Me%WaterFaces%U

            call Read_Lock(mHORIZONTALMAP_, Me%InstanceID)
            WaterFaces2DV => Me%WaterFaces%V


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWaterFaces2D


    !--------------------------------------------------------------------------

    subroutine UngetHorizontalMap2D(HorizontalMapID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer, pointer, dimension(:,:)            :: Array
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mHORIZONTALMAP_, Me%InstanceID, "UngetHorizontalMap2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHorizontalMap2D


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillHorizontalMap(HorizontalMapID, STAT)  

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer, optional, intent(OUT)              :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_              
        integer                                     :: ready_   
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalMapID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mHORIZONTALMAP_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deassociated external Instances 
                nUsers = DeassociateInstance (mGRIDDATA_,       Me%ObjTopography    )
                if (nUsers == 0) stop 'KillHorizontalMap - ModuleHorizontalMap - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillHorizontalMap - ModuleHorizontalMap - ERR02'

                deallocate(Me%WaterPoints2D                 )
                deallocate(Me%LandPoints2D                  )
                deallocate(Me%OpenPoints2D                  )
                deallocate(Me%ComputeFaces2D%U              )
                deallocate(Me%ComputeFaces2D%V              )

                deallocate(Me%BoundaryPoints2D              )
                deallocate(Me%ExteriorPoints2D              )
                deallocate(Me%ExteriorBoundaryFaces%U       )
                deallocate(Me%ExteriorBoundaryFaces%V       )
                deallocate(Me%BoundaryFaces%U               )
                deallocate(Me%BoundaryFaces%V               )
                deallocate(Me%WaterFaces%U                  )
                deallocate(Me%WaterFaces%V                  )
                
                deallocate(Me%Imin                          )
                deallocate(Me%Imax                          )

                call DeallocateInstance ()

                HorizontalMapID = 0
                STAT_           = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillHorizontalMap

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HorizontalMap), pointer          :: AuxMe
        type (T_HorizontalMap), pointer          :: PreviousMe

        !Updates pointers
        if (Me%InstanceID == FirstHorizontalMap%InstanceID) then
            FirstHorizontalMap => FirstHorizontalMap%Next
        else
            PreviousMe => FirstHorizontalMap
            AuxMe      => FirstHorizontalMap%Next
            do while (AuxMe%InstanceID /= Me%InstanceID)
                PreviousMe => AuxMe
                AuxMe      => AuxMe%Next
            enddo

            !Now update linked list
            PreviousMe%Next => AuxMe%Next

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

    subroutine Ready (Me_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: Me_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (Me_ID > 0) then
            call LocateMe(Me_ID)
            ready_ = VerifyReadLock (mHORIZONTALMAP_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateMe (MeID)

        !Arguments-------------------------------------------------------------
        integer                                     :: MeID

        !Local-----------------------------------------------------------------

        Me => FirstHorizontalMap
        do while (associated (Me))
            if (Me%InstanceID == MeID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                          &
            stop 'ModuleHorizontalMap - LocateMe - ERR01'

    end subroutine LocateMe

    !--------------------------------------------------------------------------

#ifdef _OPENMI_

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::IsWaterPoint
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_ISWATERPOINT"::IsWaterPoint
    !DEC$ ENDIF
    logical function IsWaterPoint(HorizontalMapID, i, j)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalMapID
        integer                                     :: i, j
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(HorizontalMapID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            if (Me%WaterPoints2D(i, j) == 1) then
                IsWaterPoint = .true.
            else
                IsWaterPoint = .false.
            endif
        else 
            IsWaterPoint = .false.
        end if
           
        return

    end function IsWaterPoint


#endif


end module ModuleHorizontalMap

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

