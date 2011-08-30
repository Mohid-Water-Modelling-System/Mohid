!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Digital Terrain Creator
! MODULE        : Drawing
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes / Frank Braunschweig - v4.0
! DESCRIPTION   : Module to handle Geometric objects
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

Module ModuleDrawing

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleFunctions

    implicit none
    
    private
    
    !Subroutines--------------------------------------------------------------
    public  ::    New
    private ::    Add
    public  ::    SetLimits
    public  ::    IsVisible
    public  ::    IsPointInsidePolygon
    public  ::    Intersect2D_SegPoly
    public  ::    IsPointInsideCircle
    public  ::    WriteItem
    public  ::    PointDistanceToPolygon
    public  ::    InvertVerticesOrder
    public  ::    CheckPolygonClockwise
    public  ::    ConvexPolygon

    public  ::    VertPolygonInsidePolygon
    
    public  ::    SliceCellIn4
    
    public  ::    ArrayPolygonWindow
    public  ::    CellInterSectCell
    public  ::    SegIntersectSeg
   
    private ::    NewPolygon
    private ::    NewXYZPoint
    private ::    NewXYZPoint_V2
    private ::    NewXYZPoint_V3
    private ::    NewXYZPoint_V4
    interface     New
        module procedure NewPolygon
        module procedure NewXYZPoint
        module procedure NewXYZPoint_V2 
        module procedure NewXYZPoint_V3 
        module procedure NewXYZPoint_V4 
        module procedure NewLine
    end interface New
    
    private ::    AddPolygon
    private ::    AddXYZPoint
    private ::    AddLine
    interface     Add
        module procedure AddPolygon
        module procedure AddXYZPoint
        module procedure AddLine
    end interface Add

    private ::    SetLimitsPolygon
    private ::    SetLimitsXYZ
    interface     SetLimits
        module procedure SetLimitsPolygon
        module procedure SetLimitsXYZ
    end interface SetLimits

    private ::    WriteItemPolygon
    private ::    WriteItemXYZ
    private ::    WriteItemXYZ_v2
    interface     WriteItem
        module procedure WriteItemPolygon
        module procedure WriteItemXYZ
        module procedure WriteItemXYZ_v2
    end interface WriteItem

    private :: EqualPoint
    private :: LessPoint
    private :: DotProduct2D
    private :: PerpProduct2D
        

    !Parameter-----------------------------------------------------------------
    integer(4), parameter  :: TypeX_Y_Z     =  1
    integer(4), parameter  :: TypeX_Y_Z_P   =  2    
    integer(4), parameter  :: TypeT_X_Y_Z_P =  3    
    
    
    !Types--------------------------------------------------------------
    public   T_Limits
    type     T_Limits
        real                                    :: Left         = null_real
        real                                    :: Right        = null_real
        real                                    :: Top          = null_real
        real                                    :: Bottom       = null_real
        real                                    :: MaximumValue = 999999.
        real                                    :: MinimumValue =-999999.
    end type T_Limits
    
    public   T_Point                                                
    type     T_Point                                                
        integer                                 :: I            = null_int
        integer                                 :: J            = null_int
    end type T_Point                                                                                          

    public   T_PointF
    type     T_PointF
        real                                    :: X            = null_real
        real                                    :: Y            = null_real
    end type T_PointF
   
    public   T_XYZPoints
    type     T_XYZPoints
        integer                                 :: ID           = null_int
        real,           dimension(:), pointer   :: X
        real,           dimension(:), pointer   :: Y
        real,           dimension(:), pointer   :: Z
        real,           dimension(:), pointer   :: T
        real,           dimension(:), pointer   :: P                
        logical,        dimension(:), pointer   :: Inside
        integer                                 :: FileFormat
        integer                                 :: Count
        type(T_Limits)                          :: Limits
        type (T_XYZPoints),           pointer   :: Next
    end type T_XYZPoints
    
    public   T_Segment
    type     T_Segment
        type(T_PointF), pointer                 :: StartAt
        type(T_PointF), pointer                 :: EndAt
    end type T_Segment

    public   T_Polygon
    type     T_Polygon
        integer                                 :: ID           = null_int
        type(T_Point ),  dimension(:), pointer  :: Vertices
        type(T_PointF),  dimension(:), pointer  :: VerticesF
        type(T_Limits)                          :: Limits
        integer                                 :: Count
        type(T_Polygon),              pointer   :: Next 
    end type T_Polygon


    public   T_Lines
    type     T_Lines
        integer                                 :: ID           = null_int
        real, dimension(:), pointer             :: X
        real, dimension(:), pointer             :: Y
        integer                                 :: nNodes    
        type(T_Lines) , pointer                 :: Next
    end type T_Lines

    contains
    
    !--------------------------------------------------------------------------
    
    subroutine NewPolygon(Polygons, PolygonsFileName)
        
        !Arguments-------------------------------------------------------------
        type(T_Polygon),                 pointer    :: Polygons
        character(len=*),  intent(IN)               :: PolygonsFileName

        !Local-----------------------------------------------------------------
        type(T_Polygon),                 pointer    :: Polygon
        integer                                     :: PolygonsFile = 0
        integer                                     :: ClientNumber, StartLine, EndLine
        integer                                     :: FromBlock, CurrentLineNumber, CurrentVertix
        character(len = line_length)                :: FullBufferLine
        logical                                     :: BlockFound = OFF
        real, dimension(:),              pointer    :: PointCoordinates
        integer                                     :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        call ConstructEnterData(PolygonsFile, PolygonsFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'NewPolygon - ModuleDrawing - ERR01'
    
do1 :       do
                call ExtractBlockFromBuffer(PolygonsFile,                                &
                                            ClientNumber,                                &
                                            '<beginpolygon>',                           &
                                            '<endpolygon>',                             &
                                            BlockFound,                                  &
                                            STAT = STAT_CALL)

if1 :           if      (STAT_CALL .EQ. SUCCESS_ ) then    

if2 :               if (BlockFound) then

                        call Add(Polygons, Polygon)

                        call GetExtractType(FromBlock = FromBlock)

                        call GetBlockSize(PolygonsFile,                                  &
                                          ClientNumber,                                  &
                                          StartLine,                                     &
                                          EndLine,                                       &
                                          FromBlock,                                     &
                                          STAT = STAT_CALL)

                        !ATTENTION!!! there is a bug in module EnterData regarding
                        !             the STAT variable in sub GetBlockSize, see also
                        !             modules Lagrangian and TimeSerie

                        Polygon%Count = EndLine - StartLine 
                    
                        allocate(Polygon%VerticesF(1:Polygon%Count))
                        allocate(PointCoordinates (1:2))
                    
                        CurrentVertix = 1

                        do CurrentLineNumber = StartLine + 1 , EndLine - 1 
                        
                            call GetFullBufferLine(PolygonsFile,                         &
                                                   CurrentLineNumber,                    &
                                                   FullBufferLine,                       &
                                                   STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_)stop 'NewPolygon - ModuleDrawing - ERR02'
                        
                            call GetExtractType(FromBlock = FromBlock)


                            call GetData(PointCoordinates,                               &
                                         PolygonsFile,                                   &
                                         flag = iflag,                                   &
                                         SearchType = FromBlock,                         &
                                         Buffer_Line = CurrentLineNumber,                &
                                         STAT        = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_ .and.                             &
                               STAT_CALL .ne. SIZE_ERR_) stop 'NewPolygon - ModuleDrawing - ERR03'
                            
                            Polygon%VerticesF(CurrentVertix)%X  = PointCoordinates(1)
                            Polygon%VerticesF(CurrentVertix)%Y  = PointCoordinates(2)
                       
                            CurrentVertix = CurrentVertix + 1
                    
                        end do
                                                   
                        !Close polygon
                        Polygon%VerticesF(CurrentVertix)%X  = Polygon%VerticesF(1)%X
                        Polygon%VerticesF(CurrentVertix)%Y  = Polygon%VerticesF(1)%Y

                        call SetLimits(Polygon)

                        deallocate(PointCoordinates)
                        nullify   (Polygon)
                    else

                        call Block_Unlock(PolygonsFile, ClientNumber, STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)stop 'NewPolygon - ModuleDrawing - ERR04'
                       

                        exit do1    !No more blocks

                    end if if2

                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1

                    if(STAT_CALL .ne. SUCCESS_)stop 'NewPolygon - ModuleDrawing - ERR05'

                end if if1
            end do do1
        
            call KillEnterData(PolygonsFile, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'NewPolygon - ModuleDrawing - ERR06'

    end subroutine NewPolygon
    !--------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------
    subroutine AddPolygon (Polygons, ObjPolygon)

        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer                   :: Polygons
        type (T_Polygon), pointer                   :: ObjPolygon

        !Local-----------------------------------------------------------------
        type (T_Polygon), pointer                   :: NewPolygon
        type (T_Polygon), pointer                   :: PreviousPolygon
        integer, save                               :: NextID = 1
        
        !Allocates new instance
        allocate (NewPolygon)
        nullify  (NewPolygon%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Polygons)) then
            Polygons             => NewPolygon
            ObjPolygon           => NewPolygon
        else
            PreviousPolygon      => Polygons
            ObjPolygon           => Polygons%Next
            do while (associated(ObjPolygon))
                PreviousPolygon  => ObjPolygon
                ObjPolygon       => ObjPolygon%Next
            enddo
            ObjPolygon           => NewPolygon
            PreviousPolygon%Next => NewPolygon
        endif

        !Attributes ID
        ObjPolygon%ID            = NextID
        NextID                   = NextID + 1

    end subroutine AddPolygon

    !--------------------------------------------------------------------------
   
    subroutine NewLine(Lines, LinesFileName)
        
        !Arguments-------------------------------------------------------------
        type(T_Lines),                 pointer      :: Lines
        character(len=*),  intent(IN)               :: LinesFileName

        !Local-----------------------------------------------------------------
        integer                                     :: ObjEnterData = 0
        logical                                     :: BlockFound
        integer                                     :: ClientNumber, STAT_CALL
        integer                                     :: iflag, FromBlock
        integer                                     :: ret, StartLine, EndLine
        integer                                     :: iPoint, CurrentLineNumber
        type(T_Lines), pointer                      :: CurrLine
        real, dimension(:), pointer                 :: PointCoordinates

        nullify(CurrLine)

        call ConstructEnterData(ObjEnterData, LinesFileName, STAT = ret) 
        if (ret == 0) then

            !ColorStuff
            BlockFound = .true.
            do while (BlockFound)

                call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,              &
                                            "<begin_line>", "<end_line>",            &
                                            BlockFound,StartLine, EndLine,           &
                                            STAT = STAT_CALL)
                if (BlockFound) then
                    
                    call AddLine (Lines, CurrLine)

                    CurrLine%nNodes = EndLine - StartLine - 1
                
                    iPoint = 1

                    allocate(CurrLine%X(1:CurrLine%nNodes))
                    allocate(CurrLine%Y(1:CurrLine%nNodes))
                    allocate(PointCoordinates(1:2))

                    do CurrentLineNumber = StartLine + 1 , EndLine - 1 
                    
                        call GetData(PointCoordinates,                               &
                                     ObjEnterData,                                   &
                                     flag = iflag,                                   &
                                     SearchType = FromBlock,                         &
                                     Buffer_Line = CurrentLineNumber,                &
                                     STAT        = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_ .and.                             &
                           STAT_CALL .ne. SIZE_ERR_) stop 'NewLine - ModuleDrawing - ERR03'

                        CurrLine%X(iPoint) = PointCoordinates(1)
                        CurrLine%Y(iPoint) = PointCoordinates(2)

                        iPoint = iPoint + 1

                    end do

                endif
            enddo
        end if

        call KillEnterData(ObjEnterData, ret)

    end subroutine NewLine   
   
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine AddLine (Lines, ObjLine)

        !Arguments-------------------------------------------------------------
        type (T_Lines), pointer                   :: Lines
        type (T_Lines), pointer                   :: ObjLine

        !Local-----------------------------------------------------------------
        type (T_Lines), pointer                   :: NewLine
        type (T_Lines), pointer                   :: PreviousLine
        integer, save                            :: NextID = 1
        
        !Allocates new instance
        allocate (NewLine)
        nullify  (NewLine%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Lines)) then
            Lines             => NewLine
            ObjLine           => NewLine
        else
            PreviousLine      => Lines
            ObjLine           => Lines%Next
            do while (associated(ObjLine))
                PreviousLine  => ObjLine
                ObjLine       => ObjLine%Next
            enddo
            ObjLine           => NewLine
            PreviousLine%Next => NewLine
        endif

        !Attributes ID
        ObjLine%ID            = NextID
        NextID                = NextID + 1

    end subroutine AddLine

    !--------------------------------------------------------------------------
    subroutine NewXYZPoint(XYZPoints, XYZPointsFileName, DepthReferential, AddFactor, FileFormat)
        
        !Arguments-------------------------------------------------------------
        type(T_XYZPoints),               pointer    :: XYZPoints
        character(len=*),  intent(IN)               :: XYZPointsFileName
        real,    optional, intent(IN)               :: DepthReferential, AddFactor
        integer, optional, intent(IN)               :: FileFormat

        !Local-----------------------------------------------------------------
        type(T_XYZPoints),               pointer    :: XYZPoint
        integer                                     :: XYZPointFile = 0
        integer                                     :: ClientNumber
        integer                                     :: StartLine, EndLine
        integer                                     :: FromBlock, CurrentLineNumber
        character(len = line_length)                :: FullBufferLine
        logical                                     :: BlockFound = OFF
        real, dimension(:),              pointer    :: PointCoordinates
        integer                                     :: i, STAT_CALL
        real                                        :: DepthReferential_, AddFactor_

        !Begin-----------------------------------------------------------------

        call ConstructEnterData(XYZPointFile, XYZPointsFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'NewXYZPoint - ModuleDrawing - ERR01'

        if (present(DepthReferential)) then
            DepthReferential_ = DepthReferential
        else
            DepthReferential_ = 1.
        endif

        if (present(AddFactor)) then
            AddFactor_ = AddFactor
        else
            AddFactor_ = 0.
        endif
        
do1 :       do
                call ExtractBlockFromBuffer(XYZPointFile,                                &
                                            ClientNumber,                                &
                                            '<begin_xyz>',                               &
                                            '<end_xyz>',                                 &
                                            BlockFound,                                  &
                                            STAT = STAT_CALL)

if1 :           if      (STAT_CALL .EQ. SUCCESS_ ) then    

if2 :               if (BlockFound) then
                        
                        call GetExtractType(FromBlock = FromBlock)
                        
                        call GetBlockSize(XYZPointFile,                                  &
                                          ClientNumber,                                  &
                                          StartLine,                                     &
                                          EndLine,                                       &
                                          FromBlock,                                     &
                                          STAT = STAT_CALL)

                        !ATTENTION!!! there is a bug in module EnterData regarding
                        !             the STAT variable in sub GetBlockSize, see also
                        !             modules Lagrangian and TimeSerie

                        
                        call Add(XYZPoints, XYZPoint)
                        
                        if (present(FileFormat)) then
                            XYZPoint%FileFormat = FileFormat
                        else
                            XYZPoint%FileFormat = TypeX_Y_Z
                        endif

                        if      (XYZPoint%FileFormat == TypeX_Y_Z)     then
                            allocate(PointCoordinates(1:3))
                        else if (XYZPoint%FileFormat == TypeX_Y_Z_P) then
                            allocate(PointCoordinates(1:4))
                        else if (XYZPoint%FileFormat == TypeT_X_Y_Z_P) then
                            allocate(PointCoordinates(1:5))
                        endif

                        XYZPoint%Count = EndLine - StartLine - 1
                        
                        allocate(XYZPoint%X     (1:(XYZPoint%Count))) 
                        allocate(XYZPoint%Y     (1:(XYZPoint%Count)))
                        allocate(XYZPoint%Z     (1:(XYZPoint%Count))) 
                        allocate(XYZPoint%Inside(1:(XYZPoint%Count))) 

                        if      (XYZPoint%FileFormat == TypeX_Y_Z_P) then
                            allocate(XYZPoint%P(1:(XYZPoint%Count)))                            
                        else if (XYZPoint%FileFormat == TypeT_X_Y_Z_P) then
                            allocate(XYZPoint%T(1:(XYZPoint%Count)))
                            allocate(XYZPoint%P(1:(XYZPoint%Count)))                            
                        endif
                        
                        XYZPoint%Inside = .false.

                        i = 0

                        do CurrentLineNumber = StartLine + 1 , EndLine - 1 
                            
                            i = i + 1

                            call GetFullBufferLine(XYZPointFile,                         &
                                                   CurrentLineNumber,                    &
                                                   FullBufferLine,                       &
                                                   STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_)stop 'NewXYZPoint - ModuleDrawing - ERR02'
                        
                            call GetExtractType(FromBlock = FromBlock)

                            !Read XYZ even if there is a comment after Z
                            call GetData(PointCoordinates,                               &
                                         XYZPointFile,                                   &
                                         FromBlock,                                      &
                                         Buffer_Line = CurrentLineNumber,                &
                                         STAT        = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_ .and.                             &
                               STAT_CALL .ne. SIZE_ERR_) then
                                
                                write(*,*)'Error in Line Number = ', CurrentLineNumber
                                write(*,*)'File name            = ', trim(XYZPointsFileName)
                                stop 'NewXYZPoint - ModuleDrawing - ERR03'

                            end if
                            
                            if (XYZPoint%FileFormat == TypeX_Y_Z) then
                                XYZPoint%X(i)  = PointCoordinates(1)
                                XYZPoint%Y(i)  = PointCoordinates(2)
                                XYZPoint%Z(i)  = PointCoordinates(3)                            
                            else if (XYZPoint%FileFormat == TypeX_Y_Z_P) then
                                XYZPoint%X(i)  = PointCoordinates(1)
                                XYZPoint%Y(i)  = PointCoordinates(2)
                                XYZPoint%Z(i)  = PointCoordinates(3)                            
                                XYZPoint%P(i)  = PointCoordinates(4)                            
                            else if (XYZPoint%FileFormat == TypeT_X_Y_Z_P) then
                                XYZPoint%T(i)  = PointCoordinates(1)                            
                                XYZPoint%X(i)  = PointCoordinates(2)
                                XYZPoint%Y(i)  = PointCoordinates(3)
                                XYZPoint%Z(i)  = PointCoordinates(4)                            
                                XYZPoint%P(i)  = PointCoordinates(5)                            
                            endif
                            
                            XYZPoint%Z(i)  = XYZPoint%Z(i) * DepthReferential_ + AddFactor_

                        end do
                        
                        call SetLimits(XYZPoint)
                        
                        deallocate(PointCoordinates)
                        nullify   (XYZPoint)
                    else

                        call Block_Unlock(XYZPointFile, ClientNumber, STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)stop 'NewXYZPoint - ModuleDrawing - ERR04'
                       

                        exit do1    !No more blocks

                    end if if2

                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1

                    if(STAT_CALL .ne. SUCCESS_)stop 'NewXYZPoint - ModuleDrawing - ERR05'

                end if if1
            end do do1
        
            call KillEnterData(XYZPointFile, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'NewXYZPoint - Polygon - ModuleDrawing'

    end subroutine NewXYZPoint
    !--------------------------------------------------------------------------
    
   
    !--------------------------------------------------------------------------
    subroutine NewXYZPoint_V2(XYZPoints, XVector, YVector, ZVector)
        
        !Arguments-------------------------------------------------------------
        type(T_XYZPoints),               pointer    :: XYZPoints
        real,   dimension(:),            pointer    :: XVector, YVector, ZVector

        !Local-----------------------------------------------------------------
        type(T_XYZPoints),               pointer    :: XYZPoint
        integer                                     :: i
 
        !Begin-----------------------------------------------------------------

 
        
        call Add(XYZPoints, XYZPoint)

        XYZPoint%Count = size(XVector)
        
        allocate(XYZPoint%X     (1:(XYZPoint%Count))) 
        allocate(XYZPoint%Y     (1:(XYZPoint%Count)))
        allocate(XYZPoint%Z     (1:(XYZPoint%Count))) 
        allocate(XYZPoint%Inside(1:(XYZPoint%Count))) 
        
        XYZPoint%FileFormat = TypeX_Y_Z
       
        XYZPoint%Inside = .false.

         do i = 1 , XYZPoint%Count
           
            XYZPoint%X(i)  = XVector(i)
            XYZPoint%Y(i)  = YVector(i)
            XYZPoint%Z(i)  = ZVector(i)

        end do
        
        call SetLimits(XYZPoint)
        

        nullify   (XYZPoint)

    end subroutine NewXYZPoint_V2
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine NewXYZPoint_V3(XYZPoints, XVector, YVector, ZVector, PVector)
        
        !Arguments-------------------------------------------------------------
        type(T_XYZPoints),               pointer    :: XYZPoints
        real,   dimension(:),            pointer    :: XVector, YVector, ZVector, PVector

        !Local-----------------------------------------------------------------
        type(T_XYZPoints),               pointer    :: XYZPoint
        integer                                     :: i
 
        !Begin-----------------------------------------------------------------

 
        
        call Add(XYZPoints, XYZPoint)

        XYZPoint%Count = size(XVector)
        
        allocate(XYZPoint%X     (1:(XYZPoint%Count))) 
        allocate(XYZPoint%Y     (1:(XYZPoint%Count)))
        allocate(XYZPoint%Z     (1:(XYZPoint%Count))) 
        allocate(XYZPoint%Inside(1:(XYZPoint%Count))) 
        allocate(XYZPoint%P     (1:(XYZPoint%Count)))                 
        
        XYZPoint%FileFormat = TypeX_Y_Z_P
       
        XYZPoint%Inside = .false.

         do i = 1 , XYZPoint%Count
            XYZPoint%X(i)  = XVector(i)
            XYZPoint%Y(i)  = YVector(i)
            XYZPoint%Z(i)  = ZVector(i)
            XYZPoint%P(i)  = PVector(i)
        end do
        
        call SetLimits(XYZPoint)
        
        nullify   (XYZPoint)

    end subroutine NewXYZPoint_V3
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    subroutine NewXYZPoint_V4(XYZPoints, TVector, XVector, YVector, ZVector, PVector)
        
        !Arguments-------------------------------------------------------------
        type(T_XYZPoints),               pointer    :: XYZPoints
        real,   dimension(:),            pointer    :: TVector, XVector, YVector, ZVector, PVector

        !Local-----------------------------------------------------------------
        type(T_XYZPoints),               pointer    :: XYZPoint
        integer                                     :: i
 
        !Begin-----------------------------------------------------------------

 
        
        call Add(XYZPoints, XYZPoint)

        XYZPoint%Count = size(XVector)
        
        allocate(XYZPoint%X     (1:(XYZPoint%Count))) 
        allocate(XYZPoint%Y     (1:(XYZPoint%Count)))
        allocate(XYZPoint%Z     (1:(XYZPoint%Count))) 
        allocate(XYZPoint%Inside(1:(XYZPoint%Count))) 
        allocate(XYZPoint%T     (1:(XYZPoint%Count)))         
        allocate(XYZPoint%P     (1:(XYZPoint%Count)))                 
        
        XYZPoint%FileFormat = TypeT_X_Y_Z_P
       
        XYZPoint%Inside = .false.

         do i = 1 , XYZPoint%Count
            XYZPoint%T(i)  = TVector(i)
            XYZPoint%X(i)  = XVector(i)
            XYZPoint%Y(i)  = YVector(i)
            XYZPoint%Z(i)  = ZVector(i)
            XYZPoint%P(i)  = PVector(i)
        end do
        
        call SetLimits(XYZPoint)
        
        nullify   (XYZPoint)

    end subroutine NewXYZPoint_V4
    !--------------------------------------------------------------------------    
    !--------------------------------------------------------------------------
    subroutine AddXYZPoint (XYZPoints, ObjXYZPoint)

        !Arguments-------------------------------------------------------------
        type (T_XYZPoints), pointer                 :: XYZPoints
        type (T_XYZPoints), pointer                 :: ObjXYZPoint

        !Local-----------------------------------------------------------------
        type (T_XYZPoints), pointer                 :: NewXYZPoint
        type (T_XYZPoints), pointer                 :: PreviousXYZPoint
        integer, save                               :: NextID = 1
        
        !Allocates new instance
        allocate (NewXYZPoint)
        nullify  (NewXYZPoint%Next)

        if (.not. associated(XYZPoints)) then
            XYZPoints            => NewXYZPoint
            ObjXYZPoint          => NewXYZPoint
        else
            PreviousXYZPoint     => XYZPoints
            ObjXYZPoint          => XYZPoints%Next
            do while (associated(ObjXYZPoint))
                PreviousXYZPoint => ObjXYZPoint
                ObjXYZPoint      => ObjXYZPoint%Next
            enddo
            ObjXYZPoint          => NewXYZPoint
            PreviousXYZPoint%Next=> NewXYZPoint
        endif


        !Attributes ID
        ObjXYZPoint%ID           = NextID
        NextID                   = NextID + 1

    end subroutine AddXYZPoint


    !--------------------------------------------------------------------------
    
    subroutine SetLimitsPolygon(Polygon)

        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer                   :: Polygon
        
        !Begin-----------------------------------------------------------------

        Polygon%Limits%Left   = minval(Polygon%VerticesF%X)
        Polygon%Limits%Right  = maxval(Polygon%VerticesF%X)
        Polygon%Limits%Bottom = minval(Polygon%VerticesF%Y)
        Polygon%Limits%Top    = maxval(Polygon%VerticesF%Y)

    end subroutine SetLimitsPolygon
    
    !--------------------------------------------------------------------------

    subroutine SetLimitsXYZ(XYZ)

        !Arguments-------------------------------------------------------------
        type (T_XYZPoints), pointer                 :: XYZ
        
        !Begin-----------------------------------------------------------------

        XYZ%Limits%Left     = minval(XYZ%X)
        XYZ%Limits%Right    = maxval(XYZ%X)
        XYZ%Limits%Bottom   = minval(XYZ%Y)
        XYZ%Limits%Top      = maxval(XYZ%Y)

    end subroutine SetLimitsXYZ

    !--------------------------------------------------------------------------
    
    logical function IsVisible(Polygons, Point)
        
        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer                   :: Polygons
        type (T_PointF),   pointer                  :: Point
        
        !Local-----------------------------------------------------------------
        type (T_Polygon), pointer                   :: CurrentPolygon

        !Begin-----------------------------------------------------------------

        IsVisible = .false.

        CurrentPolygon => Polygons

        do while(associated(CurrentPolygon))

            if(IsPointInsidePolygon(Point, CurrentPolygon))then
                IsVisible = .true.
                return
            end if

            CurrentPolygon => CurrentPolygon%Next
        enddo

    end function IsVisible
    !--------------------------------------------------------------------------
    
    
    
    !--------------------------------------------------------------------------
    !Description    :   The first action performed is an acceleration test. if 
    !                   point is not inside maximum limits of the polygon
    !                   then it is not inside the polygon. If point is inside
    !                   the maximum bounds of the polygon then there is a 
    !                   possibility of being inside the polygon. Thus, it is 
    !                   drawn a "semi-recta" on the X axis to the right of the
    !                   point. If the "semi-recta" intersects the polygon an 
    !                   odd ("ímpar") number of times then the point is 
    !                   inside the polygon. The intersection test is performed 
    !                   for every segment of the polygon. If the point belongs 
    !                   is a vertix of the polygon or belongs to one of the 
    !                   segments that defines the polygon then it is considered 
    !                   to be inside the polygon. There are a few cases in which
    !                   special care needs to be taken, regarding the intersection
    !                   counting. If the intersection point coincides with one of 
    !                   the vertices then one intersection is counted only if it is
    !                   the vertice with the higher Y value. This is the way to 
    !                   count only one intersection once one polygon vertix belongs
    !                   to two segments and otherwise it would be tested twice.
    !                   Another case is when the segment belongs to the "semi-recta".
    !                   In this case no intersection is counted. This in the next 
    !                   segment to be tested the "semi-recta" will intersect the first
    !                   vertice of the segment. If the segments slope is negative then
    !                   the vertice has the higher Y value and one intersection is 
    !                   counted. 
     
    logical function IsPointInsidePolygon(Point, Polygon)

        !Arguments-------------------------------------------------------------
        type (T_PointF),   pointer                  :: Point
        type (T_Polygon), pointer                   :: Polygon

        !Local-----------------------------------------------------------------
        integer                                     :: CurrentVertix
        integer                                     :: NumberOfIntersections
        type(T_Segment),  pointer                   :: Segment
        real                                        :: Slope, IntersectionPointX
        real                                        :: HigherY

        !Begin-----------------------------------------------------------------

        IsPointInsidePolygon  = .false.
        NumberOfIntersections = 0

        !Acceleration test. If point is not inside maximum limits of the polygon
        !then it is not inside the polygon
        if(Point%X < Polygon%Limits%Left .or. Point%X > Polygon%Limits%Right .or.        &
           Point%Y > Polygon%Limits%Top  .or. Point%Y < Polygon%Limits%Bottom)return
        
        allocate(Segment)
        allocate(Segment%StartAt)
        allocate(Segment%EndAt)

do1:    do CurrentVertix = 1, Polygon%Count - 1
            
            !construct segment
            Segment%StartAt%X = Polygon%VerticesF(CurrentVertix)%X
            Segment%StartAt%Y = Polygon%VerticesF(CurrentVertix)%Y
            Segment%EndAt%X   = Polygon%VerticesF(CurrentVertix + 1)%X
            Segment%EndAt%Y   = Polygon%VerticesF(CurrentVertix + 1)%Y
            
            !if point coincides with one of the segments vertices then
            !it in inside the polygon
cd1:        if((Point%X .eq. Segment%StartAt%X .and. Point%Y .eq. Segment%StartAt%Y).or. &
               (Point%X .eq. Segment%EndAt%X   .and. Point%Y .eq. Segment%EndAt%Y  ))then

                IsPointInsidePolygon = .true.
                return
            
            !if point is placed between the segment vertices in the Y axis
            elseif((Point%Y <= Segment%StartAt%Y .and. Point%Y >= Segment%EndAt%Y) .or.  &
                   (Point%Y >= Segment%StartAt%Y .and. Point%Y <= Segment%EndAt%Y))then cd1

                !if segment has a slope
cd3:            if(Segment%StartAt%Y .ne. Segment%EndAt%Y)then
                    
                    !compute slope
                    Slope = (Segment%StartAt%X - Segment%EndAt%X)  /                     &
                            (Segment%StartAt%Y - Segment%EndAt%Y)
                    
                    !compute intersection point X coordinate of the "semi-recta"
                    !with the segment 
                    IntersectionPointX = Segment%StartAt%X - Slope *                     &
                                        (Segment%StartAt%Y - Point%Y)
                    
                    !if point belongs to the segment then it is inside the polygon
cd4:                if(IntersectionPointX .eq. Point%X)then

                        IsPointInsidePolygon = .true.
                        return

                    elseif(IntersectionPointX > Point%X)then cd4
                        
                        !if the intersection point coincides with one of the
                        !vertices then one intersection is counted only if it is
                        !the vertice with the higher Y value 
cd5:                    if(Point%Y .eq. Segment%StartAt%Y .or.                           &
                           Point%Y .eq. Segment%EndAt%Y)then
                            
                            !find higher segment ends higher Y coordinate
                            HigherY = max(Segment%StartAt%Y, Segment%EndAt%Y)

cd6:                        if(Point%Y .eq. HigherY)then
                                NumberOfIntersections = NumberOfIntersections + 1
                            end if cd6

                        else cd5
                            !if the intersection point is placed to the right 
                            !of the point or if point belongs to the segment 
                            !then one intersection is counted
                            NumberOfIntersections = NumberOfIntersections + 1
                        end if cd5
                        
                    end if cd4
                
                elseif(Segment%StartAt%Y .eq. Segment%EndAt%Y)then cd3

cd7:                if(Point%Y .eq. Segment%StartAt%Y)then
                        
                        if(Segment%StartAt%X > Segment%EndAt%X)then

                            if(Point%X <= Segment%StartAt%X .and. Point%X >= Segment%EndAt%X)then
                                !if point belongs to the segment then it is inside the polygon
                                IsPointInsidePolygon = .true.
                                return
                            end if

                        elseif(Segment%StartAt%X < Segment%EndAt%X)then
                            
                            if(Point%X >= Segment%StartAt%X .and. Point%X <= Segment%EndAt%X)then
                                !if point belongs to the segment then it is inside the polygon
                                IsPointInsidePolygon = .true.
                                return
                            end if


                        elseif(Segment%StartAt%X == Segment%EndAt%X)then

                            if(Point%X == Segment%StartAt%X)then
                                !if point belongs to the segment then it is inside the polygon
                                IsPointInsidePolygon = .true.
                                return
                            end if

                        end if    
                                           
                        end if cd7

                end if cd3
                
            end if cd1
        enddo do1

        deallocate(Segment%StartAt)
        deallocate(Segment%EndAt)

        deallocate(Segment)

        !if number of intersections is odd (odd = ímpar) then
        !point is inside the polygon
        IsPointInsidePolygon = IsOdd(NumberOfIntersections)

    end function IsPointInsidePolygon
    
    !--------------------------------------------------------------------------

!  Copyright 2001, softSurfer (www.softsurfer.com)
!  This code may be freely used and modified for any purpose
!  providing that this copyright notice is included with it.
!  SoftSurfer makes no warranty for this code, and cannot be held
!  liable for any real or imagined damage resulting from its use.
!  Users of this code must verify correctness for their application.

!  Assume that classes are already given for the objects:
!     Point and Vector with
!         coordinates {float x, y;}
!         operators for:
!             == to test equality
!             != to test inequality
!             Point  = Point ± Vector
!             Vector = Point - Point
!             Vector = Vector ± Vector
!             Vector = Scalar * Vector    (scalar product)
!     Segment with defining endpoints {Point P0, P1;}
! ===================================================================

!#define TRUE    1
!#define FALSE   0
!#define SMALL_NUM  0.00000001 !  anything that avoids division overflow

!  Note: see the April 2001 Algorithm for the "perp product"
!  http://www.geometryalgorithms.com/Archive/algorithm_0111/algorithm_0111.htm#intersect2D_SegPoly()

!  intersect2D_SegPoly():
!     Input:  S = 2D segment to intersect with the convex polygon
!             n = number of 2D points in the polygon
!             V[] = array of n+1 vertex points with V[n]=V[0]
!       Note: The polygon MUST be convex and
!                 have vertices oriented counterclockwise (ccw).
!             This code does not check for and verify these conditions.
!     Return: FALSE = no intersection
!             TRUE  = a valid intersection segment exists
    logical function Intersect2D_SegPoly(Segment, Polygon)

        !Arguments ------------------------------------------------------------
        type(T_Segment),    pointer                 :: Segment
        type (T_Polygon),   pointer                 :: Polygon

        !Local-----------------------------------------------------------------
        type(T_PointF),     pointer                 :: SegmentDir, AuxP, AuxP1, AuxP2, AuxD
        real(8)                                     :: tE, tL
        real(8)                                     :: t, N, D
        integer                                     :: i
        real, parameter                             :: SMALL_NUM = 1.e-16 !  anything that avoids division overflow

        !Begin-----------------------------------------------------------------

        intersect2D_SegPoly =.true. 

        if (.not. ConvexPolygon(Polygon)) then

            write(*,*) 'Polygon must be convex'
            stop 'Intersect2D_SegPoly - Drawing - ERR10'

        endif

        if (CheckPolygonClockwise(Polygon)) call InvertVerticesOrder(Polygon)


        allocate(AuxP, AuxD, SegmentDir)

        if (EqualPoint(Segment%StartAt,Segment%EndAt)) then !  the segment S is a single point
                                                   !  test for inclusion of S.P0 in the polygon
            intersect2D_SegPoly = IsPointInsidePolygon(Segment%StartAt, Polygon)
            return
        endif

        SegmentDir = LessPoint(Segment%EndAt,Segment%StartAt)   !  the segment direction vector

        tE = 0; tL = 1;
    
        do i=1, Polygon%Count - 1 !process polygon edge V[i]V[i+1]

            AuxP2 => Polygon%VerticesF(i+1)
            AuxP1 => Polygon%VerticesF(i)
            AuxP =   LessPoint(AuxP2, AuxP1)
            AuxD =   LessPoint(Segment%StartAt, AuxP1)

            N = PerpProduct2D(AuxP,AuxD);!  = -dot(ne, S.P0-V[i])

            D = - PerpProduct2D(AuxP , SegmentDir);      !  = dot(ne, dS)

            if (abs(D) < SMALL_NUM) then      !  S is nearly parallel to this edge
                if (N < 0) then              !  P0 is outside this edge, so
                    intersect2D_SegPoly = .false.
                    exit               !  S is outside the polygon
                endif                  !  S cannot cross this edge, so
            endif                      !  ignore this edge

            t = N / D;
            if (D < 0) then        !  segment S is entering across this edge
                if (t > tE) then   !  new max tE
                    tE = t;
                    if (tE > tL)  then !  S enters after leaving polygon
                        intersect2D_SegPoly = .false.
                        exit 
                    endif
                endif

            else                   !  segment S is leaving across this edge
                if (t < tL) then   !  new min tL
                    tL = t
                    if (tL < tE) then !  S leaves before entering polygon
                        intersect2D_SegPoly = .false.
                        exit 
                    endif
                endif
            endif


                !  tE <= tL implies that there is a valid intersection subsegment
                !IS->P0 = S.P0 + tE * dS;   !  = P(tE) = point where S enters polygon
                !IS->P1 = S.P0 + tL * dS;   !  = P(tL) = point where S leaves polygon
        enddo

        deallocate(AuxP, AuxD, SegmentDir)
        nullify   (AuxP, AuxD, SegmentDir, AuxP1, AuxP2)

    end function intersect2D_SegPoly 

!---------------------------------------------------------------------------------

!This function returns if the a polygon is convex or concave 
!
!          CONVEX == .true.
!          CONCAVE == .false.
!   It is assumed that the polygon is simple
!   (does not intersect itself or have holes)

logical function ConvexPolygon(Polygon)

        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer           :: Polygon
        
        !Local-----------------------------------------------------------------
        real                                :: Aux, Aux2, x1, x2, x3, y1, y2, y3
        integer                             :: cv, i1, i2, i3

        !Begin-----------------------------------------------------------------
        
        if (Polygon%VerticesF(1)%X /= Polygon%VerticesF(Polygon%Count)%X .or.           &
            Polygon%VerticesF(1)%Y /= Polygon%VerticesF(Polygon%Count)%Y) then
            
            stop 'ConvexPolygon - Drawing - ERR10'

        endif

        ConvexPolygon = .true.

        Aux = 0.

        do cv = 1, Polygon%Count-1
            if (cv==1) then 
                i1 = Polygon%Count
            else
                i1 = cv-1
            endif
            i2 = cv
            i3 = cv + 1
            !(xi - xi-1) * (yi+1 - yi) - (yi - yi-1) * (xi+1 - xi)
            x1 = Polygon%VerticesF(i1)%X
            x2 = Polygon%VerticesF(i2)%X
            x3 = Polygon%VerticesF(i3)%X
            y1 = Polygon%VerticesF(i1)%Y
            y2 = Polygon%VerticesF(i2)%Y
            y3 = Polygon%VerticesF(i3)%Y

            Aux2 = (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2)
            
            if ( cv>1 .and. ((Aux2 > 0. .and. Aux < 0.) .or. (Aux2 < 0. .and. Aux > 0.))) then
                ConvexPolygon = .false.
            endif

            Aux = Aux2

        enddo


end function ConvexPolygon
!-------------------------------------------------------------------------


logical function CheckPolygonClockwise(Polygon)

        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer           :: Polygon
        
        !Local-----------------------------------------------------------------
        real                                :: Aux, x1, x2, x3, y1, y2, y3
        integer                             :: cv, i1, i2, i3

        !Begin-----------------------------------------------------------------
        
        if (Polygon%VerticesF(1)%X /= Polygon%VerticesF(Polygon%Count)%X .or.           &
            Polygon%VerticesF(1)%Y /= Polygon%VerticesF(Polygon%Count)%Y) then
            
            stop 'CheckPolygonClockwise - Drawing - ERR10'

        endif

        Aux = 0.

        do cv = 1, Polygon%Count-1
            if (cv==1) then 
                i1 = Polygon%Count
            else
                i1 = cv-1
            endif
            i2 = cv
            i3 = cv + 1
            !(xi - xi-1) * (yi+1 - yi) - (yi - yi-1) * (xi+1 - xi)
            x1 = Polygon%VerticesF(i1)%X
            x2 = Polygon%VerticesF(i2)%X
            x3 = Polygon%VerticesF(i3)%X
            y1 = Polygon%VerticesF(i1)%Y
            y2 = Polygon%VerticesF(i2)%Y
            y3 = Polygon%VerticesF(i3)%Y
            
            Aux = Aux + (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2)
        enddo

        if (Aux > 0) then
            CheckPolygonClockwise = .false.
        else
            CheckPolygonClockwise = .true.
        endif


    end function CheckPolygonClockwise


!-------------------------------------------------------------------------

!Thius function check if the vertixes of a polygon A are inside of polygon B
logical function VertPolygonInsidePolygon(PolygonA, PolygonB)

        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer           :: PolygonA, PolygonB
        
        !Local-----------------------------------------------------------------
        type (T_PointF),   pointer          :: Point
        integer                             :: cv

        !Begin-----------------------------------------------------------------
        VertPolygonInsidePolygon = .true. 

        do cv = 1, PolygonA%Count
            Point => PolygonA%VerticesF(cv)
            if(.not. IsPointInsidePolygon(Point, PolygonB))then
                VertPolygonInsidePolygon = .false. 
            endif    
        enddo



    end function VertPolygonInsidePolygon


!-------------------------------------------------------------------------

    subroutine InvertVerticesOrder(Polygon)

        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer           :: Polygon
        
        !Local-----------------------------------------------------------------
        type (T_Polygon), pointer           :: PolygonAux
        integer                             :: i, cv

        !Begin-----------------------------------------------------------------
        
        if (Polygon%VerticesF(1)%X /= Polygon%VerticesF(Polygon%Count)%X .or.           &
            Polygon%VerticesF(1)%Y /= Polygon%VerticesF(Polygon%Count)%Y) then
            
            stop 'InvertVerticesOrder - Drawing - ERR10'

        endif

        allocate(PolygonAux)
        
        allocate(PolygonAux%VerticesF(Polygon%Count))
        
        PolygonAux%Count = Polygon%Count


        do cv = 2, Polygon%Count-1

            i = Polygon%Count - cv + 1

            PolygonAux%VerticesF(i)%X = Polygon%VerticesF(cv)%X
            PolygonAux%VerticesF(i)%Y = Polygon%VerticesF(cv)%Y

        enddo

        do cv = 2, Polygon%Count-1
            Polygon%VerticesF(cv)%X = PolygonAux%VerticesF(cv)%X
            Polygon%VerticesF(cv)%Y = PolygonAux%VerticesF(cv)%Y
        enddo

        deallocate(PolygonAux%VerticesF)
        nullify   (PolygonAux%VerticesF)
        
        deallocate(PolygonAux)
        nullify   (PolygonAux)
        


    end subroutine InvertVerticesOrder


!-------------------------------------------------------------------------
    logical function IsPointInsideCircle(Point, CircleOrigin, Radius)

        !Arguments-------------------------------------------------------------
        type (T_PointF),   pointer                  :: Point, CircleOrigin
        real                                        :: Radius

        !Local-----------------------------------------------------------------
        !type(T_Segment),  pointer                   :: Segment
        real                                        :: Distance

        !Begin-----------------------------------------------------------------

        IsPointInsideCircle  = .false.

        !allocate(Segment)

        !Segment%StartAt%X = CircleOrigin%X
        !Segment%StartAt%Y = CircleOrigin%Y
        !Segment%EndAt%X   = Point%X
        !Segment%EndAt%Y   = Point%Y

        Distance = sqrt((CircleOrigin%X-Point%X)**2.+(CircleOrigin%Y-Point%Y)**2.)

        if(Distance <= Radius)IsPointInsideCircle = .true.

        !deallocate(Segment)



    end function IsPointInsideCircle

    
    !--------------------------------------------------------------------------
    
    
    subroutine WriteItemPolygon(Polygon, FilePath)

        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer           :: Polygon
        character(len=*), intent(IN)        :: FilePath
        
        !Local-----------------------------------------------------------------
        integer                             :: UnitID, STAT_CALL, CurrentVertix

        !Begin-----------------------------------------------------------------
        
        call UnitsManager(UnitID, OPEN_FILE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteItemPolygon - ModuleDrawing - ERR01'

    

        open(unit=UnitID, status = "unknown", file = FilePath)
        rewind(UnitID)


        write(UnitID, *)"<begin_polygon>"

        do CurrentVertix = 1, Polygon%Count
            write(UnitID, *)Polygon%VerticesF(CurrentVertix)%X, Polygon%VerticesF(CurrentVertix)%Y
        enddo
        write(UnitID, *)"<end_polygon>"


        call UnitsManager(UnitID, CLOSE_FILE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteItemPolygon - ModuleDrawing - ERR02'


    end subroutine WriteItemPolygon
    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------
    subroutine WriteItemXYZ(XYZ, FilePath)

        !Arguments-------------------------------------------------------------
        type (T_XYZPoints), pointer         :: XYZ
        character(len=*),   intent(IN)      :: FilePath
        
        !Local-----------------------------------------------------------------
        integer                             :: UnitID, STAT_CALL, CurrentPoint

        !Begin-----------------------------------------------------------------
        
        call UnitsManager(UnitID, OPEN_FILE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteItemXYZ - ModuleDrawing - ERR01'

        open(unit=UnitID, status = "unknown", file = FilePath)
        rewind(UnitID)

        write(UnitID, *)"<begin_xyz>"

        do CurrentPoint = 1, XYZ%Count
            write(UnitID, *)XYZ%X(CurrentPoint), XYZ%Y(CurrentPoint), XYZ%Z(CurrentPoint)
        enddo
        
        write(UnitID, *)"<end_xyz>"


        call UnitsManager(UnitID, CLOSE_FILE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteItemXYZ - ModuleDrawing - ERR02'


    end subroutine WriteItemXYZ

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine WriteItemXYZ_v2(XYZArray, TotalXYZ, FilePath)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer                :: XYZArray
        integer,                         intent(IN)  :: TotalXYZ
        character(len=*),                intent(IN)  :: FilePath

        
        !Local-----------------------------------------------------------------
        integer                             :: UnitID, STAT_CALL, CurrentPoint

        !Begin-----------------------------------------------------------------
        
        call UnitsManager(UnitID, OPEN_FILE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteItemXYZ - ModuleDrawing - ERR01'

        open(unit=UnitID, status = "unknown", file = FilePath)
        rewind(UnitID)

        write(UnitID, *)"<begin_xyz>"

        do CurrentPoint = 1, TotalXYZ
            write(UnitID, *)XYZArray(CurrentPoint,1), XYZArray(CurrentPoint,2), XYZArray(CurrentPoint,3)
        enddo
        
        write(UnitID, *)"<end_xyz>"


        call UnitsManager(UnitID, CLOSE_FILE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteItemXYZ - ModuleDrawing - ERR02'


    end subroutine WriteItemXYZ_v2

    !--------------------------------------------------------------------------
    
    !Description: This routine was created to compute the minimum distance from a 
    !point to a polygon or group of polygons in any given directions. The calculation
    !is made for points that do not belong to the polygon or group of polygons.
    !The aim of this routine was to compute distances to land.
    !The routine inputs are the point in question, the polygon or group of polygons
    !and a angle list (in degrees) with the desired directions.
    !The subroutine returns a list of the minimum distances from the point to the 
    !border of the polygons, in the selected directions.
    !The angle list and the result list (MinDistanceFromPoint) must be allocated in 
    !the module that calls this routine having both the same dimension.
    !
    !Method: It is computed the intersection point between the tangent to the direction 
    !and the tangent to the polygon segment. If the intersection point is between the 
    !polygon segment limits and is along the right direction then distances are evaluated  
 
    subroutine PointDistanceToPolygon (Point, Polygons, AngleList,                &
                                             MinDistanceFromPoint)
        
        !Arguments-------------------------------------------------------------
        type (T_PointF),    pointer                 :: Point
        type (T_Polygon),   pointer                 :: Polygons
        real, dimension(:), pointer                 :: AngleList
        real, dimension(:), pointer                 :: MinDistanceFromPoint      

        !Local-----------------------------------------------------------------
        integer                                     :: CurrentVertix
        type(T_Segment),  pointer                   :: Segment
        real                                        :: SegmentSlope, DirectionSlope  
        real                                        :: IntersectionPointX, IntersectionPointY
        real                                        :: Xleft, Xright, Ytop, Ybottom
        integer                                     :: Direction, EndDirection, LastVertix
        real                                        :: Angle
        real                                        :: DirectionNewPointX, DirectionNewPointY
        real                                        :: XDistance, YDistance, DistanceFromPoint
        real                                        :: DirectionX, DirectionY
        type (T_Polygon), pointer                   :: CurrentPolygon
        !Begin-----------------------------------------------------------------
        
        !How many directions to compute
        EndDirection = size(AngleList)
        
        !Nullify Distances
        do Direction = 1, EndDirection 
            MinDistanceFromPoint(Direction)=0.
        enddo
      
        !If point is not inside polygons collection
        if(.not.IsVisible(Polygons, Point)) then

            allocate(Segment)
                      
            do Direction = 1, EndDirection
                
                MinDistanceFromPoint(Direction)=999999999999.
                Angle = AngleList(Direction) * Pi / 180.
                
                !Creates a point along the selected direction to produce a direction 
                !segment
                DirectionNewPointX = Point%X + cos (Angle)
                DirectionNewPointY = Point%Y + sin (Angle)
                
                CurrentPolygon => Polygons

                do while(associated(CurrentPolygon))      

                    !Last vertix has always to be filled in AddPolygon routine
                    !If polygon(s) input is(are) already closed last vertix will
                    !be the same as previous and the evaluation at this vertix
                    !will not be made (to avoid equal consecutive points in the polygons)
                    if(CurrentPolygon%VerticesF(CurrentPolygon%Count)%X.eq.&
                        CurrentPolygon%VerticesF(CurrentPolygon%Count-1)%X.and.&
                        CurrentPolygon%VerticesF(CurrentPolygon%Count)%Y.eq.&
                        CurrentPolygon%VerticesF(CurrentPolygon%Count-1)%Y) then

                        LastVertix = CurrentPolygon%Count - 2
                    else
                        LastVertix = CurrentPolygon%Count - 1
                    endif
                    
                    !Run across polygon segments
                    do CurrentVertix = 1, LastVertix

                        Segment%StartAt%X = CurrentPolygon%VerticesF(CurrentVertix)%X
                        Segment%StartAt%Y = CurrentPolygon%VerticesF(CurrentVertix)%Y
                        Segment%EndAt%X   = CurrentPolygon%VerticesF(CurrentVertix + 1)%X
                        Segment%EndAt%Y   = CurrentPolygon%VerticesF(CurrentVertix + 1)%Y

                        Xright= max(Segment%StartAt%X, Segment%EndAt%X)
                        Xleft= min(Segment%StartAt%X, Segment%EndAt%X)
                        Ytop= max(Segment%StartAt%Y, Segment%EndAt%Y)
                        Ybottom= min(Segment%StartAt%Y, Segment%EndAt%Y)

                        !if segments have slope 
                        if (Segment%StartAt%X.ne.Segment%EndAt%X.and.                   &
                            DirectionNewPointX.ne.Point%X) then

                            SegmentSlope = (Segment%StartAt%Y - Segment%EndAt%Y)  /     &
                                           (Segment%StartAt%X - Segment%EndAt%X)

                            DirectionSlope = (DirectionNewPointY - Point%Y)  /          &
                                             (DirectionNewPointX - Point%X)

                            if (SegmentSlope.ne.DirectionSlope) then
                                                   
                                !Solving tangent equations
                                IntersectionPointX = (Segment%StartAt%Y - SegmentSlope *&
                                                      Segment%StartAt%X - Point%Y +     &
                                                      DirectionSlope * Point%X) /       &
                                                      (DirectionSlope - SegmentSlope)

                                IntersectionPointY = (DirectionSlope * (SegmentSlope *  &
                                                      Segment%StartAt%X -               &
                                                      Segment%StartAt%Y) +              &
                                                      SegmentSlope * (Point%Y -         &
                                                      DirectionSlope * Point%X)) /      &
                                                      (SegmentSlope - DirectionSlope)

                            else
                                !Segments are parallel so there is no Intersection Point.
                                !Obbligate the intersection point to go outside the polygon segment limits
                                IntersectionPointX = Xleft - 1
                                IntersectionPointY = YTop + 1

                            endif

                        !if there is one segment with vertical slope,
                        !tangent equations are to be changed (axe change)
                        elseif (Segment%StartAt%X.eq.Segment%EndAt%X.or.                &
                                DirectionNewPointX.eq.Point%X) then
                        
                            !if direction segment and polygon segment are not perpendicular
                            if (Segment%StartAt%Y.ne.Segment%EndAt%Y.and.               &
                                DirectionNewPointY.ne.Point%Y) then

                                SegmentSlope = (Segment%StartAt%X - Segment%EndAt%X)  / &
                                               (Segment%StartAt%Y - Segment%EndAt%Y)

                                DirectionSlope = (DirectionNewPointX - Point%X)  /      &
                                                 (DirectionNewPointY - Point%Y)

                                if (SegmentSlope.ne.DirectionSlope) then

                                    IntersectionPointX = (SegmentSlope * (DirectionSlope *         &
                                                            Point%Y - Point%X) +                   &
                                                            DirectionSlope * (Segment%StartAt%X -  & 
                                                            SegmentSlope * Segment%StartAt%Y)) /   &
                                                            (DirectionSlope - SegmentSlope)

                                    IntersectionPointY = (Point%X - DirectionSlope * Point%Y -     &
                                                          Segment%StartAt%X + SegmentSlope *       &
                                                          Segment%StartAt%Y) /                     &
                                                          (SegmentSlope - DirectionSlope)
                        
                                else
                                    !Segments are parallel so there is no Intersection Point.
                                    !Obbligate the intersection point to go outside the polygon segment limits
                                    IntersectionPointX = Xleft - 1
                                    IntersectionPointY = YTop + 1

                                endif
                    
                            !if direction segment and polygon segment are perpendiclar.
                            !Two cases possible
                            elseif (Segment%StartAt%Y.eq.Segment%EndAt%Y.or.            &
                                    DirectionNewPointY.eq.Point%Y) then

                                if (Segment%StartAt%X.eq.Segment%EndAt%X.and.           &
                                    DirectionNewPointY.eq.Point%Y) then

                                    IntersectionPointX = Segment%StartAt%X
                                    IntersectionPointY = Point%Y

                                elseif (Segment%StartAt%Y.eq.Segment%EndAt%Y.and.       &
                                        DirectionNewPointX.eq.Point%X) then
    
                                    IntersectionPointX = Point%X 
                                    IntersectionPointY = Segment%StartAt%Y

                                else 
                                
                                    write(*,*) 'Possibly one of the polygons defining the wave study area'
                                    write(*,*) 'have consecutive points with the same X and Y coordinates.'
                                    write(*,*) 'Check your polygon files'
                                    stop 'PointDistanceToPolygon - ModuleDrawing - ERR01' 
                                    
                                endif

                            else
                        
                                stop 'PointDistanceToPolygon - ModuleDrawing - ERR10'
                            endif

                        else

                            stop 'PointDistanceToPolygon - ModuleDrawing - ERR20'
                
                        endif                   
                
                        !if the intersection point is in the polygon segment limits "box"
                        if (IntersectionPointX.ge.Xleft.and.                            &
                            IntersectionPointX.le.Xright.and.                           &
                            IntersectionPointY.ge.Ybottom.and.                          &
                            IntersectionPointY.le.Ytop) then
                    
                            !Distances from point to polygon border projected in axes
                            XDistance = IntersectionPointX - Point%X
                            YDistance = IntersectionPointY - Point%Y
                    
                            !Distances used to check if the intersection point obtained is along
                            !the desired direction or along the opossite.
                            DirectionX = DirectionNewPointX - Point%X
                            DirectionY = DirectionNewPointY - Point%Y
                        
                            !The direction tangent can have a intersection point with the polygon
                            !in the right direction or in the opposite. Distances are computed only
                            !if the right direction is encountered.
                            !
                            !Directions in the positive x-axe (From NNE to SSE) and
                            !negative x-axe (From SSW to NNW).
                            if (DirectionX.ne.0.) then
                        
                                if ((DirectionX.gt.0..and.XDistance.ge.0.).or.          &
                                    (DirectionX.lt.0..and.XDistance.le.0.)) then                     
                            
                                    DistanceFromPoint = Sqrt(XDistance**2 + YDistance**2)
                
                                    if (DistanceFromPoint.lt.MinDistanceFromPoint(Direction)) then
                                        MinDistanceFromPoint(Direction) = DistanceFromPoint

                                    endif

                                endif
                    
                            !Directions North and South
                            elseif (DirectionX.eq.0.) then
                        
                                if ((DirectionY.gt.0..and.YDistance.ge.0.).or.          &
                                    (DirectionY.lt.0..and.YDistance.le.0.)) then
                        
                                    DistanceFromPoint = Sqrt(XDistance**2 + YDistance**2)
                
                                    if (DistanceFromPoint.lt.MinDistanceFromPoint(Direction)) then
                                        MinDistanceFromPoint(Direction) = DistanceFromPoint

                                    endif

                                endif

                            endif
                                               
                        endif
            
                    enddo 
                    
                    CurrentPolygon => CurrentPolygon%Next             

                enddo

            enddo
        
            deallocate(Segment)
    
        !else

            !do nothing. if point is inside polygon collection, go to end of routine. 
            !Distances were reset to zero in the beggining of routine.
            
        endif
        

    end subroutine PointDistanceToPolygon

    !--------------------------------------------------------------------------

    real function DotProduct2D(PointA, PointB) 
        
        !Arguments-------------------------------------------------------------
        type (T_PointF), pointer  :: PointA, PointB

        !Begin-----------------------------------------------------------------
        DotProduct2D = PointA%X * PointB%X + PointA%Y * PointB%Y

    end function DotProduct2D


    !--------------------------------------------------------------------------

    real function PerpProduct2D(PointA, PointB) 
        
        !Arguments-------------------------------------------------------------
        type (T_PointF), pointer :: PointA, PointB

        !Begin-----------------------------------------------------------------
        PerpProduct2D = PointA%X * PointB%Y - PointA%Y * PointB%X

    end function PerpProduct2D


    logical function EqualPoint(PointA, PointB) 
        
        !Arguments-------------------------------------------------------------
        type (T_PointF), pointer :: PointA, PointB

        !Begin-----------------------------------------------------------------
        EqualPoint = .false. 
        if (PointA%X == PointB%X .and. PointA%Y == PointB%Y) EqualPoint = .true. 

    end function EqualPoint


    type (T_PointF) function LessPoint(PointA, PointB) 
        
        !Arguments-------------------------------------------------------------
        type (T_PointF), pointer  :: PointA, PointB

        !Begin-----------------------------------------------------------------
        LessPoint%X = PointA%X - PointB%X 
        LessPoint%Y = PointA%Y - PointB%Y

    end function LessPoint




logical function CellInterSectCell (Cell1, Cell2)
!Arguments---------------------------------------
real(8), dimension(2,4) :: Cell1, Cell2
!Local-------------------------------------------
integer                  :: i
!Begin-------------------------------------------

    CellInterSectCell = .false.
    !if cell1 corners inside Cell2
    do i=1,4
        if (PointInsideCell(Cell1(1,i),Cell1(2,i),Cell2)) CellInterSectCell = .true.
    enddo

    !if cell2 corners inside Cell1
    do i=1,4
        if (PointInsideCell(Cell2(1,i),Cell2(2,i),Cell1)) CellInterSectCell = .true.
    enddo


end function CellInterSectCell

logical function PointInsideCell (X, Y, Cell)
!Arguments---------------------------------------
real(8)                  :: X, Y
real(8), dimension(2,4) :: Cell

!Local------------------------------------------
real(8)            :: Xmin, Ymin, Xmax, Ymax
!Begin-------------------------------------------

    PointInsideCell = .false.
    Xmin = Cell(1,1)
    Ymin = Cell(2,1)

    Xmax = Cell(1,3)
    Ymax = Cell(2,3)

    if (X>= Xmin .and. X<=Xmax) then
        if (Y>= Ymin .and. Y<=Ymax) PointInsideCell = .true.
    endif


end function PointInsideCell


subroutine SliceCellIn4 (Xo, Yo, Face, Cells4)

!Arguments----------------------------------------
real(8)                   :: Xo, Yo, Face
real(8), dimension(2,4,4) :: Cells4
!Local--------------------------------------------
integer                   :: i
!Begin--------------------------------------------


    !Lower left cell - lower left corner
    Cells4(1,1,1) = Xo
    Cells4(2,1,1) = Yo    

    
    !Upper left cell - lower left corner
    Cells4(1,1,2) = Xo 
    Cells4(2,1,2) = Yo + Face/2.   

    !Upper right cell - lower left corner
    Cells4(1,1,3) = Xo + Face/2.
    Cells4(2,1,3) = Yo + Face/2.
    
    !Lower right cell- lower left corner
    Cells4(1,1,4) = Xo + Face/2.
    Cells4(2,1,4) = Yo 


    do i=1,4
    
        Cells4(1,2,i) = Cells4(1,1,i)
        Cells4(2,2,i) = Cells4(2,1,i) + Face/2.

        Cells4(1,3,i) = Cells4(1,2,i) + Face/2.
        Cells4(2,3,i) = Cells4(2,2,i) 

        Cells4(1,4,i) = Cells4(1,3,i) 
        Cells4(2,4,i) = Cells4(2,3,i) - Face/2.
    
    enddo
    

end subroutine SliceCellIn4

!This subroutine returns the indexes (WindowOut) of the smallest sub-domain of a curvilinear grid (XX, YY) 
!that comprehends a cartesian window (WindowOut). 
!Returns also a logical if the grid intersects the cartesian window (WindowWithData)
subroutine ArrayPolygonWindow(XX, YY, WindowIn, ILB, IUB, JLB, JUB, WindowOut, WindowWithData)

!Arguments------------------------------------
real(8), dimension(:,:), intent(IN)     :: XX, YY
real(8), dimension(2,4), intent(IN)     :: WindowIn
integer                , intent(IN)     :: ILB, IUB, JLB, JUB
integer, dimension(2,4), intent(OUT)    :: WindowOut
logical,                 intent(OUT)    :: WindowWithData

!Local----------------------------------------
<<<<<<< .mine
type(T_Polygon),          pointer       :: PolygonDomain
type(T_PointF),           pointer       :: Point
integer, dimension(:), allocatable      :: Xj, Yi
integer, dimension(2,100)               :: WindowInIndex
integer                                 :: i, imin, imax, jmin, jmax, count, p, ncells, j
logical, dimension(5)                   :: WindowCornerInside, GridCornerInside
real(8), dimension(2,5)                 :: WindowInAux
real(8), dimension(2,4)                 :: DomainEnvelop
real(8)                                 :: x3, y3, x4, y4, X, Y
logical                                 :: Intersect, NoOverlap

=======
type(T_Polygon),          pointer    :: PolygonDomain
type(T_PointF),           pointer    :: Point
integer, dimension(:), allocatable   :: Xj, Yi
integer, dimension(2,8)     :: WindowInIndex
integer                     :: i, imin, imax, jmin, jmax, count, p, ncells, j
logical, dimension(5)       :: WindowCornerInside, GridCornerInside
real(8), dimension(2,5)     :: WindowInAux
real(8), dimension(2,4)     :: DomainEnvelop
real(8)                     :: x3, y3, x4, y4, X, Y
logical                     :: Intersect
>>>>>>> .r69505
!Begin----------------------------------------

!WindowIn(1,:) -> X
!WindowIn(2,:) -> Y

!WindowInIndex(1,:) -> i
!WindowInIndex(2,:) -> j

   
    WindowInIndex (:,:) = -99

    allocate(Point)
    call CreateDomainPolygon(XX, YY, ILB, IUB, JLB, JUB, PolygonDomain)
    
    DomainEnvelop(1,1) = PolygonDomain%Limits%Left
    DomainEnvelop(2,1) = PolygonDomain%Limits%Bottom
    DomainEnvelop(1,2) = PolygonDomain%Limits%Left
    DomainEnvelop(2,2) = PolygonDomain%Limits%Top
    DomainEnvelop(1,3) = PolygonDomain%Limits%Right
    DomainEnvelop(2,3) = PolygonDomain%Limits%Top
    DomainEnvelop(1,4) = PolygonDomain%Limits%Right
    DomainEnvelop(2,4) = PolygonDomain%Limits%Bottom
    
i1: if (.not. cellIntersectCell(DomainEnvelop, WindowIn)) then
        !if the  window request do not intersect the grid envelop. There is no data to returm    
        WindowWithData =.false.
        
    else i1
        !if there an intersection between the window request and the grid envelop there a high probability of 
        !the window intersect the grid
        count = 0

        if (PointInsideCell(XX(ILB,JLB),YY(ILB,JLB),WindowIn)) then
            count = count + 1                                              
            WindowInIndex(1,count ) = ILB
            WindowInIndex(2,count ) = JLB
            GridCornerInside(1) = .true. 
        else
            GridCornerInside(1) = .false. 
        endif

        if (PointInsideCell(XX(IUB,JLB),YY(IUB,JLB),WindowIn)) then
            count = count + 1                                              
            WindowInIndex(1,count) = IUB-1
            WindowInIndex(2,count) = JLB
            GridCornerInside(2) = .true. 
        else
            GridCornerInside(2) = .false.             
        endif
        
        if (PointInsideCell(XX(IUB,JUB),YY(IUB,JUB),WindowIn)) then
            count = count + 1                                              
            WindowInIndex(1,count) = IUB-1
            WindowInIndex(2,count) = JUB-1
            GridCornerInside(3) = .true. 
        else
            GridCornerInside(3) = .false.             
        endif        

        if (PointInsideCell(XX(ILB,JUB),YY(ILB,JUB),WindowIn)) then
            count = count + 1                                              
            WindowInIndex(1,count) = ILB
            WindowInIndex(2,count) = JUB-1
            GridCornerInside(4) = .true. 
        else
            GridCornerInside(4) = .false.             
        endif       
        
       
        if (count<4) then 
                
            do i=1,4
                Point%X = WindowIn(1,i)
                Point%Y = WindowIn(2,i)

                if (IsPointInsidePolygon(Point, PolygonDomain)) then
                    count = count + 1                                              
                    call LocateCellPolygonsV2(XX,YY, Point, ILB, IUB, JLB, JUB, &
                                              WindowInIndex(1,count),WindowInIndex(2,count))
                    WindowCornerInside(i) = .true.                                      
                else
                    WindowCornerInside(i) = .false.
                endif
            enddo    
        
        endif
        
        if (count==0) then
            !If there is not window corners inside the grid two of the follow conditions are possible:
                !Window  contains the entire grid
                !Window contains part of the grid
            
            !the first option is true if the boundary segments of the window do mnot intersect the grid domain
            Intersect = .false. 

            WindowInAux(:,1:4) = WindowIn(:,:)
            WindowInAux(:,5  ) = WindowIn(:,1)
            
            do i=1,4 
                !Segment 
                x3 = WindowInAux(1,i)
                y3 = WindowInAux(2,i)
                
                x4 = WindowInAux(1,i+1)
                y4 = WindowInAux(2,i+1)
                
                call IntersectionBoundCell(XX,YY, x3, y3, x4, y4, ILB, IUB, JLB, JUB, Intersect, &
                                           WindowInIndex(1,1),WindowInIndex(2,1))
                !Window bound segments intersect the grid domain                                            
                if (Intersect) exit                                           
            enddo
            
            
            if (Intersect) then
                !Window bound segments intersect the grid domain                                            
                !In this case a brute force is requeired. Are verify what domain bound points are inside window
                nCells = 0
                allocate(Xj(1:2*(IUB+1)+2*(JUB+1)),Yi(1:2*(IUB+1)+2*(JUB+1)))
                
                do p=1, PolygonDomain%Count-1
                    X = PolygonDomain%VerticesF(p)%X
                    Y = PolygonDomain%VerticesF(p)%Y
                    if (PointInsideCell(X,Y,WindowIn)) then
                        nCells     = nCells + 1
                        Xj(nCells) = j
                        Yi(nCells) = i
                    endif
                enddo    
                
                imin = minval(Yi(1:nCells))
                imax = maxval(Yi(1:nCells))

                jmin = minval(Xj(1:nCells))
                jmax = maxval(Xj(1:nCells))
                
                
            else
                NoOverlap = .false.
                do i=1,4
                    if (WindowCornerInside(i) .and. GridCornerInside(i)) then
                        NoOverlap      = .true.
                        exit
                    endif
                enddo
                
                if (NoOverlap) then
                    WindowWithData = .false. 
                else
                    !The window contains the entire grid domain
                    imin = ILB
                    imax = IUB - 1

                    jmin = JLB
                    jmax = JUB - 1
                endif
            endif
            
        else        
            !At least there are one window corner inside the grid domain
            if (count<4) then
                !If not all window corners are inside the grid domain it is necesssary to 
                !identify the boundary cell along which the window boundary segments intersect the 
                !grid border
                WindowCornerInside(5) = WindowCornerInside(1)            
        
                WindowInAux(:,1:4) = WindowIn(:,:)
                WindowInAux(:,5  ) = WindowIn(:,1)
                
           
                do i=1,4 
                    !if i outside inside and i+1 outside or vice-versa compute the intersection cell   
                    if (.not.WindowCornerInside(i).eqv.WindowCornerInside(i+1)) then
                        Count = Count + 1
                        !Segment 
                        x3 = WindowInAux(1,i)
                        y3 = WindowInAux(2,i)
                        
                        x4 = WindowInAux(1,i+1)
                        y4 = WindowInAux(2,i+1)
                        
                        call IntersectionBoundCell(XX,YY, x3, y3, x4, y4, ILB, IUB, JLB, JUB, &
                                                   Intersect, WindowInIndex(1,count),WindowInIndex(2,count))
                    endif
                enddo

                GridCornerInside(5) = GridCornerInside(1)
                
<<<<<<< .mine
                call IntersectionBoundCellV2(XX,YY, WindowInAux, ILB, IUB, JLB, JUB, &
                                             Intersect, WindowInIndex, GridCornerInside, Count)
               
=======
                call IntersectionBoundCellV2(XX,YY, WindowInAux, ILB, IUB, JLB, JUB, &
                                           Intersect, WindowInIndex, GridCornerInside, Count)
               
>>>>>>> .r69505
        
            endif
        
            imin = minval(WindowInIndex(1,1:count))
            imax = maxval(WindowInIndex(1,1:count))

            jmin = minval(WindowInIndex(2,1:count))
            jmax = maxval(WindowInIndex(2,1:count))
        
        endif
        
        WindowOut(1,1) = imin
        WindowOut(2,1) = jmin

        WindowOut(1,2) = imax
        WindowOut(2,2) = jmin

        WindowOut(1,3) = imax
        WindowOut(2,3) = jmax

        WindowOut(1,4) = imin
        WindowOut(2,4) = jmax

        WindowWithData = .true.
        
    endif  i1     

end subroutine ArrayPolygonWindow

!Locates in a grid (XX, YY) the cell (IZ, JZ) where a point is (Point)
!There is a similar subroutine in ModuleHorizontalGrid
subroutine LocateCellPolygonsV2(XX, YY, Point, ILB, IUB, JLB, JUB, IZ, JZ)

        !Arguments ---------------------------------------------------------
        real(8), dimension(:,:)              :: XX, YY
        type(T_PointF),           pointer    :: Point
        integer,                 intent(IN ) :: ILB, IUB, JLB, JUB
        integer,                 intent(OUT) :: IZ, JZ

        !Local -------------------------------------------------------------
        integer                              :: ICenter, JCenter, Iupper, ILower,       &
                                                Jupper, JLower
        logical                              :: CellFound
        type(T_Polygon),          pointer    :: Polygon
        integer                              :: i, j, pi
        integer                              :: I1, I2, I3, I4
        integer                              :: J1, J2, J3, J4
        !Begin -------------------------------------------------------------

        allocate(Polygon)
        allocate(Polygon%VerticesF(1: 2*(IUB-ILB+1+JUB-JLB+1)))

        CellFound = .false.

        Iupper = IUB
        ILower = ILB

        Jupper = JUB
        JLower = JLB

        do while (.not. CellFound)

            ICenter = int(real((Iupper - ILower)/ 2.)) + ILower
            JCenter = int(real((Jupper - JLower)/ 2.)) + JLower

            !Construct 4 polygons NW, NE, SW, SE
            do pi = 1, 4
                
                !All polygons are defined anti-Clockwise 
                if      (pi==1) then

                    if (ILower == ICenter) cycle
                    if (JLower == JCenter) cycle
                    
                    I1 = ILower;  J1 = JLower
                    I2 = ILower;  J2 = JCenter
                    I3 = ICenter; J3 = JCenter
                    I4 = ICenter; J4 = JLower

                else if (pi==2) then

                    if (ILower == ICenter) cycle

                    I1 = ILower;  J1 = JCenter;
                    I2 = ILower;  J2 = Jupper;
                    I3 = ICenter; J3 = Jupper;
                    I4 = ICenter; J4 = JCenter;

                else if (pi==3) then

                    if (JLower == JCenter) cycle

                    I1 = ICenter; J1 = JLower;
                    I2 = ICenter; J2 = JCenter;
                    I3 = IUpper ; J3 = JCenter;
                    I4 = IUpper;  J4 = JLower;

                else if (pi==4) then


                    I1 = ICenter; J1 = JCenter;
                    I2 = ICenter; J2 = JUpper;
                    I3 = IUpper ; J3 = JUpper;
                    I4 = IUpper;  J4 = JCenter;

                endif
                
                Polygon%Count = 0.

                do j = J1, J2

                    Polygon%Count = Polygon%Count + 1                   

                    Polygon%VerticesF(Polygon%Count)%X  = XX(I1,j)
                    Polygon%VerticesF(Polygon%Count)%Y  = YY(I1,j)

                end do

                do i = I2, I3

                    Polygon%Count = Polygon%Count + 1                   

                    Polygon%VerticesF(Polygon%Count)%X  = XX(i,J2)
                    Polygon%VerticesF(Polygon%Count)%Y  = YY(i,J2)

                end do

                do j = J3, J4, -1

                    Polygon%Count = Polygon%Count + 1                   

                    Polygon%VerticesF(Polygon%Count)%X  = XX(I3,j)
                    Polygon%VerticesF(Polygon%Count)%Y  = YY(I3,j)
                
                end do

                do i = I4, I1, -1

                    Polygon%Count = Polygon%Count + 1                   

                    Polygon%VerticesF(Polygon%Count)%X  = XX(i,J4)
                    Polygon%VerticesF(Polygon%Count)%Y  = YY(i,J4)
                
                end do

                Polygon%Count = Polygon%Count + 1                   

                !Close polygon
                Polygon%VerticesF(Polygon%Count)%X  = Polygon%VerticesF(1)%X
                Polygon%VerticesF(Polygon%Count)%Y  = Polygon%VerticesF(1)%Y

                Polygon%Limits%Left   = minval(Polygon%VerticesF%X)
                Polygon%Limits%Right  = maxval(Polygon%VerticesF%X)
                Polygon%Limits%Bottom = minval(Polygon%VerticesF%Y)
                Polygon%Limits%Top    = maxval(Polygon%VerticesF%Y)
        
                if (IsPointInsidePolygon(Point, Polygon)) exit

            enddo

            Iupper  = I3
            ILower  = I1

            Jupper  = J2
            JLower  = J1

            !Test if the cell was founded
            if ((Iupper - ILower) == 1 .and. (Jupper - JLower) == 1) CellFound = .true.

        end do 

        IZ = ILower
        JZ = JLower

            
        deallocate(Polygon%VerticesF)        
        nullify   (Polygon%VerticesF)        

        deallocate(Polygon)
        nullify   (Polygon)


    end subroutine LocateCellPolygonsV2

!Returns a polygon (PolygonDomain) that corresponds to the boundary of a grid (XX, YY)
subroutine CreateDomainPolygon(XX, YY, ILB, IUB, JLB, JUB, PolygonDomain)

        !Arguments ---------------------------------------------------------
        real(8), dimension(:,:)              :: XX, YY
        integer,                 intent(IN ) :: ILB, IUB, JLB, JUB
        type(T_Polygon),          pointer    :: PolygonDomain
        !Local -------------------------------------------------------------


        integer                              :: i, j
        !Begin -------------------------------------------------------------

        allocate(PolygonDomain)
        allocate(PolygonDomain%VerticesF(1: 2*(IUB-ILB+1+JUB-JLB+1)+1))

               
        PolygonDomain%Count = 0.

        do j = JLB, JUB

            PolygonDomain%Count = PolygonDomain%Count + 1                   

            PolygonDomain%VerticesF(PolygonDomain%Count)%X  = XX(ILB,j)
            PolygonDomain%VerticesF(PolygonDomain%Count)%Y  = YY(ILB,j)

        end do

        do i = ILB, IUB

            PolygonDomain%Count = PolygonDomain%Count + 1                   

            PolygonDomain%VerticesF(PolygonDomain%Count)%X  = XX(i,JUB)
            PolygonDomain%VerticesF(PolygonDomain%Count)%Y  = YY(i,JUB)

        end do

        do j = JUB, JLB, -1

            PolygonDomain%Count = PolygonDomain%Count + 1                   

            PolygonDomain%VerticesF(PolygonDomain%Count)%X  = XX(IUB,j)
            PolygonDomain%VerticesF(PolygonDomain%Count)%Y  = YY(IUB,j)
        
        end do

        do i = IUB, ILB, -1

            PolygonDomain%Count = PolygonDomain%Count + 1                   

            PolygonDomain%VerticesF(PolygonDomain%Count)%X  = XX(i,JLB)
            PolygonDomain%VerticesF(PolygonDomain%Count)%Y  = YY(i,JLB)
        
        end do

        PolygonDomain%Count = PolygonDomain%Count + 1                   

        !Close PolygonDomain
        PolygonDomain%VerticesF(PolygonDomain%Count)%X  = PolygonDomain%VerticesF(1)%X
        PolygonDomain%VerticesF(PolygonDomain%Count)%Y  = PolygonDomain%VerticesF(1)%Y

        PolygonDomain%Limits%Left   = minval(PolygonDomain%VerticesF%X)
        PolygonDomain%Limits%Right  = maxval(PolygonDomain%VerticesF%X)
        PolygonDomain%Limits%Bottom = minval(PolygonDomain%VerticesF%Y)
        PolygonDomain%Limits%Top    = maxval(PolygonDomain%VerticesF%Y)

    end subroutine CreateDomainPolygon
    
 
 
 
!Checks if a segment (x3, y3, x4, y4) interscts or not (Intersect) the boundary of a grid (XX, YY)
!If intersects returns the boundary cell where the intersection happens (iOut, jOut)
    subroutine IntersectionBoundCell(XX,YY,x3, y3, x4, y4,                             &
                                     ILB, IUB, JLB, JUB, Intersect, iOut,jOut)   
    !Arguments------------------------------------
    real(8), dimension(:,:), intent(IN)     :: XX, YY
    real(8)                                 :: x3, y3, x4, y4
    integer                , intent(IN)     :: ILB, IUB, JLB, JUB
    logical,     optional,   intent(OUT)    :: Intersect
    integer,     optional,   intent(OUT)    :: iOut,jOut


    !Local----------------------------------------
    real(8)                     :: x1,y1,x2,y2
    integer                     :: i, j

    !Begin----------------------------------------
    Intersect = .false.

    
    !West Boundary
    x1 = XX(ILB,JLB)
    y1 = YY(ILB,JLB)
    x2 = XX(IUB,JLB)
    y2 = YY(IUB,JLB)
    
    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
        do i=ILB,IUB-1        
            x1 = XX(i,JLB)
            y1 = YY(i,JLB)
            x2 = XX(i+1,JLB)
            y2 = YY(i+1,JLB)
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                Intersect = .true.
                iOut = i
                jOut = JLB
                exit
            endif                      
        enddo
        
        if (.not. Intersect) stop "IntersectionBoundCell - ERR10"
        
    endif

    !North Boundary
    x1 = XX(IUB,JLB)
    y1 = YY(IUB,JLB)
    x2 = XX(IUB,JUB)
    y2 = YY(IUB,JUB)
    
    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
        do j=JLB,JUB-1        
            x1 = XX(IUB,j  )
            y1 = YY(IUB,j  )
            x2 = XX(IUB,j+1)
            y2 = YY(IUB,j+1)
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                Intersect = .true.
                iOut = IUB-1
                jOut = j
                exit
            endif                      
        enddo
        
        if (.not. Intersect) stop "IntersectionBoundCell - ERR20"
        
    endif

    !East Boundary
    x1 = XX(ILB,JUB)
    y1 = YY(ILB,JUB)
    x2 = XX(IUB,JUB)
    y2 = YY(IUB,JUB)
    
    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
        do i=ILB,IUB-1        
            x1 = XX(i,JUB)
            y1 = YY(i,JUB)
            x2 = XX(i+1,JUB)
            y2 = YY(i+1,JUB)
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                Intersect = .true.
                iOut = i
                jOut = JUB-1
                exit
            endif                      
        enddo
        
        if (.not. Intersect) write(*,*) "IntersectionBoundCell - ERR30"
        
    endif

    !South Boundary
    x1 = XX(ILB,JLB)
    y1 = YY(ILB,JLB)
    x2 = XX(ILB,JUB)
    y2 = YY(ILB,JUB)
    
    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
        do j=JLB,JUB-1        
            x1 = XX(ILB,j  )
            y1 = YY(ILB,j  )
            x2 = XX(ILB,j+1)
            y2 = YY(ILB,j+1)
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                Intersect = .true.
                iOut = ILB
                jOut = j
                exit
            endif                      
        enddo
        
        if (.not. Intersect) stop "IntersectionBoundCell - ERR40"
        
    endif
    
    end subroutine IntersectionBoundCell   

<<<<<<< .mine

subroutine IntersectionBoundCellV2(XX,YY, WindowInAux, ILB, IUB, JLB, JUB, &
                                   Intersect, WindowInIndex, GridCornerInside, Count)


    !Arguments------------------------------------
    real(8), dimension(:,:), intent(IN)     :: XX, YY
    real(8), dimension(2,5)                 :: WindowInAux    
    integer                , intent(IN)     :: ILB, IUB, JLB, JUB
    logical,     optional,   intent(OUT)    :: Intersect
    integer, dimension(2,8)                 :: WindowInIndex    
    logical, dimension(5)                   :: GridCornerInside    
    integer                                 :: Count


    !Local----------------------------------------
    real(8)                     :: x1,y1,x2,y2,x3,y3,x4,y4
    integer                     :: i, j, g, w

    !Begin----------------------------------------
    Intersect = .false.
    

                
g1: do g=1,4
    
        if (GridCornerInside(g).eqv.GridCornerInside(g+1)) cycle
                
w1:     do w=1,4 

            !Segment 
            x3 = WindowInAux(1,w)
            y3 = WindowInAux(2,w)
            
            x4 = WindowInAux(1,w+1)
            y4 = WindowInAux(2,w+1)

        
            !West Boundary
            x1 = XX(ILB,JLB)
            y1 = YY(ILB,JLB)
            x2 = XX(IUB,JLB)
            y2 = YY(IUB,JLB)
            
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                do i=ILB,IUB-1        
                    x1 = XX(i,JLB)
                    y1 = YY(i,JLB)
                    x2 = XX(i+1,JLB)
                    y2 = YY(i+1,JLB)
                    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                        Intersect = .true.

                        Count = Count + 1

                        WindowInIndex(1,Count) = i
                        WindowInIndex(2,Count) = JLB
                    endif                      
                enddo
                
                if (.not. Intersect) stop "IntersectionBoundCell - ERR10"
                
            endif

            !North Boundary
            x1 = XX(IUB,JLB)
            y1 = YY(IUB,JLB)
            x2 = XX(IUB,JUB)
            y2 = YY(IUB,JUB)
            
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                do j=JLB,JUB-1        
                    x1 = XX(IUB,j  )
                    y1 = YY(IUB,j  )
                    x2 = XX(IUB,j+1)
                    y2 = YY(IUB,j+1)
                    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                        Intersect = .true.
                        Count = Count + 1                    
                        WindowInIndex(1,Count) = IUB-1
                        WindowInIndex(2,Count) = j
                    endif                      
                enddo
                
                if (.not. Intersect) stop "IntersectionBoundCell - ERR20"
                
            endif

            !East Boundary
            x1 = XX(ILB,JUB)
            y1 = YY(ILB,JUB)
            x2 = XX(IUB,JUB)
            y2 = YY(IUB,JUB)
            
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                do i=IUB-1,ILB,-1
                    x1 = XX(i,JUB)
                    y1 = YY(i,JUB)
                    x2 = XX(i+1,JUB)
                    y2 = YY(i+1,JUB)
                    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                        Intersect = .true.
                        Count = Count + 1                    
                        WindowInIndex(1,Count) =i
                        WindowInIndex(2,Count) =JUB-1
                    endif                      
                enddo
                
                if (.not. Intersect) write(*,*) "IntersectionBoundCell - ERR30"
                
            endif

            !South Boundary
            x1 = XX(ILB,JLB)
            y1 = YY(ILB,JLB)
            x2 = XX(ILB,JUB)
            y2 = YY(ILB,JUB)
            
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                do j=JUB-1,JLB,-1
                    x1 = XX(ILB,j  )
                    y1 = YY(ILB,j  )
                    x2 = XX(ILB,j+1)
                    y2 = YY(ILB,j+1)
                    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                        Intersect = .true.
                        Count = Count + 1
                        WindowInIndex(1,Count) =ILB
                        WindowInIndex(2,Count) =j
                    endif                      
                enddo
                
                if (.not. Intersect) stop "IntersectionBoundCell - ERR40"
                
            endif
            
        enddo w1
    enddo g1       
    
    end subroutine IntersectionBoundCellV2   


=======

subroutine IntersectionBoundCellV2(XX,YY, WindowInAux, ILB, IUB, JLB, JUB, &
                                   Intersect, WindowInIndex, GridCornerInside, Count)


    !Arguments------------------------------------
    real(8), dimension(:,:), intent(IN)     :: XX, YY
    real(8), dimension(2,5)                 :: WindowInAux    
    integer                , intent(IN)     :: ILB, IUB, JLB, JUB
    logical,     optional,   intent(OUT)    :: Intersect
    integer, dimension(2,8)                 :: WindowInIndex    
    logical, dimension(5)                   :: GridCornerInside    
    integer                                 :: Count


    !Local----------------------------------------
    real(8)                     :: x1,y1,x2,y2,x3,y3,x4,y4
    integer                     :: i, j, g, w

    !Begin----------------------------------------
    Intersect = .false.
    

                
g1: do g=1,4
    
        if (GridCornerInside(g).eqv.GridCornerInside(g+1)) cycle
                
w1:     do w=1,4 

            !Segment 
            x3 = WindowInAux(1,w)
            y3 = WindowInAux(2,w)
            
            x4 = WindowInAux(1,w+1)
            y4 = WindowInAux(2,w+1)

        
            !West Boundary
            x1 = XX(ILB,JLB)
            y1 = YY(ILB,JLB)
            x2 = XX(IUB,JLB)
            y2 = YY(IUB,JLB)
            
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                do i=ILB,IUB-1        
                    x1 = XX(i,JLB)
                    y1 = YY(i,JLB)
                    x2 = XX(i+1,JLB)
                    y2 = YY(i+1,JLB)
                    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                        Intersect = .true.

                        Count = Count + 1

                        WindowInIndex(1,Count) = i
                        WindowInIndex(2,Count) = JLB
                        exit
                    endif                      
                enddo
                
                if (.not. Intersect) stop "IntersectionBoundCell - ERR10"
                
            endif

            !North Boundary
            x1 = XX(IUB,JLB)
            y1 = YY(IUB,JLB)
            x2 = XX(IUB,JUB)
            y2 = YY(IUB,JUB)
            
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                do j=JLB,JUB-1        
                    x1 = XX(IUB,j  )
                    y1 = YY(IUB,j  )
                    x2 = XX(IUB,j+1)
                    y2 = YY(IUB,j+1)
                    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                        Intersect = .true.
                        Count = Count + 1                    
                        WindowInIndex(1,Count) = IUB-1
                        WindowInIndex(2,Count) = j
                        exit
                    endif                      
                enddo
                
                if (.not. Intersect) stop "IntersectionBoundCell - ERR20"
                
            endif

            !East Boundary
            x1 = XX(ILB,JUB)
            y1 = YY(ILB,JUB)
            x2 = XX(IUB,JUB)
            y2 = YY(IUB,JUB)
            
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                do i=ILB,IUB-1        
                    x1 = XX(i,JUB)
                    y1 = YY(i,JUB)
                    x2 = XX(i+1,JUB)
                    y2 = YY(i+1,JUB)
                    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                        Intersect = .true.
                        Count = Count + 1                    
                        WindowInIndex(1,Count) =i
                        WindowInIndex(2,Count) =JUB-1
                        exit
                    endif                      
                enddo
                
                if (.not. Intersect) write(*,*) "IntersectionBoundCell - ERR30"
                
            endif

            !South Boundary
            x1 = XX(ILB,JLB)
            y1 = YY(ILB,JLB)
            x2 = XX(ILB,JUB)
            y2 = YY(ILB,JUB)
            
            if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                do j=JLB,JUB-1        
                    x1 = XX(ILB,j  )
                    y1 = YY(ILB,j  )
                    x2 = XX(ILB,j+1)
                    y2 = YY(ILB,j+1)
                    if (SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)) then
                        Intersect = .true.
                        Count = Count + 1
                        WindowInIndex(1,Count) =ILB
                        WindowInIndex(2,Count) =j
                        exit
                    endif                      
                enddo
                
                if (.not. Intersect) stop "IntersectionBoundCell - ERR40"
                
            endif
            
        enddo w1
    enddo g1       
    
    end subroutine IntersectionBoundCellV2   


>>>>>>> .r69505
 !Checks if a segment (x1,y1,x2,y2) intersect nother segment (x3, y3, x4, y4)
    logical function SegIntersectSeg(x1,y1,x2,y2, x3, y3, x4, y4)
    !Arguments--------------------------------------------------
    real(8)                         :: x1,y1,x2,y2, x3, y3, x4, y4
    !Local------------------------------------------------------
    real(8)                         :: xi, yi, d
    !Begin------------------------------------------------------

    SegIntersectSeg = .true.

    d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);
    if (d == 0) SegIntersectSeg = .false.
    
    if (SegIntersectSeg) then
    
        xi = ((x3-x4)*(x1*y2-y1*x2)-(x1-x2)*(x3*y4-y3*x4))/d
        
        if      (abs(x1-x2)>0) then
            yi = y2 + (y1-y2)/(x1-x2) * (xi-x2)
        else if (abs(x3-x4)>0) then
            yi = y4 + (y3-y4)/(x3-x4) * (xi-x4)
        else
            SegIntersectSeg = .false.
        endif
        
        if (SegIntersectSeg) then
        
            if (abs(x1-x2) > abs(y1-y2)) then
                if (xi < min(x1,x2) .or. xi > max(x1,x2)) SegIntersectSeg = .false.
            else
                if (yi < min(y1,y2) .or. yi > max(y1,y2)) SegIntersectSeg = .false.
            endif 

            if (abs(x3-x4) > abs(y3-y4)) then
                if (xi < min(x3,x4) .or. xi > max(x3,x4)) SegIntersectSeg = .false.        
            else
                if (yi < min(y3,y4) .or. yi > max(y3,y4)) SegIntersectSeg = .false.
            endif
        endif
    endif
    
    end function SegIntersectSeg    



end module ModuleDrawing

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
