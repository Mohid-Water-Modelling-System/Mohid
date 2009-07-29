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


end module ModuleDrawing

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
