!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : GridGenerator
! PROGRAM       : GridGenerator
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2004
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Program to generate regular grids with non-constant spacing
!
!------------------------------------------------------------------------------
 
!Data file - default name 'GridGenerator.dat' (must be placed in working directory)

!   OUTPUT_FILE                 : char              -           !Output file name
!   ORIGIN_X                    : real              -           !Coordinate X of the lower left corner
!   ORIGIN_Y                    : real              -           !Coordinate Y of the lower left corner
!   GRID_ANGLE                  : real              -           !Rotation angle of the grid
!   COORD_TIP                   : int               -           !Coordinate type


!<begin_grid_xx>
!   NUMBER_OF_SEGMENTS          : int               -           !Number of segments in the XX direction
!   NODES                       : real vector       -           !Coordinates starting in 0. of the XX nodes
!   RESOLUTION                  : real vector       -           !Resolution at each of the defined XX nodes
!<end_grid_xx>

!<begin_grid_yy>
!   NUMBER_OF_SEGMENTS          : int               -           !Number of segments in the YY direction
!   NODES                       : real vector       -           !Coordinates starting in 0. of the YY nodes
!   RESOLUTION                  : real vector       -           !Resolution at each of the defined YY nodes
!<end_grid_yy>

program GridGenerator

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions

    implicit none

    !Time variables------------------------------------------------------------
    type (T_Time)                       :: InitialSystemTime, FinalSystemTime
    real                                :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)               :: F95Time

    !Input data file-----------------------------------------------------------
    character(len = PathLength)         :: DataFile = 'GridGenerator.dat'

    !Parameters----------------------------------------------------------------
    integer, parameter                  :: XX = 1
    integer, parameter                  :: YY = 2

    !Types---------------------------------------------------------------------
    type     T_Resolution
        real                            :: Initial
        real                            :: Final
    end type T_Resolution
        
    type     T_Segment
        real                            :: StartAt
        real                            :: EndAt
    end type T_Segment

    type     T_GridSegment
        integer                         :: Direction
        type(T_Resolution )             :: Resolution
        type(T_Segment    )             :: Segment
        type(T_GridSegment),    pointer :: Next
    end type T_GridSegment

    type     T_GridGenerator
        integer                         :: ObjEnterData = 0
        integer                         :: OutputUnit
        integer                         :: CoordType
        real                            :: GridAngle
        real                            :: OriginX
        real                            :: OriginY
        character(PathLength)           :: OutputFileName
        type(T_Size2D       )           :: Size
        type(T_GridSegment  ),  pointer :: FirstGridSegment
    end type T_GridGenerator

    type(T_GridGenerator),      pointer :: Me

    !Begin---------------------------------------------------------------------

    nullify(Me); allocate(Me); nullify(Me%FirstGridSegment)

    call ConstructGridGenerator
    call GenerateGrid
    call KillGridGenerator

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructGridGenerator
        
        call StartUpMohid("GridGenerator")

        call StartCPUTime

        call ReadOptions

    end subroutine ConstructGridGenerator
    
    !--------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        call ConstructEnterData (Me%ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - GridGenerator - ERR01'

        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'OUTPUT_FILE',                      &
                     ClientModule = 'GridGenerator',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - GridGenerator - ERR02'

        call GetData(Me%OriginX,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'ORIGIN_X',                         &
                     ClientModule = 'GridGenerator',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - GridGenerator - ERR03'

        call GetData(Me%OriginY,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'ORIGIN_Y',                         &
                     ClientModule = 'GridGenerator',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - GridGenerator - ERR04'

        call GetData(Me%GridAngle,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'GRID_ANGLE',                       &
                     ClientModule = 'GridGenerator',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - GridGenerator - ERR05'

        call GetData(Me%CoordType,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'COORD_TIP',                        &
                     ClientModule = 'GridGenerator',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - GridGenerator - ERR05'

        call ConstructXX

        call ConstructYY

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - GridGenerator - ERR05'

    end subroutine ReadOptions

    
    !--------------------------------------------------------------------------

    
    subroutine ConstructXX
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        type(T_GridSegment), pointer                :: NewGridSegment
        real, dimension(:),  pointer                :: NodesX, ResolutionX
        integer                                     :: nSegments, nNodes, iSegment

        !Begin-----------------------------------------------------------------
        
        BlockFound = OFF

        call ExtractBlockFromBuffer(Me%ObjEnterData,                            &
                                    ClientNumber,                               &
                                    '<begin_grid_xx>',                          &
                                    '<end_grid_xx>',                            &
                                    BlockFound,                                 &
                                    STAT = STAT_CALL)

        if(STAT_CALL .EQ. SUCCESS_ )then    

            if (BlockFound) then

                call GetData(nSegments,                                         &
                             Me%ObjEnterData, iflag,                            &
                             SearchType   = FromBlock,                          &
                             keyword      = 'NUMBER_OF_SEGMENTS',               &
                             ClientModule = 'GridGenerator',                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructXX - GridGenerator - ERR01'

                nNodes = nSegments + 1

                nullify (NodesX,           ResolutionX          ) 
                allocate(NodesX(1:nNodes), ResolutionX(1:nNodes))

                call GetData(NodesX,                                            &
                             Me%ObjEnterData, iflag,                            &
                             SearchType   = FromBlock,                          &
                             keyword      = 'NODES',                            &
                             ClientModule = 'GridGenerator',                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructXX - GridGenerator - ERR02'

                call GetData(ResolutionX,                                       &
                             Me%ObjEnterData, iflag,                            &
                             SearchType   = FromBlock,                          &
                             keyword      = 'RESOLUTION',                       &
                             ClientModule = 'GridGenerator',                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructXX - GridGenerator - ERR03'


                do iSegment = 1, nSegments
                    
                    call AddGridSegment(NewGridSegment)

                    NewGridSegment%Direction            = XX
                    NewGridSegment%Segment%StartAt      = NodesX     (iSegment    )
                    NewGridSegment%Segment%EndAt        = NodesX     (iSegment + 1)
                    NewGridSegment%Resolution%Initial   = ResolutionX(iSegment    )
                    NewGridSegment%Resolution%Final     = ResolutionX(iSegment + 1)

                end do

                deallocate(NodesX, ResolutionX); nullify(NodesX, ResolutionX)
               
            endif

        elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then 

            stop 'ConstructXX - GridGenerator - ERR04'

        end if

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructXX - GridGenerator - ERR05'

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructXX - GridGenerator - ERR06'


    end subroutine ConstructXX

    !--------------------------------------------------------------------------

    subroutine ConstructYY

        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        type(T_GridSegment), pointer                :: NewGridSegment
        real, dimension(:),  pointer                :: NodesY, ResolutionY
        integer                                     :: nSegments, nNodes, iSegment

        !Begin-----------------------------------------------------------------
        
        BlockFound = OFF

        call ExtractBlockFromBuffer(Me%ObjEnterData,                            &
                                    ClientNumber,                               &
                                    '<begin_grid_yy>',                          &
                                    '<end_grid_yy>',                            &
                                    BlockFound,                                 &
                                    STAT = STAT_CALL)

        if(STAT_CALL .EQ. SUCCESS_ )then    

            if (BlockFound) then
                
                call GetData(nSegments,                                         &
                             Me%ObjEnterData, iflag,                            &
                             SearchType   = FromBlock,                          &
                             keyword      = 'NUMBER_OF_SEGMENTS',               &
                             ClientModule = 'GridGenerator',                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructYY - GridGenerator - ERR01'

                nNodes = nSegments + 1

                nullify (NodesY,           ResolutionY          ) 
                allocate(NodesY(1:nNodes), ResolutionY(1:nNodes))

                call GetData(NodesY,                                            &
                             Me%ObjEnterData, iflag,                            &
                             SearchType   = FromBlock,                          &
                             keyword      = 'NODES',                            &
                             ClientModule = 'GridGenerator',                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructYY - GridGenerator - ERR02'

                call GetData(ResolutionY,                                       &
                             Me%ObjEnterData, iflag,                            &
                             SearchType   = FromBlock,                          &
                             keyword      = 'RESOLUTION',                       &
                             ClientModule = 'GridGenerator',                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructYY - GridGenerator - ERR03'

                do iSegment = 1, nSegments
                    
                    call AddGridSegment(NewGridSegment)

                    NewGridSegment%Direction            = YY
                    NewGridSegment%Segment%StartAt      = NodesY     (iSegment    )
                    NewGridSegment%Segment%EndAt        = NodesY     (iSegment + 1)
                    NewGridSegment%Resolution%Initial   = ResolutionY(iSegment    )
                    NewGridSegment%Resolution%Final     = ResolutionY(iSegment + 1)

                end do

                deallocate(NodesY, ResolutionY); nullify(NodesY, ResolutionY)

            endif

        elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then 

            stop 'ConstructYY - GridGenerator - ERR04'

        end if

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructYY - GridGenerator - ERR05'


    end subroutine ConstructYY 


    !--------------------------------------------------------------------------
    
    subroutine AddGridSegment (ObjGridSegment)

        !Arguments-------------------------------------------------------------
        type (T_GridSegment), pointer                   :: ObjGridSegment

        !Local-----------------------------------------------------------------
        type (T_GridSegment), pointer                   :: NewGridSegment
        type (T_GridSegment), pointer                   :: PreviousGridSegment
        
        !Allocates new instance
        allocate (NewGridSegment)
        nullify  (NewGridSegment%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Me%FirstGridSegment)) then
            Me%FirstGridSegment      => NewGridSegment
            ObjGridSegment           => NewGridSegment
        else
            PreviousGridSegment      => Me%FirstGridSegment
            ObjGridSegment           => Me%FirstGridSegment%Next
            do while (associated(ObjGridSegment))
                PreviousGridSegment  => ObjGridSegment
                ObjGridSegment       => ObjGridSegment%Next
            enddo
            ObjGridSegment           => NewGridSegment
            PreviousGridSegment%Next => NewGridSegment
        endif


    end subroutine AddGridSegment
    
    
    !--------------------------------------------------------------------------


    subroutine GenerateGrid
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------


        call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GenerateGrid - GridGenerator - ERR01'
        
        open(unit = Me%OutputUnit, file = Me%OutputFileName, status = 'unknown')

        call ComputeXX

        call ComputeYY

        write(Me%OutputUnit, *)'ILB_IUB         :', Me%Size%ILB, Me%Size%IUB
        write(Me%OutputUnit, *)'JLB_JUB         :', Me%Size%JLB, Me%Size%JUB
        write(Me%OutputUnit, *)'ORIGIN          :', Me%OriginX, Me%OriginY
        write(Me%OutputUnit, *)'GRID_ANGLE      :', Me%GridAngle
        write(Me%OutputUnit, *)'COORD_TIP       :', Me%CoordType

        call UnitsManager(Me%OutputUnit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GenerateGrid - GridGenerator - ERR01'

    
    end subroutine GenerateGrid

    !--------------------------------------------------------------------------

    subroutine ComputeXX
        
        !Local-----------------------------------------------------------------
        type (T_GridSegment), pointer                   :: GridSegment

        !Begin-----------------------------------------------------------------

        nullify(GridSegment)

        Me%Size%JLB = 1
        Me%Size%JUB = 0

        write(Me%OutputUnit, *)'<BeginXX>'
        write(Me%OutputUnit, *) 0.

        GridSegment => Me%FirstGridSegment

        do while(associated(GridSegment))

            if    (GridSegment%Direction .eq. XX)then
                
                call OneSegment(GridSegment%Resolution%Initial, &
                                GridSegment%Resolution%Final,   &
                                GridSegment%Segment%StartAt,    &
                                GridSegment%Segment%EndAt,      &
                                Me%Size%JUB)


            end if

            GridSegment => GridSegment%Next
        end do

        write(Me%OutputUnit, *)'<EndXX>'
        
        nullify(GridSegment)

    end subroutine ComputeXX
    
    !--------------------------------------------------------------------------

    subroutine ComputeYY
        
        !Local-----------------------------------------------------------------
        type (T_GridSegment), pointer                   :: GridSegment

        !Begin-----------------------------------------------------------------

        nullify(GridSegment)

        Me%Size%ILB = 1
        Me%Size%IUB = 0

        write(Me%OutputUnit, *)'<BeginYY>'
        write(Me%OutputUnit, *) 0.

        GridSegment => Me%FirstGridSegment

        do while(associated(GridSegment))

            if    (GridSegment%Direction .eq. YY)then
                
                call OneSegment(GridSegment%Resolution%Initial, &
                                GridSegment%Resolution%Final,   &
                                GridSegment%Segment%StartAt,    &
                                GridSegment%Segment%EndAt,      &
                                Me%Size%IUB)


            end if

            GridSegment => GridSegment%Next
        end do
        
        write(Me%OutputUnit, *)'<EndYY>'

        nullify(GridSegment)


    end subroutine ComputeYY
    

    !--------------------------------------------------------------------------

    subroutine KillGridGenerator

        call StopCPUTime

        call ShutdownMohid ("GridGenerator", ElapsedSeconds, TotalCPUTime)

    end subroutine KillGridGenerator
    
    !--------------------------------------------------------------------------

    subroutine StartCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)

    end subroutine StartCPUTime
    
    !--------------------------------------------------------------------------

    subroutine StopCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        
        call cpu_time(TotalCPUTime)

        ElapsedSeconds = FinalSystemTime - InitialSystemTime

    end subroutine StopCPUTime
    
    !--------------------------------------------------------------------------
    
    subroutine OneSegment(StartResolution, EndResolution, StartAt, EndAt, count)

        !Arguments--------------------------------------------------
        real                        :: StartResolution        
        real                        :: EndResolution        
        real,    intent(in)         :: StartAt, EndAt
        integer                     :: count, i

        !Local------------------------------------------------------
        real, dimension(:), pointer :: segment
        integer                     :: n
        
        !Begin------------------------------------------------------

        call gera(StartResolution, EndResolution, StartAt, EndAt, segment)

        n = size(segment) - 1 

        count = count + n

        do i = 1, n
            write(Me%OutputUnit,*) segment(i)
        end do

        deallocate(segment); nullify(segment) 

    end subroutine OneSegment

    !--------------------------------------------------------------------------

    subroutine gera(StartRes, EndRes, StartAt, EndAt, segment)
        
        !Arguments--------------------------------------------------
        real                        :: StartRes, EndRes
        real                        :: StartAt, EndAt
        real, dimension(:), pointer :: segment

        !Local------------------------------------------------------
        real                        :: r                = 0.
        real                        :: x                = 0.
        real                        :: error            = 0.
        real                        :: SumDistances     = 0.
        integer                     :: i                = 0
        integer                     :: n
        real                        :: TotalLength
        real, dimension(:), pointer :: vector
        logical                     :: Swap, EqualDistance
        real                        :: OriginalStartRes, OriginalEndRes
        real                        :: StartResolution, EndResolution       
        real                        :: dx

        !Begin------------------------------------------------------

        nullify(segment)

        OriginalStartRes = StartRes
        OriginalEndRes = EndRes

        if    (StartRes > EndRes)   then

            Swap            = .true.
            EqualDistance   = .false.
            StartResolution = OriginalEndRes
            EndResolution   = OriginalStartRes

        elseif(StartRes .eq. EndRes)then

            EqualDistance   = .true.
            StartResolution = StartRes
            EndResolution   = EndRes

        else

            StartResolution = StartRes
            EndResolution   = EndRes
            Swap            = .false.
            EqualDistance   = .false.

        end if


        if(.not. EqualDistance)then

            n = 0; TotalLength = EndAt - StartAt

            r = (TotalLength - EndResolution)/(TotalLength - StartResolution)

            SumDistances = 0.

            n = n + 1; x = EndResolution

            SumDistances = SumDistances + x

            do while(x > StartResolution)

                x = r * x; SumDistances = SumDistances + x; n = n + 1

            end do

            error =  SumDistances - TotalLength

            allocate(vector (0:n))
            allocate(segment(0:n))

            error = error/float(n)

            vector(0) = 0.;
            vector(1) = EndResolution 

            do i = 2, n

                vector(i) = vector(i-1) * r

            end do

            do i = 1, n

                vector (i) = vector(i) + vector(i-1) - error

            end do

            segment(0) = StartAt

            if(.not. Swap)then

                do i = n, 0, -1                

                    segment(i) = EndAt - vector(n-i) 

                end do

            else

                do i = 1, n

                    segment(i) = StartAt + vector(i) 

                end do

            end if

        else

            !Equal distance

            TotalLength     = EndAt - StartAt
            dx              = StartResolution
            n               = int(TotalLength/dx)

            allocate(vector (0:n))
            allocate(segment(0:n))

            dx     = TotalLength/n

            vector(0) = 0.

            do i = 1, n

                vector (i) = vector(i-1) + dx

            end do

            segment(0) = StartAt

            do i = 1, n

                segment(i) = StartAt + vector(i) 

            end do

        end if 


        deallocate(vector)
        nullify   (vector)

    end subroutine gera

end program GridGenerator
