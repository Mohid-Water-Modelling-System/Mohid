!  FromASCIIGMTtoBNL.f90 
!
!  FUNCTIONS:
!	FromASCIIGMTtoBNL      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: FromASCIIGMTtoBNL
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

	program FromASCIIGMTtoBNL

	implicit none

    type T_Line
        real,  dimension(:), pointer :: x,y
        integer                      :: Number
        logical                      :: Deleted          = .false. 
        logical                      :: Order            = .false.
    end type T_Line

	! Variables
    character(len = 100)             :: CharC, FileIn, FileOut
    real,   dimension (0:200000)       :: StartBlock, EndBlock
    real,   dimension (3000000)       :: x,y
    integer                          :: i, j, ntotal, n, jblocks, OutputType
    
    type(T_Line), dimension (:), pointer :: Line, LineAux
    real ,        dimension (:), pointer :: xaux, yaux
    real                                 :: dxy1, dxy2, dxy3, dxy4, dx, dy
    integer                              :: il, l, k, MaxN, imax
    real                                 :: SearchRadius = 0.1
    integer                              :: MinPoints = 50, MinWritePoints = 200
    logical                              :: MergeLines = .true., LinesFound
    integer                              :: StartBlockON, EndBlockON 



	! Body of FromASCIIGMTtoBNL

    
    write(*,*) 'Qual o ficheiro de input ?'
    
    read (*,*) FileIn
    !FileIn ="AtlanticNE.gmt"

    write(*,*) 'Qual o ficheiro de output ?'
    read (*,*) FileOut
    !FileOut="AtlanticNE.lin"


    open(1,file=trim(FileIn))
    open(2,file=trim(FileOut))

    write(*,*) 'Qual o formato de saída?'
    write(*,*) 'Surfer   - 1'
    write(*,*) 'MohidGis - 2'

    read (*,*) OutputType 
    !OutputType = 2
    if (OutputType /= 1 .and. OutputType /= 2) stop 'Wrong output type'
    
    i = 1
    j = 1

    read(1,*)

    do 
        read(1,'(A100)',end=10) CharC

        if (CharC(1:1) =='>') then 
            StartBlock(j  ) = i
            EndBlock  (j-1) = i - 1
            j = j + 1
        else
            read(CharC,*) x(i), y(i)
            i = i+1
        endif

    enddo
10  continue

	ntotal = i-2
    jblocks = j-1
    EndBlock  (jblocks) = ntotal

    allocate(Line(1:jblocks))


    do j = 1, jblocks

        if (OutputType == 2) then
            EndBlock(j) = EndBlock(j) - 1
        endif

        Line(j)%Number = EndBlock(j) - StartBlock(j) + 1

        allocate(Line(j)%x(1:Line(j)%Number),Line(j)%y(1:Line(j)%Number))

        n = Line(j)%Number

        Line(j)%x(1:n) = x(StartBlock(j): EndBlock(j))
        Line(j)%y(1:n) = y(StartBlock(j): EndBlock(j))

    enddo

    if (MergeLines) then

        ! merge lines
        do j = 1, jblocks

            n = Line(j)%Number
        
            if (n > MinPoints .and. .not. Line(j)%Deleted) then

                LinesFound = .true.

                do while (LinesFound)

                do k = 1, jblocks

                    LinesFound = .false.    

                    if (k/=j .and. .not. Line(k)%Deleted) then

                        l = Line(k)%Number

                        if (l < MinPoints) cycle

                        dx = Line(j)%x(1) - Line(k)%x(1)
                        dy = Line(j)%y(1) - Line(k)%y(1)

                        dxy1 = sqrt(dx*dx + dy*dy)

                        dx = Line(j)%x(1) - Line(k)%x(l)
                        dy = Line(j)%y(1) - Line(k)%y(l)

                        dxy2 = sqrt(dx*dx + dy*dy)

                        dx = Line(j)%x(n) - Line(k)%x(1)
                        dy = Line(j)%y(n) - Line(k)%y(1)

                        dxy3 = sqrt(dx*dx + dy*dy)

                        dx = Line(j)%x(n) - Line(k)%x(l)
                        dy = Line(j)%y(n) - Line(k)%y(l)

                        dxy4 = sqrt(dx*dx + dy*dy)

                        if (min(dxy1, dxy2, dxy3, dxy4) < SearchRadius) then

                            LinesFound = .true.

                            allocate (xaux(1:n+l), yaux(1:n+l))

                            Line(k)%Deleted = .true.

                            !do not invert points order 
                            if      (dxy3 < SearchRadius) then

                                xaux(1  :n  ) = Line(j)%x(1:n); 
                                yaux(1  :n  ) = Line(j)%y(1:n);
                                xaux(n+1:l+n) = Line(k)%x(1:l); 
                                yaux(n+1:l+n) = Line(k)%y(1:l);

                            else if (dxy2 < SearchRadius) then

                                xaux(1  :l  ) = Line(k)%x(1:l);
                                yaux(1  :l  ) = Line(k)%y(1:l);
                                xaux(l+1:l+n) = Line(j)%x(1:n);
                                yaux(l+1:l+n) = Line(j)%y(1:n);

                            !invert points order
                            else if (dxy1 < SearchRadius) then

                                do il = 1, l
                                    xaux(1 :l) = Line(k)%x(l-il+1); 
                                    yaux(1 :l) = Line(k)%y(l-il+1);                            
                                enddo
                                xaux(l+1:l+n) = Line(j)%x(1:n); 
                                yaux(l+1:l+n) = Line(j)%y(1:n);

                            else if (dxy4 < SearchRadius) then

                                do il = 1, l
                                    xaux(n+1:n+l) = Line(k)%x(l-il+1); 
                                    yaux(n+1:n+l) = Line(k)%y(l-il+1);                            
                                enddo
                                xaux(1:n) = Line(j)%x(1:n);
                                yaux(1:n) = Line(j)%y(1:n);

                            endif

                            deallocate(Line(j)%x, Line(j)%y)

                            n  = l + n
                            Line(j)%Number = n

                            allocate(Line(j)%x(1:n),Line(j)%y(1:n))

                            Line(j)%x(1:n) = xaux(1:n)
                            Line(j)%y(1:n) = yaux(1:n)

                            exit

                        endif

                    endif

                enddo
                enddo
            endif

        enddo

    endif

    !order by size
    allocate(LineAux(jblocks))
    n = 0

    do i = 1, jblocks
        MaxN = -999
        do j = 1, jblocks
            if (.not. Line(j)%Deleted .and. Line(j)%Number > MinWritePoints .and. .not. Line(j)%Order) then

                if (Line(j)%Number > MaxN) then
                    MaxN = Line(j)%Number
                    iMax = j
                endif

            endif

        enddo    
        if (MaxN > -999) then
            LineAux(i)       =  Line(iMax)
            Line(iMax)%Order = .true.
            n = n + 1
        endif
    enddo



    do j = 1, n

        if (OutputType == 1) then
            write(2,*)  LineAux(j)%Number, ",1"
        endif

        if (OutputType == 2) then
            write(2,*) "!Polygon number= ",j
            write(2,*) "!Number of points= ",LineAux(j)%Number
            write(2,*) '<begin_line>'
        endif

        do i=1,LineAux(j)%Number

            write(2,*)  LineAux(j)%x(i), LineAux(j)%y(i)

        enddo
  
        if (OutputType == 2) then
            write(2,*) '<end_line>'
        endif

    enddo

    close(1)
    close(2)

	end program FromASCIIGMTtoBNL

