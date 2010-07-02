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

	! Variables
    character(len = 100)    :: CharC, FileIn, FileOut
    real, dimension (0:2000) :: StartBlock, EndBlock
    real, dimension (300000) :: x,y
    integer                 :: i, j, ntotal, n, jblocks, OutputType


	! Body of FromASCIIGMTtoBNL

    
    write(*,*) 'Qual o ficheiro de input ?'
    
    read (*,*) FileIn
    !FileIn ="MedWest.gmt"

    write(*,*) 'Qual o ficheiro de output ?'
    read (*,*) FileOut
    !FileOut="MedWest.lin"


    open(1,file=FileIn)
    open(2,file=FileOut)

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

    do j = 1, jblocks

        n = EndBlock(j) - StartBlock(j) + 1

        if (OutputType == 1) then
            write(2,*) n, ",1"
        endif

        if (OutputType == 2) then
            write(2,*) '<begin_line>'
            EndBlock(j) = EndBlock(j) - 1
        endif


        do i=StartBlock(j), EndBlock(j)

            write(2,*) x(i)-360,y(i)

        enddo

        if (OutputType == 2) then
            write(2,*) '<end_line>'
        endif


    enddo

    close(1)
    close(2)

	end program FromASCIIGMTtoBNL

