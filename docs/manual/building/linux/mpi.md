# MPI

Is an acronym standing for **Message Passing Interface**. MPI are programing directives which allow a program to launch child processes in different *nodes* (other computers) through the network and communicate between them. This allows [parallel processing](../../development/programming/parallel-processing.md). To Build a MOHID project for MPI refer to [Compiling Mohid with MPI](../../development/build/compiling-mohid-with-mpi.md).

## MPICH

**MPICH** are the libraries containing the **MPI** directives.

### Linux installation

- To build the **MPICH** libraries in linux, get them and:

```
>./configure FC=ifort --enable-f90
>make
>make testing
>make install
```

- To compile Mohid in linux in MPI refer to [this wiki](../../development/build/compiling-mohid-with-mpi.md).
- To compile any program in linux simply type

```
>mpif90 -i-static foo.f90
```

- To use an MPI program on a Linux machine simply type these lines:

```
mpd &
mpiexec -n 4 ./MohidWater
mpdallexit
```

The three lines do the following: i) install the mpi daemon, ii) run MohidWater in 4 processes, iii)once the program is finished, kill all daemons.

## Samples

### cpi

This little C program calculates pi.

```
#include "mpi.h"
#include <stdio.h>
#include <math.h>

double f(double);

double f(double a)
{
   return (4.0 / (1.0 + a*a));
}

int main(int argc,char *argv[])
{
   int    n, myid, numprocs, i;
   double PI25DT = 3.141592653589793238462643;
   double mypi, pi, h, sum, x;
   double startwtime = 0.0, endwtime;
   int    namelen;
   char   processor_name[MPI_MAX_PROCESSOR_NAME];

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   MPI_Get_processor_name(processor_name,&namelen);

   fprintf(stdout,"Process %d of %d is on %s\n",
           myid, numprocs, processor_name);
   fflush(stdout);

   n = 10000;                  /* default # of rectangles */
   if (myid == 0)
       startwtime = MPI_Wtime();
   MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

   h   = 1.0 / (double) n;
   sum = 0.0;
   /* A slightly better approach starts from large i and works back */
   for (i = myid + 1; i <= n; i += numprocs)
   {
       x = h * ((double)i - 0.5);
       sum += f(x);
   }
   mypi = h * sum;

   MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   if (myid == 0) {
       endwtime = MPI_Wtime();
       printf("pi is approximately %.16f, Error is %.16f\n",
              pi, fabs(pi - PI25DT));
       printf("wall clock time = %f\n", endwtime-startwtime);
       fflush(stdout);
   }

   MPI_Finalize();
   return 0;
}
```

### f90pi

The same little program in fortran 90:

```
!**********************************************************************
!   pi3f90.f - compute pi by integrating f(x) = 4/(1 + x**2)
!
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in top-level directory.
!
!   Each node:
!    1) receives the number of rectangles used in the approximation.
!    2) calculates the areas of it's rectangles.
!    3) Synchronizes for a global summation.
!   Node 0 prints the result.
!
!  Variables:
!
!    pi  the calculated result
!    n   number of points of integration.
!    x           midpoint of each rectangle's interval
!    f           function to integrate
!    sum,pi      area of rectangles
!    tmp         temporary scratch space for global summation
!    i           do loop index
!****************************************************************************
program main

use mpi

double precision  PI25DT
parameter        (PI25DT = 3.141592653589793238462643d0)

double precision  mypi, pi, h, sum, x, f, a
integer n, myid, numprocs, i, rc
!                                 function to integrate
f(a) = 4.d0 / (1.d0 + a*a)

call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
print *, 'Process ', myid, ' of ', numprocs, ' is alive'

sizetype   = 1
sumtype    = 2

do
   if ( myid .eq. 0 ) then
      write(6,98)
98    format('Enter the number of intervals: (0 quits)')
      read(5,99) n
99    format(i10)
   endif

   call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!                                 check for quit signal
   if ( n .le. 0 ) exit

!                                 calculate the interval size
   h = 1.0d0/n

   sum  = 0.0d0
   do i = myid+1, n, numprocs
      x = h * (dble(i) - 0.5d0)
      sum = sum + f(x)
   enddo
   mypi = h * sum

!                                 collect all the partial sums
   call MPI_REDUCE(mypi,pi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                   MPI_COMM_WORLD,ierr)

!                                 node 0 prints the answer.
   if (myid .eq. 0) then
       write(6, 97) pi, abs(pi - PI25DT)
97     format('  pi is approximately: ', F18.16, &
              '  Error is: ', F18.16)
   endif

enddo
call MPI_FINALIZE(rc)
stop
end
```

## Other references

- [Parallel processing](../../development/programming/parallel-processing.md).
- [Setup a MOHID simulation using MPI](../../development/setup-a-mohid-simulation-using-mpi.md).
- [Compiling Mohid with MPI](../../development/build/compiling-mohid-with-mpi.md).
- [MPIConfig.dat](../../water/mpiconfig-dat.md).
