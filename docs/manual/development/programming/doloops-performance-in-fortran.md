# DoLoops performance in Fortran

What is the best performance that Fortran can give when computing do-loops over large matrices? A single triple-do-loop test-case was implemented with matrix size ranging from 1M to 216M four-byte units. The test-case below shows:

- a 400% performance increase when looping with the order (k,j,i) instead of looping with the order (i,j,k),
- that a 300% performance increase occur when using openmp directives in a quad-core processor (i7-870),
- that changing the chunk size or alternating between static and dynamic scheduling yield less than 10% differences in performance. Best performance is obtained, for this test-case, with a small dynamic chunk.

## Limiting number of threads before running

If one wishes to limit the maximum number of threads to 2 when running openmp-able executables:

```
> set OMP_NUM_THREADS=2
> MohidWater_omp.exe
```

## Test-cases

### Simple triple do-loop

#### Description

Performs a triple-do-loop with a simple computation over a cubic matrix. Size is chosen by the user through standard input.

#### Hardware

- Intel Core i7 - 870
- 8 GB Ram

#### Code

- The main program

```
program DoloopsOpenmp

   use moduleDoloopsOpenmp, only: makeloop

   implicit none

   integer, dimension(:,:,:), pointer      :: mycube
   integer                                 :: M = 1
   real                                    :: elapsedtime
   real                                    :: time = 0.0

   do while (M < 1000)

       write(*,*) 'Insert the cube size M (or insert 1000 to exit): '
       read(*,*) M

       if (M > 999) exit

       allocate(mycube(1:M,1:M,1:M))

       !Tic()
       time = elapsedtime(time)

       call makeloop(mycube)

       !Toc()
       time = elapsedtime(time)

       write(*,10) time
       write(*,*)

       deallocate(mycube)
       nullify(mycube)

   end do

   10 format('Time elapsed: ',F6.2)

end program DoloopsOpenmp

!This function computes the time
real function elapsedtime(lasttime)

   integer     :: count, count_rate
   real        :: lasttime

   call system_clock(count, count_rate)

   elapsedtime = count * 1.0 / count_rate - lasttime

end function elapsedtime
```

---

- The module

```
module moduleDoloopsOpenmp

   use omp_lib

   implicit none

   private

   public :: makeloop

contains

   subroutine makeloop(cubicmatrix)

       !Arguments --------------
       integer, dimension(:,:,:), pointer  :: cubicmatrix

       !Local variables --------
       integer                             :: i, j, k, lb, ub

       lb = lbound(cubicmatrix,3)
       ub = ubound(cubicmatrix,3)

       !$OMP PARALLEL PRIVATE(i,j,k)
       !$OMP DO
       do k = lb, ub
       do j = lb, ub
       do i = lb, ub

           cubicmatrix(i,j,k) = cubicmatrix(i,j,k) + 1

       end do
       end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL

   end subroutine makeloop

end module moduleDoloopsOpenmp
```

#### Results

http://content.screencast.com/users/GRiflet/folders/Jing/media/49c2beae-32f4-4da0-96b1-ddae90bd2d02/2010-10-21\_1735.png

- Full results

http://content.screencast.com/users/GRiflet/folders/Jing/media/9e999f22-7204-481a-84be-495ffe35e9f7/2010-10-21\_1734.png

- Looking only at the results with STATIC/DYNAMIC/CHUNK variations.

```
----------------------------
DO (i,j,k) / NO CHUNK
----------------------------
Table A.1 - Debug do(i,j,k)
Size	Time
100		 0.04
200		 0.37
300		 1.58
400		 7.60
500		19.66
600		41.65

Table A.2 - Debug openmp without !$OMP PARALLEL directives do(i,j,k)
Size	Time
100		 0.04
200		 0.37
300		 1.58
400		 7.27
500		19.34
600		41.34

Table A.3 - Debug openmp with one !$OMP PARALLEL DO directive do(i,j,k)
Size	Time
100		 0.02
200		 0.19
300		 0.70
400		 1.86
500		 4.05
600		 7.83

----------------------------
DO (k,j,i) / NO CHUNK
----------------------------

Table B.1 - Debug do(k,j,i)
Size	Time
100		 0.04
200		 0.31
300		 1.22
400		 3.36
500		 7.55
600		14.88

Table B.2 - Debug openmp without !$OMP PARALLEL directives do(k,j,i)
Size	Time
100		 0.04
200		 0.31
300		 1.21
400		 3.36
500		 7.82
600		15.07

Table B.3 - Debug openmp with one !$OMP PARALLEL DO directive do(k,j,i)
Size	Time
100		 0.02
200		 0.09
300		 0.36
400		 0.94
500		 2.04
600		 3.89

----------------------------
DO (k,j,i) / STATIC CHUNK = (UBOUND - LBOUND) / NTHREADS + 1
----------------------------

Table C.3 - Debug openmp with one !$OMP PARALLEL DO directive do(k,j,i)
Size	Time
100		 0.02
200		 0.15
300		 0.42
400		 1.03
500		 2.12
600		 3.97

----------------------------
DO (k,j,i) / STATIC CHUNK = 10
----------------------------

Table D.3 - Debug openmp with one !$OMP PARALLEL DO directive do(k,j,i)
Size	Time
100		 0.02
200		 0.16
300		 0.43
400		 1.04
500		 2.18
600		 4.05

----------------------------
DO (k,j,i) / DYNAMIC CHUNK = 10
----------------------------

Table E.3 - Debug openmp with one !$OMP PARALLEL DO directive do(k,j,i)
Size	Time
100		 0.01
200		 0.10
300		 0.36
400		 0.93
500		 2.01
600		 3.89

----------------------------
DO (k,j,i) / DYNAMIC CHUNK = (UBOUND - LBOUND) / NTHREADS + 1
----------------------------

Table F.3 - Debug openmp with one !$OMP PARALLEL DO directive do(k,j,i)
Size	Time
100		 0.02
200		 0.09
300		 0.39
400		 1.04
500		 2.14
600		 4.00
```

#### Conclusions

- do(k,j,i) Vs do(i,j,k) ==> 2 to 4 times faster!
- dynamic small chunks, or no chunk at all yield 10% increased performance over large dynamic chunks. Probably better off with no-chunk.
- More test-cases representing different scenarios of do-loops may yield different choices of CHUNK/scheduling.
- Single precision computation over large numbers (such as summing the entries in a large matrix) yield significant errors. Furthermore, the results yielded are different between openmp and no-openmp.
- Double precision computation yields the correct results. The results are the same between openmp and no-openmp.

---

### SetMatrix3D\_Constant

#### Description

This subroutine is in the ModuleFunctions of MohidWater. In the context of MohidWater, parallelizing this subroutine yields up to 15% increase in performance. However, in the little test program, the same OMP directives yield quite good results (under 1/3 of the simulation time or a 200% increase in performance).

#### Hardware

- Core i7-870
- 8 GB RAM

#### Code

```
   real function SetMatrixValues3D_R8_Constant (Matrix, Valueb, MapMatrix)

       !Arguments-------------------------------------------------------------
       real, dimension(:, :, :), pointer               :: Matrix
       real, intent (IN)                               :: Valueb
       integer, dimension(:, :, :), pointer, optional  :: MapMatrix

       !Local-----------------------------------------------------------------
       integer                                         :: i, j, k
       integer                                         :: ilb, iub, jlb, jub, klb, kub

       !Begin-----------------------------------------------------------------

       ilb = lbound(Matrix,1)
       iub = ubound(Matrix,1)
       jlb = lbound(Matrix,2)
       jub = ubound(Matrix,2)
       klb = lbound(Matrix,3)
       kub = ubound(Matrix,3)

       !griflet: omp slowdown
       if (present(MapMatrix)) then
           !$OMP PARALLEL DO PRIVATE(i,j,k)
           do k = klb, kub
           do j = jlb, jub
           do i = ilb, iub
               if (MapMatrix(i, j, k) == 1) then
                   Matrix (i, j, k) = Valueb
               endif
           enddo
           enddo
           enddo
           !$OMP END PARALLEL DO
       else
           !$OMP PARALLEL DO PRIVATE(i,j,k)
           do k = klb, kub
           do j = jlb, jub
           do i = ilb, iub
               Matrix (i, j, k) = Valueb
           enddo
           enddo
           enddo
           !$OMP END PARALLEL DO
       endif

       SetMatrixValues3D_R8_Constant = sumMatrix3D(Matrix)

   end function SetMatrixValues3D_R8_Constant
```

#### Results

http://content.screencast.com/users/GRiflet/folders/Jing/media/136f62db-a24f-4315-9eae-f4bcd9110d4e/2010-10-22\_1553.png

## MOHID

The MOHID parallelization is a complex matter because:

- Time keeping is hard to keep with due to the fact that CPU time is the sum of the computation time of each thread.
- Do loops that parallelize very well in small programs add a lot of overhead in big programs like MOHID and actually tend to decrease performance.

This wiki-entry comments which loops are efficiently parallelized in each module.

### Hardware

- Intel Core i7 - 870
- 8 GB Ram

### Compiler options

- Here are the different compiler options used throughout this test-case, as seen from visual studio 2008:

http://i.imgur.com/KcDFt.png http://i.imgur.com/hNxrI.png http://i.imgur.com/638x0.png http://i.imgur.com/xEPnq.png http://i.imgur.com/utpq1.png http://i.imgur.com/5F6e4.png http://i.imgur.com/qglDB.png http://i.imgur.com/fJZBQ.png http://i.imgur.com/nYcMC.png http://i.imgur.com/qVw0M.png http://i.imgur.com/KpTcx.png http://i.imgur.com/5bzy6.png http://i.imgur.com/1bhkA.png http://i.imgur.com/V3U3d.png http://i.imgur.com/0bDqn.png

### PCOMS test-case

- Here's the present situation with the codeplex build (20101029). The chart below depicts the PCOMS test-case performance with the growing number of threads (up to 8). Maximum performance gains are roughly 15% for the 4 threads. Since the i7-870 is a 4 core machine, it makes sense that 4 threads perform better than 5 or more, or than 3 or less. Also note that a single-threaded openmp code is slower by 9% than a no-openmp code.

http://i.imgur.com/dhSg5.png

#### Without openmp compiler option

A 3 hour run of the PCOMS is made, which takes around **870s** without parallelization.

Here's an excerpt of the outwatch log.

```
                         Main ModifyMohidWater                          863.33
              ModuleFunctions THOMASZ                                    46.06
              ModuleFunctions SetMatrixValues3D_R8_Constant              40.09
              ModuleFunctions SetMatrixValues3D_R8_FromMatri             16.69
              ModuleFunctions InterpolateLinearyMatrix3D                  8.42
```

Another 3 hour run with the above compiler settings takes, rougly, **425s** without parallelization: http://i.imgur.com/R3ZEc.png

#### With openmp compiler option, with current code from Codeplex (codename: Angela)

A 1 hour run of the PCOMS is made, and takes around **400s** with parallelization with 8 threads.

- All threads (8):

http://content.screencast.com/users/GRiflet/folders/Jing/media/c8b434c0-0ab0-4931-9491-4e018c3ce566/2010-10-27\_1705.png

Here's an excerpt of the outwatch log:

```
                         Main ModifyMohidWater                          346.76
              ModuleFunctions SetMatrixValues3D_R8_Constant               5.02
              ModuleFunctions SetMatrixValues3D_R8_FromMatri              1.83
              ModuleFunctions InterpolateLinearyMatrix3D                  1.10
              ModuleFunctions SetMatrixValues2D_R8_Constant               0.13
              ModuleFunctions InterpolateLinearyMatrix2D                  0.08
```

- 1 Thread only (set OMP\_NUM\_THREADS=1) takes more than **460s**:

http://i.imgur.com/zY3K3.png

- 2 Threads only (set OMP\_NUM\_THREADS=2) take less than **380s**:

http://i.imgur.com/CpFvT.png

- 3 Threads only (set OMP\_NUM\_THREADS=3) take less than **370s**:

http://i.imgur.com/TKBrC.png

- 4 Threads only (set OMP\_NUM\_THREADS=3) take less than **360s**:

http://i.imgur.com/VNFuQ.png~

- 5 Threads only (set OMP\_NUM\_THREADS=5) take more than **370s**:

http://i.imgur.com/6dOxk.png

- 6 Threads only (set OMP\_NUM\_THREADS=6) take more than **375s**:

http://i.imgur.com/n55Vo.png

- 7 Threads only (set OMP\_NUM\_THREADS=7) take more than **385s**:

http://i.imgur.com/693cw.png

- 8 Threads only (set OMP\_NUM\_THREADS=8) take more than **390S**:

http://i.imgur.com/3yjtc.png

#### With openmp, but without any openmp directives

Here's an except of the outwatch log:

```
                         Main ModifyMohidWater                          936.25
              ModuleFunctions THOMASZ                                    46.38
              ModuleFunctions SetMatrixValues3D_R8_Constant              40.12
              ModuleFunctions SetMatrixValues3D_R8_FromMatri             16.63
              ModuleFunctions InterpolateLinearyMatrix3D                  8.24
              ModuleFunctions THOMAS_3D                                   0.57
```

### ModuleFunctions

- After parallelizing the Module Functions only.

http://i.imgur.com/BlHbM.png

Here's an excerpt of the outwatch log:

```
                         Main ModifyMohidWater                          945.70
              ModuleFunctions THOMASZ                                    46.89
              ModuleFunctions SetMatrixValues3D_R8_Constant              17.66
              ModuleFunctions InterpolateLinearyMatrix3D                  8.31
              ModuleFunctions SetMatrixValues3D_R8_FromMatri              6.82
              ModuleFunctions THOMAS_3D                                   0.67
```

Parallelizing only the moduleFunctions yields a localized gain in most of the parallelized subroutines, except for the THOMASZ and the THOMAS\_3D.

### ModuleGeometry

### ModuleMap
