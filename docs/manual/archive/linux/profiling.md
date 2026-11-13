# Profiling

Is is possible to profile code compiled with the Intel Fortran compiler using gprof [

$$
1
$$

](http://www.cs.utah.edu/dept/old/texinfo/as/gprof_toc.html). Intel Fortran's current release (9.1) requires Linux to generate profiling information.

## Quick dive

To generate profiling information follow these steps:

1. ifort -p foo.F90 -o foo
2. ./foo
3. gprof foo

The first step will generate a executable with profiling information. At the end of the second step a file named gmon.out will be writen alongside the executable file; this file contains the profiling information generated during the program run. The third step screens the profile.

Note: the -p option must be issued both to compile and link the program.

## Profiling Mohid

To profile Mohid open the Makefiles used on [\*NIX platforms](nix-platforms.md). Find variable CC in every makefile and replace "ifort" with "ifort -p". make clean, make and execute Mohid. At the end of the run a file named gmon.out will be created in the executable's folder. Explore gmon.out issuing the command gprof Mohid.

## Usefull links

NCSA [

$$
2
$$

](http://www.ncsa.uiuc.edu/UserInfo/Resources/Hardware/IBMp690/Doc/Timing_Profiling.html).

GNU prof [

$$
3
$$

](http://www.cs.utah.edu/dept/old/texinfo/as/gprof_toc.html).
