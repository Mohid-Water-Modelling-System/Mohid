# Visual Studio Configuration

In Visual Studio, a configuration contains all the compiler and linker settings pre-saved, in order to build a solution (compile + link). It also may contain pre and post-build steps, such as copying the target executable.

## Mohid configuration policy

Mohid configurations used to be only "Release" and "Debug". But soon, they multiplied with different options, such as double, openmp, mpi, openmi.

This means that each configuration MOHIDWater or MOHIDLand (or any other MOHID tool) executable will behave differently. In order to allow the user of pre-built binaries, each executable will have a suffix that identifies which platform/configuration was it that built it.

Hence, here's the current list of configuration/platforms suffixes:

- debug\_single
- debug\_double
- debug\_single\_openmp
- debug\_double\_openmp
- debug\_single\_openmi
- debug\_double\_phreeqc
- release\_single
- release\_double
- release\_single\_openmp
- release\_double\_openmp
- release\_double\_openmp\_phreeqc
- x64\_debug\_single
- x64\_debug\_double
- x64\_debug\_single\_openmp
- x64\_debug\_double\_openmp
- x64\_debug\_single\_openmi
- x64\_debug\_double\_phreeqc
- x64\_release\_single
- x64\_release\_double
- x64\_release\_single\_openmp
- x64\_release\_double\_openmp
- x64\_release\_double\_openmp\_phreeqc
