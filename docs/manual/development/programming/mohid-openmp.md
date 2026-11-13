# Mohid OpenMP

This wiki entry is about the project to parallelize Mohid with [OpenMP directives](parallel-processing.md).

## How to pre-set the number of threads to use

### Windows

In the Module model it has to be added the following keyword

```
OPENMP_NUM_THREADS            :  N
```

where N is the number of threads you want to use. It is recommend that N be equal to, or less than, the number of physical cores of your machine. (For example, in a quad core machine, like the intel i7 970, N = 4 gives the best performance.)

The you can launch MOHIDWater

```
exe> MohidWater.exe
```

## Official benchmarking

Benchmarks must be identified with the processor type (intel, amd), the architecture (32bit, 64bit), the number of cores (#2, #4) and the number of threads (#4, #8, ...). The simulation must be given an id such as PCOMS, PCOMS-bio, Azores, etc... The benchmarks must be performed with the data contained in the MOHID Outwatches logs, once with the Model run without openmp, and once with the model run with openmp.

### Benchmarking formula

The benchmark formula is:

```
Speedup = A / B, where A is time with openmp and B is time without openmp. If Speedup < 1, then speed is up; if Speedup > 1, then speed is down.
```

### Official results

#### PCOMS 20100801 quad-core 4 threads

Module SubRoutine Time w/o Time w Speed-up

```
ModuleModel   RunOneModel  727.679     193.4     0.87
```

#### Azores 20100801 six-core 12 threads

Module SubRoutine Time w/o Time w Speed-up

```
ModuleModel   RunOneModel  8792.656	7627.37   0.27
```

Obs: the Module Life single loop was successfully parallelized.

## Other links

- [Parallel processing](parallel-processing.md)
