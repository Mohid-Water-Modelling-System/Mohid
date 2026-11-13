# Basic OpenMP syntax

In this page is provided an introduction with the basic concepts in OpenMP programming. The OpenMP language reference manual ([http://www.openmp.org](http://www.openmp.org)) should be consulted for further details.

OpenMP processing is made by a set of threads: the **Master** thread and the **Workers** threads.

OpenMP parallelization is accomplished by defining **Parallel regions** and by creating the threads (one for each available core): the processing inside a parallel region is divided between the Workers and the Master, which then proceed in parallel (simultaneously), instead of the unique thread existent in unparallel processing.

Each thread can have a set of **Private variables**, which are affected only by this thread. The choice of the private variables, which are explicitly defined in programming, is a central part of the OpenMP programming. Should be made private all the variables whose values are altered by each thread processing and that affect other threads processing.

An obvious choice of these private variables are the loop variables when these are used to alter positions in matrixes or vectors in each iteration of the loop.

However it should be noted that in private variables the values are undefined in enter and exit of the parallel region and that by default these variables have no storage association with the variables outside the parallel region (this default behavior can, however, be altered by specific OpenMP directives).

The instructions in OpenMP are provided by a set of **Directives**. These refer to several actions such as the definition of Parallel regions, Work sharing and Data attributes. Directives are specific for the underlying programming language being used, either C or Fortran.

In Fortran directives are case insensitive. The basic syntax in as following:

```
sentinel directive [clause [[,] clause] ...]
```

**Sentinel** is !$OMP in either fixed or free format. The continuation of directives (from one code line to another) are according with the underlying language: & in Fortran case.

**Clauses** are used to specify additional information for the directive. Important clauses are:

- PRIVATE (private variables list);

- NOWAIT: specify that threads will not syncronize (i.e. wait for each other) at the end of a specific construct (e.g. a DO loop) within a parallel region; this is specified at the end of a parallel construct;

A parallel region is typically specified as:

```
!$OMP PARALLEL PRIVATE(PrivateVariable1,...,PrivateVariableN)
```

Inside a parallel region the work is distributed by the threads through **Work Sharing constructs**.

For every processing made inside a parallel region must be assigned a Work Sharing construct. If this is not verified execution errors can occur.

An important Work Sharing construct is the DO loop, specified as:

```
!$OMP DO
... (do loop in Fortran)
!$OMP END DO
```

In this construct the iterations of the enclosed Fortran DO loop are distributed by the threads.

The way the work is distributed over the threads is managed by the SCHEDULE clause of the DO construct. An important option of SCHEDULE is DYNAMIC: provided a Chunk number of iterations each thread will process this fixed number of iterations. After finishing a Chunk each thread will began another available Chunk till all iterations are completed:

```
!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
...
!$OMP END DO
```

The distribution of Chunk among the threads requires synchronization for each assignment. This causes a overhead that could be important. The DYNAMIC option is advisable when each iteration involves an amount of work which is not predictable. This could be the case when IF constructs are present inside the loop containg extra processing done only in specific cases. Also if the threads arrive to the DO loop at different times, e.g. if they come from a previous DO loop with a NOWAIT end clause.

When the amount of work required in each DO loop iteration is predictable and the same for all threads it is advisable to use the STATIC option of SCHEDULE. This is particularly useful for DO loops which are unique in a parallel region:

```
!$OMP DO SCHEDULE(STATIC, CHUNK)
...
!$OMP END DO
```

Important characteristics of the Work Sharing constructs are that they do not create new threads, must contemplate all the existing threads or none at all and there is no barrier on entry (any available thread is not required to wait for the others) but a barrier on exit exists. The barrier on exit can be removed by the referred NOWAIT clause.

Work Sharing directives can appear outside the lexical extent of a parallel region (sequential lines of code appearing inside the parallel region). If they are dynamically linked with a parallel region (e.g. they appear in processing a subroutine call) they are orphaned. If, however, they appear outside this dynamic connection with a parallel region (and also not in the lexical extent of the region) they are ignored and the enclosed code is performed by one thread only and no parallelization is processed.

Inside a Work Sharing construct can be defined a critical region, if it is conveninent that only one thread at a time processes a specific code. This can be used to avoid problems with input/output operations potentially occuring by several threads processing at same time (problems can occur because memory locations are being assessed at same time). This is commanded as following:

```
!$OMP CRITICAL [Name]
... (code to be processed sequentially)
!$OMP END CRITICAL [Name]
```

In this case, no barriers exist in the entry and exit and all threads process the enclosed code although one at a time. If more than one critial region is defined in the code then every critical region should have a different name or the execution outcome could become undeterminated. Caution should be given to the naming of critical sections: if these names do not conflict with variable names they do conflict with names of subroutines and common blocks.

Within a Work Sharing construct can also be specified that a portion of the code is processed only by one thread, e.g. to read information from a file common to all threads. This is done with the following syntax:

```
!$OMP SINGLE [clause ...]
... (code to be processed by only one thread)
!$OMP END SINGLE
```

When one wants this single thread to be the Master then following syntax is used:

```
!$OMP MASTER
... (code to be processed by only Master thread)
!$OMP END MASTER
```

As with CRITICAL, no barriers exist in SINGLE/MASTER in the entry and exit.

In these contructs the fact that no barrier exists at the exit may cause problems if the single thread processed code section is intended to be dealt with previously from the subsequent code. In this situation a barrier can be introduced at the end of the SINGLE/MASTER construct:

```
!$OMP BARRIER
```

Under this directive threads wait at the barrier point till all threads reach it.

An important rule about the BARRIER use is that this directive cannot be nested in a Work Sharing construct such a DO construct.

Sometimes some variables used by threads cannot be made private. One notorious case consists in variables performing sums of values obtained in iterations of a cycle, which when summing integers are sometimes called «counters». As these variables should accumulate values obtained in each iteration these cannot be made private. This would imply that a critical section should be introduced for the actualization of the variable which would make the code sequential. To avoid this there is an OpenMP clause called REDUCTION which when applied to a DO Work Sharing construct has the following syntax:

```
!$OMP DO SCHEDULE(DYNAMIC, CHUNK) REDUCTION(+ : Counter)
```

Where in this example «Counter» is the reduction variable.

The REDUCTION operation has as restriction that the reduction variable only is actualized once inside the construct and must be a shared variable.

In practice in the program execution a copy of the variable is made for each thread, which actualizes it as a local variable, and in the end of the construct the thread variables are added to the shared variable from which they originate. The knowledge of this process is important because in the case of the accumulation of real values this process can originate differences relative to the unparallel result due to rounded values.

***Parallelization overheads:***

Although parallelization involves potentially gains in computer resources use it also involves overheads which can sometimes be important.

The creation of the thread team (Workers) at the beginning of each parallel region is a source of such overheads. Because of this it is preferrable in each program/subrotine to be parallelized to create only one parallel region and then deal with the code that must be done only by one thread with OpenMP directives than to create several parallel regions.

Another source of overhead is the synchronization between threads in Work Sharing constructs.

These overheads cause that in some times the parallelization solution may not be performing better than the non parallelization solution, especially if the computational burden of the problem is low.

This is the case of DO constructs when each cycle iteration has a small computational burden (e.g. equaling one matrix's value to other matrix's value) and is to be expected in nested loops when the parallelized DO loop is inner.

Generally as the scale of the problem increases the overheads will be less significant and parallelization more advantageous.

In very small scale problems parallelization may provide more comsumption of computational resources than the non parallelization.

***Race condition:***

In some case execution of a parallel program may be suspended because different threads are competing for the actualization of some variable. This is called a race condition.

The programmer should take care to avoid this situation as its debugging is very difficult. Most commonly the execution will be suspended and the processors would become idle without any error message. When a vast part of the code is parallelized it is then almost impossible to find the cause of the error unless one reckons possible race condition situations.

A way to avoid this is to use critical regions or a REDUCTION clause.

***Overheads and thread efficiency***

http://i.imgur.com/aZXGh.png

http://i.imgur.com/ff08Z.png

http://i.imgur.com/F1kLM.png

http://i.imgur.com/PEmH4.png

***References:***

- Chandra, Rohit, 2001, Parallel Programming in OpenMP, Academic Press.
- [Chapman, Barbara, Using OpenMP, Portable Shared Memory Parallel Programming, The MIT Press](http://drop.io/8gzzkv3).
