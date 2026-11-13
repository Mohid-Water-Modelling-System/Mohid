# Mohid Programming

MOHID is programmed in ANSI FORTRAN 95 using an [object oriented approach](README.md#object-oriented-approach#object-oriented-approach) and [hierarchical structure](../../water/modules/modules.md) which enables consisting of several [modules](../../water/modules/modules.md) to develop several numerical models and [support tools](../../tools/mohid-support-tools.md) with a high level of code re-use. To this structure we call [Mohid Framework](README.md#mohid-framework#mohid-framework).

Due to the high number of developers programming in MOHID, several [programming guidelines](programming-guidelines.md) have been defined in order to [change the code](changing-the-code.md) in a straighforward, systematic and sustentable way.

MOHID is licensed under GPL and the code can be downloaded from the download area. Currently, MOHID has been [compiled](../build/compiling.md) and run under [several architectures](supported-architectures.md).

## Main Program

MOHID models have a main program that is responsible for starting the computation, calling the Main Module Constructor, updating time and calling the Main Module Modifier and after the time has come to final step, calling the Main Module Destructor.

## Modules

The Main Module is the responsible for calling the construction, modify and destruction of all others. In case of MOHID Land is also the resposible for exchaging information between MOHID Land Modules. In case of MOHID Water there are Interface Modules created for this.

Every Module is responsible to solve specific processes and has the following structure:

```
Constructor where the data files are read and matrixes are created and initiaalized. This is called only once at the construction phase.
Modifier that is called each time the time changes and is resposible for updating the matrixes for the new time step.
Destructor that is where all the matrixes are released and connections between modules destroyed. This is called only once afte the simulation time is finished.
Selector that is how other modules can get the matrixes values from the Module
```

## Sequence of Call

The first phase is the construction phase and the Main Program constructor calls the Main Module consstructor that calls the other Modules constructors (one at a time, the next is only called when the last has finished). On each Module all the routines inside the constructor are called, one by one, in sequence.

When all the Modules are constructed, is the Main Program Modifier that starts and the time starts and it calls the Main Module tath calls the other Modules modifiers. Once again they are called in sequence and only one finishes the other starts. When all the Modules have finished the fist time step Modifier, the Main Program updates time for the next time and so on.

When all the time steps have been done, the Main Program calls the Main Program Destrutor that calls each Module Destructor...

## Routines

Each Module has a set of routines that are the responsible for the work done. The routines that are called in the Constructor are only called onde at the start of the model when the time was not yet started and so on.

## Variables

The variales defined inside each routine are only local variales and are only accessible inside the routine.

The variables that start by Me%... are global variables that are accesible by all routines. The Me%.. is defined at the start of the files by the Types and all the types exist to organize variables by type or origin. When a Type is referenced inside another type it becomes a subtype of the first type (e.g. Me%Type1%Type2%...)

## Programming

To add stuff to the model, the user needs to be aware of the above structure and how is the sequence of calls and to add stuff needs to add routines or change routines by adding variables (local and/or global). Follow the referred rules for programming in MOHID!

The user has at his diposition the two basic features that exist in every programming language:

```
"do loops" that loop trough the matrixes to compute work for each cell, or loop trough blocks for example when readng data files
"if statements" and "select case" to distinguish betwwen options
```

## Related links

See also [Parallel processing](parallel-processing.md), [Compiling](../build/compiling.md), [Profiling](../../archive/linux/profiling.md), [Programming in Fortran issues](programming-in-fortran-issues.md)
