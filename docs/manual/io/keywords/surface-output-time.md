# SURFACE OUTPUT TIME

## Overview

SURFACE\_OUTPUT\_TIME is a [keyword](keywords.md) which is used to identify the frequency of results outputs in [HDF5](../../tools/hdf5.md) format **for the surface layer**. It must be defined in the [Hydrodynamic](../../water/modules/hydrodynamic.md) file.

## Usage

This keyword enables to read an array of numbers which specify the output times and frequency. The first number represents the time (in seconds) after the beginning of the simulation in which the first output will be written. The following numbers represent the second, third, fourth...(and so on) times (in seconds) after the beginning of the simulation in which the second, third, fourth...(and so on) outputs will be written. The last number represents the time step (in seconds) with which the following outputs will be made.

## Sample

This example will produce outputs in the initial instant of the simulation and afterwards with a time step of 3600 seconds.

```
SURFACE_OUTPUT_TIME      : 0.  3600.
```

This example will produce an output in the initial instant of the simulation, a second output 21600 seconds after the beginning of the simulation, a third 86400 seconds after the beginning of the simulation and afterwards with a time step of 3600 seconds.

```
SURFACE_OUTPUT_TIME      : 0.  21600. 86400. 3600.
```

For every property in [WaterProperties](../../water/modules/waterproperties.md) file that you want to surface output into the HDF, define the following keyword:

```
OUTPUT_SURFACE_HDF   : 1
```

## Modules

Here is a list of modules in which this keyword can be used in the correspondent data files:

- [Module Hydrodynamic](../../water/modules/hydrodynamic.md)
- [Module WaterProperties](../../water/modules/waterproperties.md)

## See also

Here's a list of related keywords

- [OUTPUT\_TIME](output-time.md)
- [RESTART\_FILE\_OUTPUT\_TIME](restart-file-output-time.md)
- [RESTART\_FILE\_OVERWRITE](restart-file-overwrite.md)
