# OUTPUT TIME

## Overview

OUTPUT\_TIME is a [keyword](keywords.md) which is common to several MOHID [modules](../../water/modules/modules.md). It is normally used to identify the frequency of results outputs in [HDF5](../../tools/hdf5.md) format.

## Usage

This keyword enables to read an array of numbers which specify the output times and frequency. The first number represents the time (in seconds) after the beginning of the simulation in which the first output will be written. The following numbers represent the second, third, fourth...(and so on) times (in seconds) after the beginning of the simulation in which the second, third, fourth...(and so on) outputs will be written. The last number represents the time step (in seconds) with which the following outputs will be made.

## Sample

This example will produce outputs in the initial instant of the simulation and afterwards with a time step of 3600 seconds.

```
OUTPUT_TIME      : 0.  3600.
```

This example will produce an output in the initial instant of the simulation, a second output 21600 seconds after the beginning of the simulation, a third 86400 seconds after the beginning of the simulation and afterwards with a time step of 3600 seconds.

```
OUTPUT_TIME      : 0.  21600. 86400. 3600.
```

## Other information

Here is a list of modules in which this keyword can be used in the correspondent data files:

- [Module Hydrodynamic](../../water/modules/hydrodynamic.md)
- [Module WaterProperties](../../water/modules/waterproperties.md)
- [Module Turbulence](../../water/modules/turbulence.md)
- [Module Lagrangian](../../water/modules/lagrangian.md)
- [Module Assimilation](../../water/modules/assimilation.md)
- [Module InterfaceWaterAir](../../water/modules/interfacewaterair.md)
- [Module InterfaceSedimentWater](../../water/modules/interfacesedimentwater.md)
- [Module Atmosphere](../../water/modules/atmosphere.md)
- [Module SedimentProperties](../../water/modules/sedimentproperties.md)
- [Module Consolidation](../../water/modules/consolidation.md)
- [Module Waves](../../water/modules/waves.md)
- [Module Sand](../../water/modules/sand.md)

## See also

[SURFACE\_OUTPUT\_TIME](surface-output-time.md)
