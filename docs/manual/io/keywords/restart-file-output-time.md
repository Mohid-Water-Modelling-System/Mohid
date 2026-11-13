# RESTART FILE OUTPUT TIME

## Overview

**RESTART\_FILE\_OUTPUT\_TIME** is a [keyword](keywords.md) which is used to create an image of the state of a given MOHID run. The image of the state is used as input when performing a continuation of calculation. Thus, it allows to start a new MOHID run with the [CONTINUOUS](continuous.md) keyword. The image of the state of a MOHID run is a set of up to 5 `.fin` files. One per Module. Each file contains the module name and eventually the date-time of the state coded in the filename. Example: `WaterProperties_1_20070516-120000.fin5`

## Usage

This keyword enables to read an array of numbers which specify the output times and frequency. The first number represents the time (in seconds) after the beginning of the simulation in which the first output will be written. The following numbers represent the second, third, fourth...(and so on) times (in seconds) after the beginning of the simulation in which the second, third, fourth...(and so on) outputs will be written. The last number represents the time step (in seconds) with which the following outputs will be made. The outputs can yield a date suffix of the type `YYYYMMDD-hhmmss.fin` if [RESTART\_FILE\_OVERWRITE](restart-file-overwrite.md) is set to false. Otherwise, by default, the outputs are overwritten. Thus, at the end of the run, only the last output is available.

### Syntax

```
RESTART_FILE_OVERWRITE : first_output_value second_output_value last_output_is_frequency
```

is to be inserted at the root of the file (i.e. not inside a begin/end block).

### Modules

Here is a list of modules in which this keyword should be used in the correspondent data files:

- [Module Hydrodynamic](../../water/modules/hydrodynamic.md)
- [Module WaterProperties](../../water/modules/waterproperties.md)
- [Module InterfaceSedimentWater](../../water/modules/interfacesedimentwater.md)
- [Module Turbulence](../../water/modules/turbulence.md) \*
- [Module Lagrangian](../../water/modules/lagrangian.md)

"\*" Only when the Module GOTM is used.

## Sample

This example will produce outputs in the initial instant of the simulation and afterwards with a time step of 3600 seconds.

```
RESTART_FILE_OUTPUT_TIME      : 0.  3600.
```

This example will produce an output in the initial instant of the simulation, a second output 21600 seconds after the beginning of the simulation, a third 86400 seconds after the beginning of the simulation and afterwards with a time step of 3600 seconds.

```
RESTART_FILE_OUTPUT_TIME      : 0.  21600. 86400. 3600.
```

In order to overwrite the outputs, thus keeping only the last output, insert the following keyword

```
RESTART_FILE_OVERWRITE   : 1
```

If that value is set to 0 then the name of the file would have the following format modulename\_runnumber\_date\_time.fin. In the other hand when this keyword had the value of 1 it would be named as modulename\_runnumber.fin. To find the corresponding date and time for the last saved instant, it should be checked the file Display\_Runnumber.log in the exe folder.

## See also

Here's a list of related keywords

- [RESTART\_FILE\_OVERWRITE](restart-file-overwrite.md)
- [CONTINUOUS](continuous.md)
- [OUTPUT\_TIME](output-time.md)
