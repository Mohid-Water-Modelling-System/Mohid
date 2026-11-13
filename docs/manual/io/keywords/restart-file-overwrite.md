# RESTART FILE OVERWRITE

## Overview

**RESTART\_FILE\_OVERWRITE** is a [keyword](keywords.md) which is used in conjunction with [RESTART\_FILE\_OUTPUT\_TIME](restart-file-output-time.md) to create an image of the state of a given MOHID run. The image of the state of a MOHID run is a set of up to 5 `.fin` files. One per Module. Each file contains the module name and eventually the date-time of the state coded in the filename. Example: `WaterProperties_1_20070516-120000.fin5`. This keyword in particular states whether to keep a record of every state output by appending a date-time suffix, whether to keep only the last output with no date-time suffix. A date-time suffix is composed like this: *YYYYMMDD-hhmmss*.

## Usage

Simply add in every file of the required modules (see list below) the keyword with the true or false value.

### Syntax

```
RESTART_FILE_OVERWRITE : first_output_value second_output_value last_output_is_frequency
```

is to be inserted at the root of the file (i.e. not inside a begin/end block).

### Default value

By default, the outputs will be overwritten.

```
RESTART_FILE_OVERWRITE : 1
```

### Modules

Here is a list of modules in which this keyword should be used in the correspondent data files:

- [Module Hydrodynamic](../../water/modules/hydrodynamic.md)
- [Module WaterProperties](../../water/modules/waterproperties.md)
- [Module InterfaceSedimentWater](../../water/modules/interfacesedimentwater.md)
- [Module Turbulence](../../water/modules/turbulence.md)
- [Module Lagrangian](../../water/modules/lagrangian.md)

### Dependent keywords

Here's a list of mandatory keywords

- [RESTART\_FILE\_OUTPUT\_TIME](restart-file-output-time.md)

## Sample

This example will produce outputs in the initial instant of the simulation and afterwards with a time step of 3600 seconds. It will also add a time suffix for each outputs.

```
RESTART_FILE_OUTPUT_TIME: 0.  3600.
RESTART_FILE_OVERWRITE  : 0
```

## See also

Here's a list of related keywords

- [RESTART\_FILE\_OUTPUT\_TIME](restart-file-output-time.md)
- [CONTINUOUS](continuous.md)
- [OUTPUT\_TIME](output-time.md)
