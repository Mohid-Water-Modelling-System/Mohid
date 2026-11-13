# Module Assimilation

## Overview

Module Assimilation is responsible for handling information used in [relaxation](../relaxation.md) schemes. It mainly reads the reference fields and the relaxation coefficients fields used to relax the solution to the reference solution. The standard relaxation scheme implemented in MOHID follow the following equation:

$\frac{\partial \phi}{\partial t} = \frac{\phi_{ref} - \phi}{\partial \tau}$

where $\phi$ is the relaxed variable computed by the model, $\phi_{ref}$ is the reference solution and $\tau$ relaxation coefficient.

[Module Hydrodynamic](hydrodynamic.md) and [Module WaterProperties](waterproperties.md) use this module to read and prepare the information. The relaxation algorithms are actually programmed in Module Hydrodynamic and Module WaterProperties, making Module Assimilation somewhat a "reading" module. This type of approach is normally used when [nesting models](../../reference/nesting-models.md) or using a climatological reference solution.

This module relates to simplified data assimilation schemes. For information regarding advanced data assimilation in MOHID please check ModuleSequentialAssimilation.

## Concepts

### Property

A property consists of a reference field and a relaxation coefficients field. It can be 2D (e.g. water level) or 3D (e.g. temperature).

### Reference field

The reference field is a prescribed solution which is used to force the solution computed by the model. The reference field must be given in the same units as the computed field. Examples of properties are velocity U, velocity V, water level, temperature, salinity, etc. A reference field must be defined in the U, V or Z grid, via keyword **TYPE\_ZUV** and it can be set using [Module FillMatrix](../../base/module-fillmatrix.md).

### Relaxation coefficients

The [relaxation](../relaxation.md) coefficients are used to define the *rate* at which the solution computed by the model is relaxed to the prescribed solution. The units of these coefficients are seconds.

To generate in a easy way a relaxation coefficients file use the [AssimilationZones](../../tools/assimilationzones.md) tool.

## Other features

### Cold relaxation period

- **COLD\_RELAX\_PERIOD**
- *COLD\_ORDER*

### Outputs

## User manual

### Assimilation data file

A property is identified by the [keyword](../../io/keywords/keywords.md) **NAME** inside a [block](../../io/keywords/blocks.md) delimited by the *<beginproperty>* and *<endproperty>* tags. The dimension of the property is specified by keyword **DIMENSION** (e.g. for a 3D property: ***DIMENSION  : 3***).

A reference field is defined inside a block delimited by the *<<begin\_field>>* and *<<end\_field>>* tags.

A relaxation coefficients field is defined inside a block delimited by the *<<begin\_coef>>* and *<<end\_coef>>* tags.

A relaxation coefficients and reference fields must be defined in the U, V or Z grid, via keyword **TYPE\_ZUV** and can be set using [Module FillMatrix](../../base/module-fillmatrix.md).

```
<beginproperty>
NAME                    : velocity U
UNITS                   : m/s
DIMENSION               : 3
OUTPUT_HDF              : 1

<<begin_field>>
DEFAULTVALUE            : 0
INITIALIZATION_METHOD   : HDF
FILE_IN_TIME            : HDF
FILENAME                : ../../GeneralData/Reference_vel_U.hdf5
TYPE_ZUV                : z
<<end_field>>

<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : u
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : ASCII_FILE
FILENAME                : ../../GeneralData/RelaxationCoefs_U.dat
<<end_coef>>

<endproperty>

<beginproperty>
NAME                    : velocity V
UNITS                   : m/s
DIMENSION               : 3
OUTPUT_HDF              : 1

<<begin_field>>
DEFAULTVALUE            : 0
INITIALIZATION_METHOD   : HDF
FILE_IN_TIME            : HDF
FILENAME                : ../../GeneralData/Reference_vel_V.hdf5
TYPE_ZUV                : z
<<end_field>>

<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : v
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : ASCII_FILE
FILENAME                : ./../GeneralData/RelaxationCoefs_V.dat

<<end_coef>>

<endproperty>
```

```
<beginproperty>
NAME                    : temperature
UNITS                   : ºC
DIMENSION               : 3D
OUTPUT_HDF              : 1

<<begin_field>>
DEFAULTVALUE            : 0
INITIALIZATION_METHOD   : HDF
FILE_IN_TIME            : HDF
FILENAME                : ../../GeneralData/Reference_T.hdf5
TYPE_ZUV                : z
<<end_field>>

<<begin_coef>>
DEFAULTVALUE            : 1e32
TYPE_ZUV                : z
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : ASCII_FILE
FILENAME                : ../../GeneralData/RelaxationCoefs_T.dat

<<end_coef>>

<endproperty>
```
