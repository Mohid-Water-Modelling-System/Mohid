# New: MOHID Water Vectorial and Angle Properties

## Overview

This is a guide on configuring [MOHID Water](../water/mohid-water.md) new way of defining vectorial properties. And a recall on angle properties computation.

### Referentials

Input Vectorial properties in MOHID have a X (zonal component) and Y (meridian component).

Input Angle properties in MOHID in ModuleAtmosphere and ModuleWaves have nautical referential (0º is wind or wave from N towards S, 90º is from East towards West and so on).

MOHID computes fluxes using components in cell referential (aligned with U and V faces) and angles are also refered to cell trigonometric circle.

## Old Way of definig vectorial properties and angle propeties

MOHID treated until now (Dez 2015) vectorial properties as two scalar properties, defined separately.

This an example of wind velocity X and wind velocity Y in Atmosphere\_X.dat file until now:

```
<beginproperty>
NAME                      : wind velocity X
UNITS                     : m/s
DESCRIPTION               : wind velocity X
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : -1.0
REMAIN_CONSTANT           : 1
OUTPUT_HDF                : 0
<endproperty>

<beginproperty>
NAME                      : wind velocity Y
UNITS                     : m/s
DESCRIPTION               : wind velocity Y
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 0.0
REMAIN_CONSTANT           : 1
OUTPUT_HDF                : 0
<endproperty>
```

This approach is straightforward to implement as Module FillMatrix treats all properties the same, reads the options, fills the matrix and returns to calling modules.

However it creates one problem: how to deal with vectorial transformation from user input referential to the MOHID cell referential (if grid rotation and distortion exists).

The solution used until now was that each module would do its own transformations what could lead to conversion errors or lack of conversions.

The same issue occurs with angle properties where each module would do its own transformations what could lead to same problems.

## The New Way of definig vectorial properties

The New age has arrived!

Only vectorial propeties have changed the input format. Angle properties are scalar properties so in terms of input there was no change.

In order to avoid the descentralization and repetition of the transformation process around modules, the vectorial property concept was created in Module FillMatrix.

And both vectorial and angle tranformation are made inside Module FillMatrix and the resulting matrixes returned to modules are always in MOHID cell referential.

An example of the new vectorial input is in Module Atmosphere:

```
<beginproperty>
NAME                      : wind velocity
UNITS                     : m/s
DESCRIPTION               : wind velocity
DEFAULTVALUE              : -1.0 0.0
REMAIN_CONSTANT           : 1
OUTPUT_HDF                : 0
<endproperty>
```

That replaces the two properties wind velocity X and wind velocity Y.

The keyword DEFAULTVALUE can now read two or three components (if 3D property) istead of only one as before.

### What is the Input/Output conventions used by all Modules?

Basically user inputs vectorial properties in zonal/meridional components and angles in nautical referential (atmosphere and waves) and MOHID in all modules uses internaly MOHID cell referential given by Module FillMatrix (after transformation from user referential).

To write output each module call again Module FillMatrix to get the fields in user input referential. So input and output is always in the zonal/meridional components (vectorial) and angles in nautical referential (atmosphere and waves)

### What vectorial properties does the new approach changes in terms of data file definition?

Module Atmosphere

1. wind velocity X and wind velocity Y are now only defined by wind velocity property

Module InterfaceWaterAir

1. wind stress X and wind stress Y are now only defined by wind stress property

Module Waves

1. radiation stress X and radiation stress Y are now only defined by radiation stress property

### What kind of input type is allowed for the new vectorial properties?

for now only constant, timeserie, ascii\_file and hdf inputs are allowed. However new developments are welcome for other types of input.

### How does ModuleFillMatrix now works?

Module FillMatrix now behaves similarly as before where each property could have only one field of type time serie, ascii file, hdf. Now it can have more (two in 2D properties). The Module FillMatrix searches how many components each property has and does the same old process on all of them, but saving the results in different fields (e.g. component X and Y).

The modules that ask properties fields to fillmatrix also save information in different matrixes and calls to fillmatrix are very similar as before but with one more field (the second component).

### Remark

Now the keyword DEFAULTVALUE in a property can read more than one value so DO NOT add text comments or other things in the same line after the defined values

This is correct inside a 2D property:

```
 !comments go here
 DEFAULTVALUE              : -1.0 0.0
```

```
 DEFAULTVALUE              : -1.0 0.0   !this is not correct commenting and can stop the model at construction phase
```

## Examples

### Module Atmosphere wind velocity defined by ASCII grid

In case of ascii grid the user needs to provide two different files and two default values for x and y component (in this order).

```
<beginproperty>
NAME                      : wind velocity
UNITS                     : m/s
DESCRIPTION               : wind velocity
INITIALIZATION_METHOD     : ASCII_FILE
FILENAME_X                : ..\General Data\Boundary Conditions\ASCIIGrid_-1.dat
FILENAME_Y                : ..\General Data\Boundary Conditions\ASCIIGrid_0.dat
DEFAULTVALUE              : -1.0 0.0
REMAIN_CONSTANT           : 1
TIME_SERIE                : 1
OUTPUT_HDF                : 1
<endproperty>
```

### Module Atmosphere wind velocity defined by time series

In case of time serie the user can

1) provide one file but two different collumns and two default values for x and y component (in this order).

```
<beginproperty>
NAME                      : wind velocity
UNITS                     : m/s
DESCRIPTION               : wind velocity
FILE_IN_TIME              : TIMESERIE
FILENAME                  : ..\General Data\Boundary Conditions\Velocity.srm
DATA_COLUMN_X             : 2
DATA_COLUMN_Y             : 3
DEFAULTVALUE              : -1.0 0.0
REMAIN_CONSTANT           : 1
TIME_SERIE                : 1
OUTPUT_HDF                : 1
<endproperty>
```

OR 2) can provide two files and respectve collumns and two default values for x and y component (in this order).

```
<beginproperty>
NAME                      : wind velocity
UNITS                     : m/s
DESCRIPTION               : wind velocity
FILE_IN_TIME              : TIMESERIE
FILENAME_X                : ..\General Data\Boundary Conditions\VelX.srm
FILENAME_Y                : ..\General Data\Boundary Conditions\VelY.srm
DATA_COLUMN_X             : 2
DATA_COLUMN_Y             : 2
DEFAULTVALUE              : -1.0 0.0
REMAIN_CONSTANT           : 1
TIME_SERIE                : 1
OUTPUT_HDF                : 1
<endproperty>
```

### Module Atmosphere wind velocity defined by HDF

In case of hdf the user can

1) provide one file but two different hdf field names and two default values for x and y component (in this order).

```
<beginproperty>
NAME                      : wind velocity
UNITS                     : m/s
DESCRIPTION               : wind velocity
FILE_IN_TIME              : HDF
FILENAME                  : ..\General Data\Boundary Conditions\Atmosphere_11.hdf5
HDF_FIELD_NAME_X          : wind velocity X
HDF_FIELD_NAME_Y          : wind velocity Y
DEFAULTVALUE              : -1.0 0.0
REMAIN_CONSTANT           : 1
TIME_SERIE                : 1
OUTPUT_HDF                : 1
<endproperty>
```

OR 2) can provide two files and respectve hdf field names and two default values for x and y component (in this order).

```
<beginproperty>
NAME                      : wind velocity
UNITS                     : m/s
DESCRIPTION               : solar radiation
FILE_IN_TIME              : HDF
FILENAME_X                : ..\General Data\Boundary Conditions\WindX.hdf5
FILENAME_Y                : ..\General Data\Boundary Conditions\WindY.hdf5
HDF_FIELD_NAME_X          : wind velocity X
HDF_FIELD_NAME_Y          : wind velocity Y
DEFAULTVALUE              : -1.0 0.0
REMAIN_CONSTANT           : 1
TIME_SERIE                : 1
OUTPUT_HDF                : 1
<endproperty>
```

### Module InterfaceWateAir wind stress defined to be computed

This has the same options as shown in atmosphere ones so only a block is presented to ilustrate, in this case as a computed property.

```
<beginproperty>
NAME                      : wind stress
UNITS                     : N/m2
DESCRIPTION               : wind stress
DEFAULTVALUE              : 0.0 0.0
REMAIN_CONSTANT           : 0
TIME_SERIE                : 0
OUTPUT_HDF                : 0
<endproperty>
```

### Module Waves wave stress defined constant

This has the same options as shown in atmosphere ones so only a block is presented to ilustrate, in this case a constant one.

```
<begin_radiationstress>
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : -1.0 0.0
REMAIN_CONSTANT           : 1
TIME_SERIE                : 1
OUTPUT_HDF                : 1
<end_radiationstress>
```
