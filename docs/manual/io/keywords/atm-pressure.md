# ATM PRESSURE

## Overview

**ATM\_PRESSURE** defines whether to include or not an reference level pressure field.

## Usage

- Activate Atmospheric pressure **keyword** in the [module Hydrodynamic](../../water/modules/hydrodynamic.md),
- Active [ATM\_PRESSURE\_TYPE](atm-pressure-type.md) keyword in the [module Hydrodynamic](../../water/modules/hydrodynamic.md),
- Include a new property in the [module Atmosphere](../../water/modules/atmosphere.md) called **atmospheric pressure** or

**mslp**, or both (though only one will actually be used). The new property can an ASCII file, a constant value, an HDF file, ...

- That's it!

### Syntax

Describe here the correct syntax of the keyword in the configuration files:

```
ATM_PRESSURE : VALUE
```

Where value can be

- **0** i.e. no reference level pressure field to be used,
- **1** i.e. use a reference level pressure field defined in [module Atmosphere](../../water/modules/atmosphere.md).

### Type

**ATM\_PRESSURE** is a logical/boolean parameter.

### Default value

By default, no atmospheric pressure will used:

```
ATM_PRESSURE : 0
```

### Modules

Here are the modules where the keyword is supposed to be inserted:

- [Module Hydrodynamic](../../water/modules/hydrodynamic.md)

Here some other related modules:

- [Module Atmosphere](../../water/modules/atmosphere.md)

### Dependent keywords

Here's a list of mandatory keywords to use with the **ATM\_PRESSURE** keyword:

- [ATM\_PRESSURE\_TYPE](atm-pressure-type.md)

## Sample

Snippet taken from the [module Hydrodynamic](../../water/modules/hydrodynamic.md) configuration file:

```
ATM_PRESSURE            : 1
ATM_PRESSURE_TYPE       : 2
WIND                    : 2
WIND_SMOOTH_PERIOD      : 432000.
RAMP                    : 1
INERTIAL_PERIODS        : 6

INITIAL_ELEVATION       : 1
INITIAL_ELEVATION_VALUE : 2.08
```

Snippets taken from the [module Atmosphere](../../water/modules/atmosphere.md) configuration file:

```
<beginproperty>
NAME                    : atmospheric pressure
UNITS                   : Pa
DESCRIPTION             : Meteo-IST MM5 atmospheric pressure
DEFAULTVALUE            : 0.
FILE_IN_TIME            : HDF
FILENAME                : ../../GeneralData/MM5/Portugal/MM5_20060729_Portugal.hdf5
OUTPUT_HDF              : 1
<endproperty>
```

```
<beginproperty>
NAME                    : mslp
UNITS                   : Pa
DESCRIPTION             : Meteo-IST MM5 mean sea-level pressure
DEFAULTVALUE            : 0.
FILE_IN_TIME            : HDF
FILENAME                : ../../GeneralData/MM5/Portugal/MM5_20060729_Portugal.hdf5
OUTPUT_HDF              : 1
<endproperty>
```

NOTE: Only one of the above properties will actually be used. They were written both for illustrative purposes only.

## See also

- [ATM\_PRESSURE\_TYPE](atm-pressure-type.md)
- WIND
- WIND\_SMOOTH\_PERIOD
- RAMP
- INERTIAL\_PERIODS
- [Mohid Keywords Online DataBase](http://www.mohid.com/IHS)
