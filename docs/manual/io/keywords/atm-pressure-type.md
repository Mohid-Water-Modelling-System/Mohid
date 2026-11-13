# ATM PRESSURE TYPE

## Overview

**ATM\_PRESSURE\_TYPE** defines the type of Atmospheric pressure field to include. It allows currently allows two flavours: [Atmospheric pressure](http://en.wikipedia.org/wiki/Mean_Sea_Level_Pressure) and [Mslp](http://en.wikipedia.org/wiki/Mean_Sea_Level_Pressure) or Mean Sea-Level pressure.

## Usage

- Activate Atmospheric pressure [keyword](atm-pressure.md) in the [module Hydrodynamic](../../water/modules/hydrodynamic.md),
- Active **ATM\_PRESSURE\_TYPE** keyword in the [module Hydrodynamic](../../water/modules/hydrodynamic.md),
- Include a new property in the [module Atmosphere](../../water/modules/atmosphere.md) called **atmospheric pressure** or

**mslp**, or both (though only one will actually be used). The new property can be an ASCII file, a constant value, an HDF file, ...

- That's it!

### Syntax

Describe here the correct syntax of the keyword in the configuration files:

```
ATM_PRESSURE_TYPE : VALUE
```

Where value can be

- **0** i.e. no atmospheric pressure to be used,
- **1** i.e. use **atmospheric pressure** in the [module Atmosphere](../../water/modules/atmosphere.md) configuration,
- **2** i.e. use **mslp** in the [module Atmosphere](../../water/modules/atmosphere.md) configuration file;

### Default value

By default, the atmospheric pressure field will used:

```
ATM_PRESSURE_TYPE : 1
```

### Modules

Here are the modules where the keyword is supposed to be inserted:

- [Module Hydrodynamic](../../water/modules/hydrodynamic.md)

Here some other related modules:

- [Module Atmosphere](../../water/modules/atmosphere.md)

### Dependent keywords

Here's a list of mandatory keywords to use with the **ATM\_PRESSURE\_TYPE** keyword:

- [ATM\_PRESSURE](atm-pressure.md)

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

- [ATM\_PRESSURE](atm-pressure.md)
- WIND
- WIND\_SMOOTH\_PERIOD
- RAMP
- INERTIAL\_PERIODS
- [Mohid Keywords Online DataBase](http://www.mohid.com/IHS)
