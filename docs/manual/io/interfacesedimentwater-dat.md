# InterfaceSedimentWater.dat

```
OUTPUT_TIME               : 0. 3600.
OUTPUT_HDF                : 1
```

```
TIME_SERIE                : 1
TIME_SERIE_LOCATION       : ..\General Data\TimeSeriesLocation.dat
```

```
CONSOLIDATION : 1
```

```
COMPUTE_RUGOSITY : 0
```

```
<begin_rugosity>
DEFAULTVALUE                : 0.000833
<end_rugosity>
```

```
SHEAR_STRESS_METHOD           : 2
WAVETENSION                   : 0
OUTPUT_SHEAR_STRESS           : 1
```

```
<begin_critical_shear_erosion>
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 0.2
<end_critical_shear_erosion>
```

```
<begin_critical_shear_deposition>
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 0.1
<end_critical_shear_deposition>
```

```
<begin_erosion_rate>
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 0.00005
<end_erosion_rate>
```

```
<beginproperty>
NAME                      : cohesive sediment
UNITS                     : kg/m2
DESCRIPTION               : cohesive sediment mass
IS_PARTICULATE            : 1
WATER_FLUXES              : 1
SEDIMENT_FLUXES           : 1
SEDIMENT_WATER_FLUXES	   : 1
EROSION                   : 1
DEPOSITION                : 1
MASS_LIMITATION           : 1
MASS_MIN                  : 0.
OLD                       : 0
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 0.
TIME_SERIE                : 1
BOX_TIME_SERIE            : 0
OUTPUT_HDF                : 1
<endproperty>
```

```
<beginproperty>
NAME                      : SAND_1
UNITS                     : kg/m2
DESCRIPTION               : SAND CLASS
DYNAMIC                   : 1
IS_PARTICULATE            : 1
NON_COHESIVE	           : 1
WATER_FLUXES              : 1
SEDIMENT_FLUXES           : 1
SEDIMENT_WATER_FLUXES	   : 1
EROSION                   : 1
DEPOSITION                : 1
MASS_LIMITATION           : 1
MASS_MIN                  : 0.
OLD                       : 0
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 0
TIME_SERIE                : 1
BOX_TIME_SERIE            : 0
OUTPUT_HDF                : 1
<endproperty>
```

```
<beginproperty>
NAME                      : SAND_2
UNITS                     : kg/m2
DESCRIPTION               : SAND CLASS
DYNAMIC                   : 1
IS_PARTICULATE            : 1
NON_COHESIVE	           : 1
WATER_FLUXES              : 1
SEDIMENT_FLUXES	   : 1
SEDIMENT_WATER_FLUXES	   : 1
EROSION                   : 1
DEPOSITION                : 1
MASS_LIMITATION           : 1
MASS_MIN                  : 0.
OLD                       : 0
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 0
TIME_SERIE                : 1
BOX_TIME_SERIE            : 0
OUTPUT_HDF                : 1
<endproperty>
```

```
<beginproperty>
NAME                      : SAND_3
UNITS                     : kg/m2
DESCRIPTION               : SAND CLASS
DYNAMIC                   : 1
IS_PARTICULATE            : 1
NON_COHESIVE              : 1
WATER_FLUXES              : 1
SEDIMENT_FLUXES	   : 1
SEDIMENT_WATER_FLUXES	   : 1
EROSION                   : 1
DEPOSITION                : 1
MASS_LIMITATION           : 1
MASS_MIN                  : 0.
OLD                       : 0
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 0
TIME_SERIE                : 1
BOX_TIME_SERIE            : 0
OUTPUT_HDF                : 1
<endproperty>
```
