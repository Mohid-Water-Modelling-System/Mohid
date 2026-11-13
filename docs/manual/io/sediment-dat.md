# Sediment.dat

```
IN_SEDIMENT               : ..\General Data\Sediment_Column_Thickness.dat
```

```
OUTPUT_TIME               : 0. 600.
RESIDUAL                  : 0
```

```
TIME_SERIE                : 1
TIME_SERIE_LOCATION       : ..\General Data\TimeSerie.dat
```

```
CONTINUOUS                : 0
```

```
BATHYMETRY_EVOLUTION      : 1
MORPHOLOGICAL_FACTOR      : 100.
```

```
BEDLOAD_METHOD            : 1
BEDSLOPE                  : 1

GEOMETRY_EVOLUTION        : 1
MIN_THICKNESS             : 0.1
MAX_THICKNESS             : 0.5
```

```
CONSOLIDATION             : 1
CONSOLIDATION_RATE1       : 0.9
CONSOLIDATION_RATE2       : 0.002
```

```
MAX_DENSITY_COHESIVE                  : 650.
MIN_DENSITY_COHESIVE                  : 150.
PERMEABILITY_EFFECTIVE_STRESS_DENSITY : 220.
```

```
<begincohesive>
NAME                      : cohesive sediment
UNITS                     : %
DESCRIPTION               : cohesive sediment
FILE_IN_TIME              : NONE
INITIALIZATION_METHOD     : LAYERS
LAYERS_VALUES             : 0. 0. 0. 1. 0.
DEFAULTVALUE              : 0.
D50                       : 4e-6
<begindrydensity>
INITIALIZATION_METHOD     : CONSTANT
DEFAULTVALUE              : 150.
<enddrydensity>
<endcohesive>
```

```
<beginsandclass>
NAME                      : SAND_1
UNITS                     : %
DESCRIPTION               : SAND_1
DYNAMIC                   : 1
FILE_IN_TIME              : NONE
INITIALIZATION_METHOD     : LAYERS
LAYERS_VALUES             : 0. 0. 1. 0. 0.
DEFAULTVALUE              : 0.
D50                       : 0.2e-3
<endsandclass>
```

```
<beginsandclass>
NAME                      : SAND_2
UNITS                     : %
DESCRIPTION               : SAND_2
DYNAMIC                   : 1
FILE_IN_TIME              : NONE
INITIALIZATION_METHOD     : LAYERS
LAYERS_VALUES             : 0. 1. 0. 0. 0.
DEFAULTVALUE              : 0.
D50                       : 0.4e-3
<endsandclass>
```

```
<beginsandclass>
NAME                      : SAND_3
UNITS                     : %
DESCRIPTION               : SAND_3
DYNAMIC                   : 1
FILE_IN_TIME              : NONE
INITIALIZATION_METHOD     : LAYERS
LAYERS_VALUES             : 1. 0. 0. 0. 0.
DEFAULTVALUE              : 0.
D50                       : 0.6e-3
<endsandclass>
```
