# Point Origin

To define a point origin in the lagrangian module you need to define the coordinates where the tracers are to be released. This can be done in several ways by adding the keywords:

```
POSITION_COORDINATES    : -10.0 43
```

where the values are the X and Y geographical coordinates,

```
POSITION_CELLS          : 45.5 21.5
```

are the values are the I and J cells' indexes.

NOTE: The above is not an exhaustive list of ways for how to insert a point origin.

## Example Lagrangian\_1.dat

```
OUTPUT_TIME             : 0 3600
OUTPUT_MAX_TRACER       : 1
OUTPUT_CONC             : 2

DT_PARTIC               : 60

<BeginOrigin>
ORIGIN_NAME             : Off Vigo
GROUP_ID                : 1
EMISSION_SPATIAL        : Point
EMISSION_TEMPORAL       : Instantaneous
BOTTOM_EMISSION         : 0
OLD                     : 0

POINT_VOLUME            : 1
NBR_PARTIC              : 8

FLOAT                   : 1

MOVEMENT                : SullivanAllen
VARVELHX                : 0.1
VARVELH                 : 0.03
TURB_V                  : Constant
VARVELVX                : 0.01
VARVELV                 : 0.003

!longitude -10.0ºW, latitude 43.0ºN
POSITION_COORDINATES    : -10.0 43
!POSITION_CELLS          : 45.5 21.5
DEPTH_METERS            : 0.

START_PARTIC_EMIT       : 2007 04 25 09 33 20

<EndOrigin>
```
