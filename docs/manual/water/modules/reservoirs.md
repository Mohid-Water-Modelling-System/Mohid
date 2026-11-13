# Module Reservoirs

## Overview

Module Reservoirs allows to compute reservois dynamic that change river flows and concentrations. Module Reservoirs recieves flows and concentrations from Module Drainage Network and gives flows and concentrations to the same module.

## Main Processes

The reservoirs do not have a explicit geometry and are located in a given Drainage Network node (that will interact flows and concetrations with it.

### Hydrology

The reservoirs minimum parameter is maximum volume (keyword MAX\_VOLUME) and optionally minimum volume (keyword MIN\_VOLUME, default is 0.0). Below minimum volume the reservoir does not discharge and above maximum volume the reservoir discharges the remainder so that volume does not get higher than maximum.

#### Inflows

The inflows are imposed by Drainage Network at the nodes where the reservoir are located.

#### Discharges

Via Module Discharges any discharge can be imposed (positive e.g. WWTP, or negative e.g. water abstraction for irigation).

#### Surface Fluxes

Surface fluxes as precipitation and evaporation are not computed yet but routines already exist for that purpose

#### Outflows

The outflows can be computed by defining several options:

**Outflow Unmanaged Reservoir**

1) If the maximum volume is the only information available, the reservoir will not discharge until complete full and thereafter will discharge the volume above maximum divided by time step.

2) There is the option to define unmanaged reservoir as a weir (by keyword IS\_WEIR : 1) and the reservoir outflow is defined by the weir equation (also defined in ModuleDischarges). Then the crest elevation and width need to be defined.

**Outflow Managed Reservoirs**

1) If beside the maximum volume only an environmental flow (keyword MIN\_OUTFLOW) is defined, then above minimum volume the environmental flow will be discharged until it reaches the maximum volume and thereafter will discharge the volume above maximum divided by time step.

2) If an operation curve is given (keyword OPERATION\_TYPE), on the form of:

1. level - outflow - OPERATION\_TYPE : 1
2. level - outflow = percentage of inflow - OPERATION\_TYPE : 2
3. percentage of res volume - outflow - OPERATION\_TYPE : 3
4. percentage of res volume - outflow = percentage of inflow - OPERATION\_TYPE : 4
5. percentage of res volume - outflow = percentage of max outflow - OPERATION\_TYPE : 5

then the curve is read and the outflow computed based on actual reservoir volume or level. The outflow can be given directly from volume or level or the curve can define a percentage of inflow, where the outflow will be a percentage of the inflow at that instant or a percentage of maximum outflow (usually at maximum volume).

The curves are given by <<beginoperation>> <<endoperation>> block where the first column is level (elevation) or percentage of reservoir volume and the second column is outflow (m3/s) or percentage of inflow or percentage of max outflow (0-1) and between defined points the model linearly interpolates.

If at any given instant the volume/level of the reservoir is below the lowest point of curves, it will impose environmental flow if available or zero. If at any given instant the volume/level of the reservoir is above the highest point of curves, it will impose the last point outflow.

**Reminder:** if level is used an accumulated volumes curve is required (volume - elevation)) through a <<beginaccvolumecurve>> <<endaccvolumecurve>> block where first column is reservoir volume and second is elevation.

3) There is also the option to impose the outflow (e.g. from measured data) where this is done via Module Discharges but the discharge needs to be negative and a keyword in discharge present (IS\_OUTFLOW)

This management option will be given priority and will be used even that environmental flow or operation curves are defined.

### Water Properties

Concentrations are computed based on two available methods (keyword PROP\_COMPUTE\_METHOD):

1. Instant Mixing - PROP\_COMPUTE\_METHOD : 1
2. Mixing occurs at retention time - PROP\_COMPUTE\_METHOD : 2

The instant mixing method computes the new concentration as the mixing in all reservoir would occur in the given time step being a mass balance between the existing and the entering mass divided by total volume

The second method uses the first as the target concentration and computes the concentration difference and the reservoir retention time (inflow /reservoir volume).

If the retention time is lower than the time step then the new concentration is the target concentration defined by instant mixing method. If not, the new concentration will be the old plus a percentage of the difference. And the percentage is given by time step / retention time.

## Data File

### KEYWORDS

```
KEYWORD                           TYPE    DEFAULT
RESERVOIR_FILE                    char    ...       !File with reservoirs properties
DISCHARGES                        0/1     0         !Connect/Disconnect discharges
SURFACE_FLUXES                    0/1     0         !Connect/Disconnect surface fluxes
BOTTOM_FLUXES                     0/1     0         !Connect/Disconnect bottom fluxes
PROP_COMPUTE_METHOD               int     1         !1 instant mixing; 2 retention time full mixing

CONTINUOUS                        0/1     0         !is a continuation?

INITIAL_VOLUME_DEFAULT_METHOD     int     1         !1 start perentage full
START_PERCENTAGE_FULL             real    0.0       !read if INITIAL_VOLUME_DEFAULT_METHOD : 1. percentage of max volume at start

<beginreservoir>
ID                                 int    ...      !Reservoir ID
NAME                               char   ...      !Reservoir Name
DN_NODE_ID                         int    ...      !Reservoir Locaton in Drainage Network Node ID
COORDINATES                        real(2) ...     !Reservoir Location in coordinates
GRID_I                             int     ...     !Reservoir Location in grid coordinates I
GRID_J                             int     ...     !Reservoir Location in grid coordinates J
MIN_VOLUME                         real    0.0     !Minimum reservoir volume - below no outflow
MAX_VOLUME                         real    ...     !Maximum reervoir volume - above outflows remainder
OPERATION_TYPE                     int      4      !level - outflow                                                - OPERATION_TYPE : 1
                                                   !level - outflow = percentage of inflow                         - OPERATION_TYPE : 2
                                                   !percentage of res volume - outflow                             - OPERATION_TYPE : 3
                                                   !percentage of res volume - outflow = percentage of inflow      - OPERATION_TYPE : 4
                                                   !percentage of res volume - outflow = percentage of max outflow - OPERATION_TYPE : 5
<<beginoperation>>                                 !Block to define operation curves
<<endoperation>>

<<beginaccvolumecurve>>                            !Block to define volume accumlated curves
<<endaccvolumecurve>>

MIN_OUTFLOW                        real    0.0     !Enviromental flow
MAX_OUTFLOW                        real    inf     !Maximum outflow (limited by reservoir hudraulic structures and mandatory in case of OPERATION_TYPE : 5)
<endreservoir>
```

### Sample

Reservoir.dat

```
RESERVOIR_FILE            : ..\General Data\Digital Terrain\ReservoirsOperation.dat
DISCHARGES                : 0
SURFACE_FLUXES            : 0
BOTTOM_FLUXES             : 0
PROP_COMPUTE_METHOD       : 1  !1 instant mixing; 2 retention time full mixing

CONTINUOUS                : 0

INITIAL_VOLUME_DEFAULT_METHOD         : 1    !1 start perentage full
START_PERCENTAGE_FULL        : 0.0  !read if INITIAL_VOLUME_DEFAULT_METHOD : 1. percentage of max volume at start

OUTPUT_TIME               : 0 1800
TIME_SERIE                : 1
TIME_SERIE_LOCATION       : ..\General Data\TimeSeries\ReservoirsTimeSeriesLocation.dat
```

Reservoirs File

```
<beginreservoir>
ID                        : 1
NAME                      : Reservoir Unmanaged
DN_NODE_ID                : 4
COORDINATES               : 82 32
GRID_I                    : 2
GRID_J                    : 54
MAX_VOLUME                : 10800.0
<endreservoir>
```

```
<beginreservoir>
ID                        : 2
NAME                      : Reservoir Managed
DN_NODE_ID                : 4
COORDINATES               : 825 325
GRID_I                    : 4
GRID_J                    : 9
!MIN_VOLUME                : 0.0
MAX_VOLUME                : 10800.0
OPERATION_TYPE            : 4
<<beginoperation>>
0.0 0.0
0.1 0.1
0.2 0.3
0.3 0.5
0.5 0.8
0.8 1.0
1.0 1.0
<<endoperation>>
<<beginaccvolumecurve>>
0 120
100 130
500 140
2000 150
1000 155
10800 157
<<endaccvolumecurve>>
!MIN_OUTFLOW               : 0.25
!MAX_OUTFLOW              :  5.0
<endreservoir>
```
