# Module Lagrangian

## Overview

Lagrangian transport models are very useful to simulate localized processes with sharp gradients (submarine outfalls, sediment erosion due to dredging works, hydrodynamic calibration, oil dispersion, etc.). MOHID’s Lagrangian module uses the concept of lagrangian tracers. The most important property of a tracer is its position (x,y,z). For a physicist a tracer can be a water mass, for a geologist it can be a sediment particle or a group of sediment particles and for a chemist it can be a molecule or a group of molecules. A biologist can spot phytoplankton cells in a tracer (at the bottom of the food chain) as well as a shark (at the top of the food chain), which means that a model of this kind can simulate a wide spectrum of processes. The movement of the tracers can be influenced by the velocity field from the hydrodynamic module, by the wind from the surface module, by the spreading velocity from oil dispersion module and by random velocity. At the present stage the model is able to simulate oil dispersion, water quality processes and sediment transport. To simulate oil dispersion the lagrangian module interacts with the [oil dispersion module](oil.md), to simulate water quality the lagrangian module uses the [water quality module](waterquality.md). Sediment transport can be associated directly to the tracers using the concept of [settling velocity](../../reference/settling-velocity.md). Another feature of the lagrangian transport model is the ability to calculate residence time. This can be very useful when studying the exchange of water masses in bays or estuaries.

## Concepts

Like referred above, the Lagrangian module uses the concept of tracer. The tracers are characterized by there spatial coordinates, volume and a list of properties (each with a given concentration). The properties can be the same one like the ones described in the water properties module or coliform bacteria. Each tracer has associated a time to perform the random movement. The tracers are “born” at origins. Tracers which belong to the same origin have the same list of properties and use the same parameters for random walk, coliform decay, etc. Origins can differ in the way they emit tracers.

There are three different ways to define origins in space:

- a **[Point Origin](../point-origin.md)** emits tracers at a given point;
- a **[Box Origin](../../io/formats/boxes.md)** emits tracers over a given area;
- a **Accident Origin** emit tracers in a circular form around a point;

There are two different ways in which origins can emit tracers in time:

- a **Continuous Origin** emits tracers during a period of time;
- a **Instantaneous Origin** emits tracers at one instant;
- a **Moving Origin** emits tracers during a period of time along a defined track;

Origins can be grouped together in Groups. Origins which belong to the same group are grouped together in the output file, so it is more easy to analyze the results.

## Main processes

### Tracer Movement

### Turbulent Diffusion

### Residence Time

#### [Monitor boxes](../../io/formats/boxes.md)

### Oil spills

### Water quality processes

To use the Water Quality processes in the lagrangian module it is necessary to define a WaterQuality.dat file different from the one used for the eulerian calculations. Defining its location by the keyword WQM\_DATA\_FILE.

The waterquality files need to include the different parameters separately in two files in the case that both waterquality processes are calculated.

#### Larvae

To model larvae as lagrangian properties it is necessary to configure the properties larvae, age and oxygen in the lagrangian file

```
<BeginOrigin>
ORIGIN_NAME           : Generic
EMISSION_SPATIAL      : Box
EMISSION_TEMPORAL     : Instantaneous
MOVEMENT              : NotRandom
ADVECTION             : 1
BOX_NUMBER            : 1
BOXVOLINIC            : 10000
OLD                   : 0
WQM_DATA_FILE         : ..\..\Test\data\WaterQuality_3.dat

<<BeginProperty>>
NAME                  : larvae
UNITS                 : kg/m3
CONCENTRATION         : 10e4
MIN_CONCENTRATION     : 0.0
AMBIENT_CONC          : 0.0
<<EndProperty>>

<<BeginProperty>>
NAME                  : age
UNITS                 : ---
CONCENTRATION         : 5.6
MIN_CONCENTRATION     : 0.
AMBIENT_CONC          : 0.0
<<EndProperty>>

<<BeginProperty>>
NAME                  : oxygen
UNITS                 : mg/l
CONCENTRATION         : 10.
MIN_CONCENTRATION     : 0.0
AMBIENT_CONC          : 10.0
<<EndProperty>>

<EndOrigin>
```

#### Ecology

#### First order decayment (coliform bacteria)

### Sediment transport

### Contaminants/Partition coefficients

## Other features

## User manual

### Setup

Add in **nomfich.dat**

```
PARTIC_DATA                   : ../../WestIberia_ET2_2K4/Portugal/data/Lagrangian_1.dat
PARTIC_HDF                    : ../../WestIberia_ET2_2K4/Portugal/res/Lagrangian_1.hdf
PARTIC_FIN                    : ../../WestIberia_ET2_2K4/Portugal/res/Lagrangian_1.fin
```

Add in **model.dat**

```
LAGRANGIAN                    : 1
```

Create **Lagrangian.dat**

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

POSITION_COORDINATES    : -10.0 43
!POSITION_CELLS          : 45.5 21.5
DEPTH_METERS            : 0.

START_PARTIC_EMIT       : 2007 04 25 09 33 20

<EndOrigin>
```

## References

- [MOHID keyword Database](../../io/keywords/keyword-list.md)
