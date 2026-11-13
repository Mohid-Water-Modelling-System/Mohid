# Module LagrangianGlobal

## Overview

The LangrangianGlobal Module is a deep upgrade of the [Lagrangian Module](lagrangian.md). These upgrade has three main goals:

1 - Run only one lagrangian model in a multi-nesting implementation. This way particles can go from grid to grid without any problems (see example of particles moving between grids). In the Lagrangian module particles are destroyed when they leave the grid where they were emitted. This is the reason way a lagrangian model is runned for each grid in the Lagrangian Module;

2 - Centralise all the interpolation and particle location (in a [grid](../../io/formats/grid.md)) methods in the [HorizontalGrid Module](README.md#module-horizontalgrid#module-horizontalgrid) and (see source code [

$$
1
$$

](http://mohid.codeplex.com/SourceControl/changeset/view/64163#1208106));

3 - Avoid redundances in the particle origin definitions. This way is possible to compress the input file (Lagrangian\_x.dat).

**If the mohid user wants to activate this module it needs in the mohid compilation process to predefined a preprocessor symbol called \_LAGRANGIAN\_GLOBAL\_.**

## Related projects

This development was made by [Hidromod](http://www.Hidromod.com) in the framework of the follow projects:

- [EASY](http://www.maretec.mohid.com/projects/easy/webgis/easy/easymap.aspx) - Interreg project - Hidromod subcontract by [IST](http://www.utl.ist.pt) - Debry tracjetory forecast in ocean/coastal domains

- [CoWaMa](http://www.aquaplansistemas.es/es/soluciones/hidrometeorologia/cowama-gestion-calidad-aguas-de-ba%C3%B1o) Hidromod subcontract by [Clabsa](http://www.Clabsa.es) and [Lyonnaise-des-Eaux](http://www.lyonnaise-des-eaux.fr/) - Forecast of fecal contamination in beaches

- [Lenvis](http://www.hidromod.com/Lenvis/LenvisMohidOnline/LENVIS/MohidOnline/) - FP7 - Hidromod [consortium member](http://www.lenvis.eu/Partners.aspx) - Web model service focus in simulating the impact of emergency fecal contamination discharges;

- [Argomarine](http://www.argomarine.eu/public) - FP7 project - Hidromod subcontract by [Algarve University](http://www.ualg.pt) - FP7 - Implementation of a oil spill forecast model in the Toscanny region (Italy)

- [EASYCO](http://www.project-easy.info) - Interreg project - Hidromod subcontract by [IST](http://www.utl.ist.pt) - Data analysis and oil spill modelling via web

## Main upgrades

### Goal 1 - Run only one Lagrangian Model

The LagrangianGlobal module is a module very similar to the [Lagrangian module](lagrangian.md). The input data file keywords are exactly the same. The modules give the same lagrangian result if the user only runs one model (no nesting). In the case of the first module when are run several nesting levels. The user only define one lagrangian input data file in the data file of the first nesting level (first model in the tree.dat). Each lagrangian tracer will use the hydrodynamic field with the higher priority depending on the tracer position. By default the priority is define inverting the tree.dat order. The user can specify the oder using a block define in the lagrangian input file where the user can define the model priority:

```
<BeginModelPriority>
Model name x
Model name y
<EndModelPriority>
```

The model are order by decreasing priority.

### Goal 2 - Centralize all interpolation and particle location methods

In the [Lagrangian module](lagrangian.md) the interpolation and particle location methods were developed in a independent way of similar methods implemented in the [HorizontalGrid Module](README.md#module-horizontalgrid#module-horizontalgrid) for the [hydrodynamic](hydrodynamic.md) and [water properties](waterproperties.md) nesting. Now the lagrangian mdoule uses the same interpolation and location methods as the Eulerian modules. The subroutines from the HorizontalGrid Module used in the LagrangianGlobal Module are:

- GetXYInsideDomain

Check if a location X,Y (eg. particle location) is inside a specific model domain;

- GetXYCellZ

For a location X,Y (eg. particle location) and a specific model domain returns the correspondent cell (i,j) and even the relative position of the location inside the cell;

- GetCellZ\_XY

Do the inverse process ofsuborutine GetXYCellZ

- InterpolXYPoint

Interpolate a generic property for a location X,Y.

All this subroutines work in all [grid](../../io/formats/grid.md) types support by Mohid (uniform variable and constant spatial step, uniform rotate and curvilinear).

For the curvilinear case was necessary to develop a quite complex method to [calculate the relative position of a particle in a cell](../../reference/relative-position-of-a-particle-in-a-cell.md).

### Goal 3 - Decrease redundances in Lagrangian\_x.dat

To decrease this redundance the "clone origin" concept was created. There are several [Lagrangian Origins](lagrangian.md) options.

## Ilustrative examples of new options

### Clone Origin concept

### Fecal contamination variable discharges

The keywords used to controle the discretization (number of particles) in a operational system focus in the forecast of water quality based in Lagrangian discharges are the follow:

The maximum volume allowed for each particle (in this case 1 m3)

```
POINT_VOLUME                   : 1
```

Time between particle emission in a origin (in this case 100 s)

```
DT_EMIT                        : 100
```

Time the particle takes to double its volume due to turbulence(in this case 7200 s or 2 h)

```
TVOL200                        : 7200.
```

Factor that multiply by the particle initial volume gives the volume limit above which the particle is killed.

```
VOLFAC                         : 100.
```

True (or 1) if the particles volume is compute automaticly function of a flow time serie

```
ESTIMATE_MIN_VOL               : 1
```

In this case at least for each run for each origin at least 1000 particles will be emitted along the run.

```
MAX_PART                       : 1000
```

#### Particle volume computation procedure per Origin

**Step 1**

If the option ESTIMATE\_MIN\_VOL is ON (ESTIMATE\_MIN\_VOL : 1) a first estimation of the particles volumes for each run (eg. 1 day simulation) and for each origin is computed:

Particle Volume estimated = max (Total Discharge Volume in origin X along the run / MAX\_PART, POINT\_VOLUME)

**Step 2**

After for each emission iteration is computed the volume to be discharge. The period of each emission iteration is from Now to Now + DT\_EMIT

Along this period is compute the total volume emitted in a specific emission iteration (Total Volume by Iteration). This volume is computed based in the flow time serie associated to the origin.

**Step 3**

For each emission iteration is compute the actual particle volume emitted

**Step 3.1**

Number of particles emitted per iteration and per origin = int (Total Volume per Iteration and per Origin / Particle Volume Estimated) + 1

**Step 3.2**

Actual Particle Volume Emitted = Total Volume per Iteration and per Origin / Number of particles emitted per iteration and per origin

The follow example used as default in the [CoWaMa](http://www.aquaplansistemas.es/es/soluciones/hidrometeorologia/cowama-gestion-calidad-aguas-de-ba%C3%B1o) system.

```
<BeginGeneralKeywords>
 OUTPUT_TIME                    : 0 3600
```

```
 RESTART_FILE_OUTPUT_TIME       : 0 3600
 RESTART_FILE_OVERWRITE         : 0
```

```
 DT_PARTIC                      : 100
```

```
 PARTIC_BOX                     : ..\..\LEVEL1_VALIDATIONECOLI\..\Lagrangian\Level2Validation\data\BOXRESTARTFILE_Lagrangian.dat
```

```
 TIME_SERIE_LOCATION            : ..\..\LEVEL1_VALIDATIONECOLI\..\Lagrangian\Level2Validation\data\TimeSeriesLocationLag.dat
```

```
 OUTPUT_CONC                    : 2
```

```
 LW_PERCENTAGE                  : .4
 LW_EXTINCTION_COEF             : .3333333
 SW_PERCENTAGE                  : .6
 SW_EXTINCTION_TYPE             : 3
```

```
 IGNORE_ON                      : 1
<EndGeneralKeywords>
```

```
<BeginModelPriority>
Level3ValidationEcoli
Level2Validation
<EndModelPriority>
```

```
<BeginOrigin>
GROUP_ID                       : 1
ORIGIN_NAME                    : Default
DEFAULT                        : 1
```

```
EMISSION_SPATIAL               : Point
EMISSION_TEMPORAL              : Continuous
EMISSION_ON                    : 1
OLD                            : 0
POSITION_COORDINATES           : 0. 0.
DEPTH_METERS                   : 0.
MOVEMENT                       : SullivanAllen
VARVELHX                       : 0.2
VARVELH                        : 0.1
TURB_V                         : Constant
VARVELVX                       : 0.02
VARVELV                        : 0.01
```

```
NBR_PARTIC                     : 1
```

The maximum volume allowed for each particle (in this case 1 m3)

```
POINT_VOLUME                   : 1
```

Time between particle emission in a origin (in this case 100 s)

```
DT_EMIT                        : 100
```

Time the particle takes to double its volume due to turbulence(in this case 7200 s or 2 h)

```
TVOL200                        : 7200.
```

Factor that multiply by the particle initial volume gives the volume limit above which the particle is killed.

```
VOLFAC                         : 100.
```

True (or 1) if the particles volume is compute automaticly function of a flow time serie

```
ESTIMATE_MIN_VOL               : 1
```

In this case at least for each run for each origin at least 1000 particles will be emitted along the run.

```
MAX_PART                       : 1000
```

```
FLOW                           : -9e16
FLOW_VARIABLE                  : 1
FLOW_COLUMN                    : 2
DISCHARGE_FILE                 : *****.***
!    T90_VAR_METHOD =  Canteras             = 1
!    T90_VAR_METHOD =  Chapra               = 2
!    T90_VAR_METHOD =  TimeSerie            = 3
```

```
<<BeginProperty>>
NAME                           : fecal coliforms
UNITS                          : MPN/100ml
CONCENTRATION                  : 0
CONC_VARIABLE                  : 1
CONC_COLUMN                    : 2
DISCHARGE_FILE                 : *****.***
T90_VARIABLE                   : 0
T90_VAR_METHOD                 : 0
T90                            : 10800
T90_FILE                       : T902.dat
T90_COLUMN                     : 2
OUTPUT_HDF                     : 1
TIME_SERIE                     : 0
<<EndProperty>>
<<BeginProperty>>
NAME                           : escherichia coli
UNITS                          : MPN/100ml
CONCENTRATION                  : 0
CONC_VARIABLE                  : 1
CONC_COLUMN                    : 2
DISCHARGE_FILE                 : *****.***
T90_VARIABLE                   : 0
T90_VAR_METHOD                 : 3
T90                            : 10800
T90_FILE                       : ..\..\LEVEL1_VALIDATIONECOLI\..\Lagrangian\Level2Validation\data\DailyT90.srw
T90_COLUMN                     : 2
OUTPUT_HDF                     : 1
TIME_SERIE                     : 0
<<EndProperty>>
<<BeginProperty>>
NAME                           : temperature
UNITS                          : ?C Summer temperature
CONCENTRATION                  : 0
OUTPUT_HDF                     : 1
TIME_SERIE                     : 0
<<EndProperty>>
<<BeginProperty>>
NAME                           : salinity
UNITS                          : ppm
CONCENTRATION                  : 0
OUTPUT_HDF                     : 1
TIME_SERIE                     : 0
<<EndProperty>>
<<BeginProperty>>
NAME                           : cohesive sediment
UNITS                          : mg/l
CONCENTRATION                  : 0
OUTPUT_HDF                     : 1
TIME_SERIE                     : 0
<<EndProperty>>
<EndOrigin>
```

```
<BeginOrigin_Clone>
CLONE                          : Default
DEFAULT                        : 0
GROUP_ID                       : 1
ORIGIN_NAME                    : Discharge 1
EMISSION_SPATIAL               : Point
EMISSION_TEMPORAL              : Continuous
DISCHARGE_FILE                 : ..\..\LEVEL1_VALIDATIONECOLI\..\Lagrangian\Level2Validation\data\DISCHARGE_Discharge 1.dat
POSITION_COORDINATES           : 1.3 41.11
OLD                            : 0
DEPTH_METERS                   : 1000000
<<BeginProperty>>
NAME                           : salinity
UNITS                          : PSU
CONCENTRATION                  : 0.2
<<EndProperty>>
<<BeginProperty>>
NAME                           : temperature
UNITS                          : ºC
CONCENTRATION                  : 19.5
<<EndProperty>>
<<BeginProperty>>
NAME                           : cohesive sediment
UNITS                          : mg/l
CONCENTRATION                  : 100
<<EndProperty>>
<<BeginProperty>>
NAME                           : fecal coliforms
UNITS                          : MPN/100ml
CONCENTRATION                  : 100000
DISCHARGE_FILE                 : ..\..\LEVEL1_VALIDATIONECOLI\..\Lagrangian\Level2Validation\data\DISCHARGE_FECAL_COLIFORMSDischarge 1.dat
<<EndProperty>>
<<BeginProperty>>
NAME                           : escherichia coli
UNITS                          : MPN/100ml
CONCENTRATION                  : 100000
DISCHARGE_FILE                 : ..\..\LEVEL1_VALIDATIONECOLI\..\Lagrangian\Level2Validation\data\DISCHARGE_E_COLIDischarge 1.dat
<<EndProperty>>
<EndOrigin_Clone>
```
