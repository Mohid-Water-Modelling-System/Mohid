# Continuous

## Overview

[CONTINUOUS](continuous.md) is the keyword that allows to pursue calculus from a previous run. It relies on the *.fin* and *.fin5* files that are created at the end of each run, or at the time intervals determined by the keyword [RESTART\_FILE\_OUTPUT\_TIME](restart-file-output-time.md). In some modules, another keyword is used: OLD.

## Usage

Make sure the above keyword is placed correctly in the *hydrodynamic.dat* and the *turbulence.dat*. (Other files may be required for this. Also, in some cases, in some files this may not be required.)

- In the Mohid GUI create a new run within the run you wish to continue. This will create a new run with an identical configuration of the previous run.
- Edit the *model.dat* and change the time. Make sure the starting time matches the previous run ending time.

### Syntax

- Add the following keyword to the *hydrodynamic.dat* and the *turbulence.dat* files

```
CONTINUOUS : 1
```

- Add the following keyword to each property of the *waterproperties.dat* or *lagrangian.dat*:

```
OLD        : 1
```

- Make sure you deactivate the SLOWSTART and the RAMP options.
- If you're not using the GUI then add these lines to the nomfich.dat files

```
AIRW_INI                      : ../../Simulation/res/InterfaceWaterAir_1.fin
BOT_INI                       : ../../Simulation/res/InterfaceSedimentWater_1.fin
IN_CNDI                       : ../../Simulation/res/Hydrodynamic_1.fin
EUL_INI                       : ../../Simulation/res/WaterProperties_1.fin
TURB_INI                      : ../../Simulation/res/GOTM_1.fin
```

### Modules

Here's the list of modules affected by this keyword

- [Module Hydrodynamic](../../water/modules/hydrodynamic.md)
- [Module Turbulence](../../water/modules/turbulence.md)
- [Module WaterProperties](../../water/modules/waterproperties.md)
- [Module InterfaceSedimentWater](../../water/modules/interfacesedimentwater.md)
- [Module InterfaceWaterAir](../../water/modules/interfacewaterair.md)

## How to restart from a crashed run

- Get the last .fin saved date in the file Error\_and\_Messages\_x.log. It should be something like this:

```
WARNING    INTERNAL   Hydrodynamic restart file saved       1998: 1: 1: 0: 0: 0
```

- Create a new run inside the run that has crashed with the starting time exactly like the saved .fin file.
- Open the new run's nomfich\_x.dat file and make sure that the keywords IN\_CNDI, EUL\_INI, BOT\_INI, TURB\_INI point to the previous(crashed) run .fin files.
    -   IN\_CNDI
    -   EUL\_INI (if running with water properties)
    -   BOT\_INI (if running with interface sediment-water properties)
    -   TURB\_INI(if running with [GOTM turbulence module](../../water/modules/gotm.md))
- Run the new run

## Advanced options

In the *hydrodynamic.dat* file you must take care about some options when you want to start your simulations from a previous run, namely remove all features where ramp options are used.

- When using a ramp period for wind forcing change from **WIND  : 2** to **WIND  : 1** and remove keyword WIND\_SMOOTH\_PERIOD from the data file.
- When using a slow start in imposing tide, i.e., if in the previous run the SLOWSTART keyword is activated then remove it from the data file.
- When using a cold period in the Flather boundary condition, i.e. if in the previous run the FLATHER\_COLD\_PERIOD keyword is activated then remove it.
- When using baroclinic ramp period, i.e. if the RAMP option was on in the previous run (**RAMP  : 1**) set it to **RAMP  : 0**. Remove keyword INERTIAL\_PERIODS from the data file.

In the *assimilation.dat* input data file you must also remove any feature regarding ramping options. If in the previous run you're using keywords COLD\_RELAX\_PERIOD and COLD\_ORDER then remove it from the data file.

## See also

- [RESTART\_FILE\_OUTPUT\_TIME](restart-file-output-time.md)
- [RESTART\_FILE\_OVERWRITE](restart-file-overwrite.md)
