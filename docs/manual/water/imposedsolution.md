# ImposedSolution

The idea behind the imposed solution, is the ability to load pre-calculated hydrodynamic and waterproperties values from hdf5 files. Typically is used to run lagrangian tracers in off-line mode, to run submodels from [hydrodynamic](modules/hydrodynamic.md) and [waterproperties](modules/waterproperties.md) hdf5 output windows, to calculate [residual fields](../reference/residual.md) of velocity and water level, as well as the [energy](energy.md) timeserie.

## Sample configuration

The geometry\_X.dat should correspond to the one found in the HDF5 file, and in case that the window file does not include a bottom layer then a keyword should be included in the file indicating that is a window file

```
WINDOW                    : 0/1   (Default 0, if 1 does not check the consistency of the bottom layer)
```

In addition, the hydrodynamic, interface sediment water and water properties data file should be configured following the inditacions below. A typical file for the hydrodyamic\_x.dat would be similar to:

```
EVOLUTION                 : Imposed Solution

<begin_waterlevel>
NAME                      : water level
UNITS                     : m
DIMENSION                 : 2D
DEFAULTVALUE              : 0
INITIALIZATION_METHOD     : hdf
FILE_IN_TIME              : hdf
FILENAME                  : ...        (location of the hdf5 file to be imposed)
<end_waterlevel>

<begin_velocity_u>
NAME                      : velocity U
UNITS                     : m/s
DIMENSION                 : 3D
DEFAULTVALUE              : 0
INITIALIZATION_METHOD     : hdf
FILE_IN_TIME              : hdf
FILENAME                  : ...        (location of the hdf5 file to be imposed)
<end_velocity_u>

<begin_velocity_v>
NAME                      : velocity V
UNITS                     : m/s
DIMENSION                 : 3D
DEFAULTVALUE              : 0
INITIALIZATION_METHOD     : hdf
FILE_IN_TIME              : hdf
FILENAME                  : ...        (location of the hdf5 file to be imposed)
<end_velocity_v>
```

The InterfaceSedimentWater\_X.dat should be configured as:

```
<begin_rugosity>
DEFAULTVALUE              : 0.0
<end_rugosity>
```

The Waterproperties\_X.dat file has to be configured using the following indications for each property:

```
<beginproperty>
NAME                           : ...        (Water property name)
UNITS                          : ...        (Water property units)
DEFAULTVALUE                   : ...        (Water property default value. Not used but mandatory)
ADVECTION_DIFFUSION            : 0
DIMENSION                      : 2D/3D      (Dimension of the property)
INITIALIZATION_METHOD          : hdf
FILE_IN_TIME                   : hdf
FILENAME                       : ...        (location of the hdf5 file to be imposed)
<endproperty>

<beginproperty>
...
...
<endproperty>
```

## See also

- [Evolution](evolution.md)
