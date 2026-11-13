# Nesting models

Here is a guide on configuring [MOHID Water](../water/mohid-water.md) with nested models.

## In MOHID GUI

- Create a new project in [MOHID GUI](../tools/gui/mohid-gui.md),
- Create a Model by inserting a new Simulation,
    -   Setup your model
- Create a Submodel by inserting a new Simulation from the current Simulation. The new simulation should fold into the previous one.
    -   Setup the model

## Setting up the files

The whole issue is how to define the open-boundary conditions for the submodel and the relaxation condition. Basically we suggest to use a Flather radiation condition at the open boundaries that will radiate the water level and the barotropic flux in conjunction with a flow relaxation scheme spanning a few cells from the open-boundaries to the interior for the velocities (*u* and *v*).

### Keywords in Hydrodynamic.dat

In the submodel hydrodynamic file

```
SUBMODEL          : 1

!Flather radiation scheme using the father model as a reference solution.
RADIATION         : 2
LOCAL_SOLUTION    : 5 ! Submodel + assimilation(relaxes the solution to the father(courser grid) domain
Other options for LOCAL_SOLUTION include:
1 - NoLocalSolution_,      2 - Submodel_, 3 - AssimilationField_, 4 - Gauge_,
5 - AssimilaPlusSubModel_, 6 - GaugePlusSubModel_, 7 -AssimilaGaugeSubModel_

!Activate the Flow relaxation scheme:
DATA_ASSIMILATION : 1
BRFORCE           : 1 ! relaxation of the baroclinic force
Only if these keywords are activated will the assimilation be read!.
```

### Keywords in Assimilation.dat

In the submodel assimilation file (assimilation.dat):

If the Keywords SUBMODEL, BRFORCE and DATA\_ASSIMILATION are active and LOCAL\_SOLUTION is 2 then the block <<begin\_field>> will not be read because the assimilation will be done using the father model results.

Example with LOCAL\_SOLUTION = 5

```
<beginproperty>
NAME                    : velocity U
UNITS                   : m/s
DIMENSION               : 3D
OUTPUT_HDF              : 1

COLD_RELAX_PERIOD       : 432000
COLD_ORDER              : 5

<<begin_field>>
DEFAULTVALUE            : 0
INITIALIZATION_METHOD   : HDF
FILE_IN_TIME            : HDF
FILENAME                : ../../GeneralData/Reference_vel_U.hdf5
TYPE_ZUV                : z
<<end_field>>

<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : u
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : ASCII_FILE
FILENAME                : ../../GeneralData/RelaxationCoefs_U.dat

<<end_coef>>

<endproperty>

<beginproperty>
NAME                    : velocity V
UNITS                   : m/s
DIMENSION               : 3D
OUTPUT_HDF              : 1

COLD_RELAX_PERIOD       : 432000
COLD_ORDER              : 5

<<begin_field>>
DEFAULTVALUE            : 0
INITIALIZATION_METHOD   : HDF
FILE_IN_TIME            : HDF
FILENAME                : ../../GeneralData/Reference_vel_V.hdf5
TYPE_ZUV                : z
<<end_field>>

<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : v
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : ASCII_FILE
FILENAME                : ../../GeneralData/RelaxationCoefs_V.dat

<<end_coef>>

<endproperty>
```

```
Example 2: LOCAL_SOLUTION = 2 (only relaxes to father domain)
```

```
<beginproperty>
NAME                    : barotropic velocity U
UNITS                   : m/s
DIMENSION               : 2D
```

```
OUTPUT_HDF              : 1
```

```
<<begin_field>>
DEFAULTVALUE            : 0  ! will not be read!
TYPE_ZUV                : z
<<end_field>>
```

```
<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : u
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : SPONGE
<<end_coef>>
<endproperty>
```

```
<beginproperty>
NAME                    : barotropic velocity V
UNITS                   : m/s
DIMENSION               : 2D

OUTPUT_HDF              : 1
```

```
<<begin_field>>
DEFAULTVALUE            : 0! will not be read!
TYPE_ZUV                : z
<<end_field>>
```

```
<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : u
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : SPONGE
<<end_coef>>
<endproperty>
```

```
<beginproperty>
NAME                    : water level
UNITS                   : m
DIMENSION               : 2D
OUTPUT_HDF              : 1

<<begin_field>>
DEFAULTVALUE            : 0! will not be read!
TYPE_ZUV                : z
<<end_field>>
```

```
<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : u
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : SPONGE
<<end_coef>>
<endproperty>
```

```
<beginproperty>
NAME                    : velocity U
UNITS                   : m/s
DIMENSION               : 3D
OUTPUT_HDF              : 1
```

```
!COLD_RELAX_PERIOD       : 43200.
!COLD_ORDER              : 5
```

```
<<begin_field>>
DEFAULTVALUE            : 0! will not be read!
TYPE_ZUV                : z
<<end_field>>
```

```
<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : u
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : SPONGE
<<end_coef>>
<endproperty>
```

```
<beginproperty>
NAME                    : velocity V
UNITS                   : m/s
DIMENSION               : 3D
OUTPUT_HDF              : 1
```

```
!COLD_RELAX_PERIOD       : 43200.
!COLD_ORDER              : 5
```

```
<<begin_field>>
DEFAULTVALUE            : 0! will not be read!
TYPE_ZUV                : z
<<end_field>>
```

```
<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : v
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : SPONGE
<<end_coef>>
<endproperty>
```

```
<beginproperty>
NAME                    : vertical z
UNITS                   : m
DIMENSION               : 3D
OUTPUT_HDF              : 1
```

```
!COLD_RELAX_PERIOD       : 43200.
!COLD_ORDER              : 5
```

```
<<begin_field>>
DEFAULTVALUE            : 0! will not be read!
TYPE_ZUV                : z
<<end_field>>
```

```
<<begin_coef>>
DEFAULTVALUE            : 1e9
TYPE_ZUV                : z
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 1
INITIALIZATION_METHOD   : SPONGE
<<end_coef>>
<endproperty>
```
