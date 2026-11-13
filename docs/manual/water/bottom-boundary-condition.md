# Bottom boundary condition

For each water property the user needs to define if it wants or not to impose sediment-water fluxes of particulate or dissolved matter. An example of a input file InterfaceWaterAir\_x.dat where bottom boundary condition is defined is presented below.

## WaterProperties\_x.dat

For each property is possible to impose a bottom boundary condition. By default no bottom boundary condition is imposed (BOTTOM\_FLUXES : 0).

```
    <beginproperty>
     BOTTOM_FLUXES     : 1
    <endproperty>
```

## InterfaceSedimentWater\_x.dat

Bottom absolute rugosity in meters

```
   <begin_rugosity>
   DEFAULTVALUE                 : 0.0025
   <end_rugosity>
```

Critical shear stress erosion in Pa

```
   <begin_critical_shear_erosion>
   INITIALIZATION_METHOD        : CONSTANT
   DEFAULTVALUE                 : 0.4
   <end_critical_shear_erosion>
```

Critical shear stress deposition in Pa

```
   <begin_critical_shear_deposition>
   INITIALIZATION_METHOD        : CONSTANT
   DEFAULTVALUE                 : 0.2
   <end_critical_shear_deposition>
```

Erosion rate in kg/m2/s

```
   <begin_erosion_rate>
   INITIALIZATION_METHOD        : CONSTANT
   DEFAULTVALUE                 : 0.00005
   <end_erosion_rate>
```

Definition of the bottom boundary of the particulate property cohesive sediment

```
   <beginproperty>
   NAME                         : cohesive sediment
   UNITS                        : kg/m2
   DESCRIPTION                  : cohesive sediment mass
   OLD                          : 1
   PARTICULATE                  : 1
   WATER_FLUXES                 : 1
   EROSION                      : 1
   DEPOSITION                   : 1
   MASS_LIMITATION              : 1
   MASS_MIN                     : 1e-8
   MIN_VALUE                    : 1e-4
   FILE_IN_TIME                 : NONE
   INITIALIZATION_METHOD        : ASCII_FILE
   FILENAME                     : ..\..\GeneralData\BottomSed\BottomSed.dat
   REMAIN_CONSTANT              : 0
   DEFAULTVALUE                 : 0.01
   TIME_SERIE                   : 0
   BOX_TIME_SERIE               : 0
   OUTPUT_HDF                   : 1
   <endproperty>
```
