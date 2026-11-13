# Atmospheric boundary condition

For each water property the user needs to define if it want or not to impose atmospheric fluxes. Here's an example of atmospheric properties that need to be defined -if the user wants to compute heat fluxes with standard bulk formulas-. The atmospheric properties are defined in the Atmospheric\_x.dat input file. The heat (bulk formulas), mass and momentum fluxes are defined in the InterfaceWaterAir\_x.dat. The example of an input file InterfaceWaterAir\_x.dat is presented below.

## WaterProperties\_x.dat

For each property it's possible to impose an atmospheric (or surface) boundary condition. By default no atmospheric boundary condition is imposed (SURFACE\_FLUXES : 0).

```
    <beginproperty>
     SURFACE_FLUXES     : 1
    <endproperty>
```

## Atmospheric\_x.dat

It's important to have the output of the atmospheric properties fields so the user can confirm what is being imposed

```
   OUTPUT_TIME             : 0. 86400
```

Zonal wind velocity

```
   <beginproperty>
   NAME                    : wind velocity X
   UNITS                   : m/s
   DESCRIPTION             : wind velocity X
   DEFAULTVALUE            : 0.
   FILE_IN_TIME            : HDF
   FILENAME                : D:\Aplica\AguasAlgarve\DadosBase\MM5\AlgarveCoast\MM5Out_D3D2_July2004.hdf5
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Meridional wind velocity

```
   <beginproperty>
   NAME                    : wind velocity Y
   UNITS                   : m/s
   DESCRIPTION             : wind velocity Y
   DEFAULTVALUE            : 0.
   FILE_IN_TIME            : HDF
   FILENAME                : D:\Aplica\AguasAlgarve\DadosBase\MM5\AlgarveCoast\MM5Out_D3D2_July2004.hdf5
   OUTPUT_HDF              : 1
   <endproperty>
```

!!!!!!Since 2016 code, zonal and meridional wind velocity are defined in one single block:

```
   <beginproperty>
   NAME                    : wind velocity
   UNITS                   : m/s
   DESCRIPTION             : wind velocity
   DEFAULTVALUE            : -1.0 0.0
   FILE_IN_TIME            : HDF
   FILENAME_X              : ..\..\GeneralData\Atmosphere\MM5_20070606_Portugal.hdf5
   FILENAME_Y              : ..\..\GeneralData\Atmosphere\MM5_20070606_Portugal.hdf5
   HDF_FIELD_NAME_X        : wind velocity X
   HDF_FIELD_NAME_Y        : wind velocity Y
   REMAIN_CONSTANT         : 0
   TIME_SERIE              : 1
   OUTPUT_HDF              : 0
   <endproperty>
```

Atmospheric pressure

```
   <beginproperty>
   NAME                    : atmospheric pressure
   UNITS                   : Pa
   DESCRIPTION             : NOAA reanalysis atmospheric pressure at sea level
   DEFAULTVALUE            : 0.
   FILE_IN_TIME            : HDF
   FILENAME                : D:\Aplica\AguasAlgarve\DadosBase\MM5\AlgarveCoast\MM5Out_D3D2_July2004.hdf5
   OUTPUT_HDF              : 1
   <endproperty>
```

Solar radiation

```
   <beginproperty>
   NAME                    : solar radiation
   UNITS                   : W/m^2
   DESCRIPTION             : meteoIST Solar Radiation
   ALBEDO                  : 0.05  ![%]
   FILE_IN_TIME            : HDF
   FILENAME                : D:\Aplica\AguasAlgarve\DadosBase\MM5\AlgarveCoast\MM5Out_D3D2_July2004.hdf5
   DEFAULTVALUE            : 0.0
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Air temperature

```
   <beginproperty>
   NAME                    : air temperature
   UNITS                   : ºC
   DESCRIPTION             : Temperature
   DEFAULTVALUE            : 0.
   FILE_IN_TIME            : HDF
   FILENAME                : D:\Aplica\AguasAlgarve\DadosBase\MM5\AlgarveCoast\MM5Out_D3D2_July2004.hdf5
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Relative humidity

```
   <beginproperty>
   NAME                    : relative humidity
   UNITS                   : fraction
   DESCRIPTION             : Humidity
   DEFAULTVALUE            : 0.
   FILE_IN_TIME            : HDF
   FILENAME                : D:\Aplica\AguasAlgarve\DadosBase\MM5\AlgarveCoast\MM5Out_D3D2_July2004.hdf5
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Cloud cover

```
   <beginproperty>
   NAME                    : cloud cover
   UNITS                   : fraction
   DESCRIPTION             : cloud cover
   DEFAULTVALUE            : 0.
   FILE_IN_TIME            : NONE
   !FILENAME                : D:\Aplica\AguasAlgarve\DadosBase\MM5\AlgarveCoast\MM5Out_D3D2_July2004.hdf5
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

## InterfaceWaterAir\_x.dat

Zonal wind stress

```
   <beginproperty>
   NAME                    : wind stress X
   UNITS                   : N/m2
   DESCRIPTION             : calculated wind stress X
   DEFAULTVALUE            : 0
   REMAIN_CONSTANT         : 0
   TYPE_ZUV                : z
   OUTPUT_HDF              : 1
   TIME_SERIE              : 0
   <endproperty>
```

Meridional wind stress

```
   <beginproperty>
   NAME                    : wind stress Y
   UNITS                   : N/m2
   DESCRIPTION             : calculated wind stress Y
   DEFAULTVALUE            : 0
   TYPE_ZUV                : z
   REMAIN_CONSTANT         : 0
   OUTPUT_HDF              : 1
   TIME_SERIE              : 0
   <endproperty>
```

!!!!!!Since 2016 code, Zonal and Meridional wind stress are defined in one single block:

```
   <beginproperty>
   NAME                      : wind stress
   UNITS                     : N/m2
   DESCRIPTION               : calculated wind stress X and Y
   DEFAULTVALUE              : 0.0 0.0
   CDWIND                    : 0.0015
   REMAIN_CONSTANT           : 0
   OUTPUT_HDF                : 0
   TIME_SERIE                : 0
   <endproperty>
```

Latent heat

```
   <beginproperty>
   NAME                    : latent heat
   UNITS                   : W/m^2
   DESCRIPTION             : Calculated latent heat
   FILE_IN_TIME            : NONE
   REMAIN_CONSTANT         : 0
   DEFAULTVALUE            : 0.
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Sensible heat

```
   <beginproperty>
   NAME                    : sensible heat
   UNITS                   : W/m^2
   DESCRIPTION             : Calculated sensible heat
   FILE_IN_TIME            : NONE
   REMAIN_CONSTANT         : 0
   DEFAULTVALUE            : 0.
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Solar radiation at sea level

```
   <beginproperty>
   NAME                    : surface radiation
   UNITS                   : W/m^2
   DESCRIPTION             : Calculated infrared radiation
   ALBEDO                  : 0.05 ![%]
   FILE_IN_TIME            : NONE
   REMAIN_CONSTANT         : 0
   DEFAULTVALUE            : 0.
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Net long wave radiation = (downward - upward) long wave radiation

```
   <beginproperty>
   NAME                    : net long wave radiation
   UNITS                   : W/m^2
   DESCRIPTION             : Calculated net long wave radiation
   FILE_IN_TIME            : NONE
   REMAIN_CONSTANT         : 0
   DEFAULTVALUE            : 0.
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Downward long wave radiation

```
   <beginproperty>
   NAME                    : downward long wave radiation
   UNITS                   : W/m^2
   DESCRIPTION             : Calculated downward long wave radiation
   FILE_IN_TIME            : HDF
   FILENAME                : D:\Aplica\AguasAlgarve\DadosBase\MM5\AlgarveCoast\MM5Out_D3D2_July2004.hdf5
   DEFAULTVALUE            : 0.
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Upward long wave radiation

```
   <beginproperty>
   NAME                    : upward long wave radiation
   UNITS                   : W/m^2
   DESCRIPTION             : Calculated upward long wave radiation
   FILE_IN_TIME            : NONE
   REMAIN_CONSTANT         : 0
   DEFAULTVALUE            : 0.
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```

Non solar flux = Net long wave radiation + Latent + Sensible

```
   <beginproperty>
   NAME                    : non solar flux
   UNITS                   : W/m^2
   DESCRIPTION             : Calculated infrared radiation
   FILE_IN_TIME            : NONE
   REMAIN_CONSTANT         : 0
   DEFAULTVALUE            : 0.
   TIME_SERIE              : 0
   OUTPUT_HDF              : 1
   <endproperty>
```
