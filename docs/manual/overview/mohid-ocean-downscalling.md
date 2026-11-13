# Mohid Ocean Downscalling

## Download

The area extracted from NOAA/GFS and MyOcean/RTOFS must be bigger than that the area of MOHID level domain choosen to impose the surface and ocean boundary conditions.

### NOOA/GFS

To download GFS model solution go to [https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs). Then select your product of interest.

### CMEMS

The Copernicus Marine Environment Monitoring Service (CMEMS; [http://marine.copernicus.eu/](http://marine.copernicus.eu/)) is the follow up of the my ocean project. To download this solution access the catalogue [http://marine.copernicus.eu/services-portfolio/access-to-products/](http://marine.copernicus.eu/services-portfolio/access-to-products/) and select the GLOBAL\_ANALYSIS\_FORECAST\_PHY\_001\_024. You will need to register in order to download your product of interest.

### RTOFS

To download te RTOFS global model solution go to [http://polar.ncep.noaa.gov/global/data\_access.shtml](http://polar.ncep.noaa.gov/global/data_access.shtml) or [http://nomads.ncep.noaa.gov:9090/dods/rtofs](http://nomads.ncep.noaa.gov:9090/dods/rtofs). The download in the second one is made via opendap. select the product and thenclick in download data.

## Conversion

### Bathymetry transition

Running tool SmoothBathymNesting.exe

Example of the input file SmoothBathymNesting.dat

```
!File  bathymetry of the external 3D solution
 FATHER_BATIM                : Father.dat
```

```
!Bathymetry File of MOHID solution
 SON_BATIM                   : Son.dat
```

```
!New Bathymetry File
 NEW_SON_BATIM               : NewSon.dat
```

```
 <begin_coef>
```

```
 !Name of generic property
 NAME                        : generic property
```

```
 !Type of initialization used
 INITIALIZATION_METHOD       : sponge

 !0-external 3D solution, 1-MOHID solution
 DEFAULTVALUE                : 1

 !sponge output
 SPONGE_OUT                  : 0
```

```
 !sponge cells number
 SPONGE_CELLS                : 10

  !1-exponential, 2-linear
 SPONGE_EVOLUTION            : 2
```

```
<end_coef>
```

### Grib-Netcdf-Hdf5

#### NOAA/GFS

STEP 1: Grib to NetCDF

To convert grib files to NetCDF files create a batch file and write the follow comand line "java -Xmx1024m -classpath netcdfAll-4.3.jar ucar.nc2.dataset.NetcdfDataset -in gfs\_4\_20130415\_0000\_003.grb2 -out gfs\_4\_20130415\_0000\_003.nc > NCDumpLog.txt". The input file must be in the same folder of batch file.

STEP 2: NetCDF to HDF5

Running tool [ConversionHDF5.exe](../tools/converttohdf5.md)

Example of the file ConvertToHDF5.dat (see [keyword definition for converting from netcdf CF to mohid hdf5](../tools/converttohdf5.md))

```
<begin_file>
```

```
   ACTION              : CONVERT NETCDF CF TO HDF5 MOHID
   HDF5_OUT            : 1
   OUTPUTFILENAME      : outfile.hdf5
   NETCDF_OUT          : 1
   OUTPUT_NETCDF_FILE  : outfile.nc
```

```
   WINDOW_OUT : 1173 1323 1428 1643
   !is optional but usefull to generate smaller area for interpolation. In case of already extracted a small area durin the download data then change to WINDOW_OUT : 0.
```

```
 <<begin_time>>
   NETCDF_NAME         : time
 <<end_time>>
```

```
 <<begin_grid>>
   NETCDF_NAME_LAT     : lat
   NETCDF_NAME_LONG    : lon
   NETCDF_NAME_MAPPING : Temperature_height_above_ground
   MAPPING_LIMIT       : -10000
 <<end_grid>>
```

```
   PROPERTIES_NUMBER   : 8
```

```
 <<begin_field>>
   NETCDF_NAME         : u-component_of_wind_height_above_ground
   NAME                : wind velocity X
   UNITS               : m/s
   DESCRIPTION         : MOHID
   DIM                 : 2
 <<end_field>>
```

```
 <<begin_field>>
   NETCDF_NAME         : v-component_of_wind_height_above_ground
   NAME                : wind velocity Y
   UNITS               : m/s
   DESCRIPTION         : MOHID
   DIM                 : 2
 <<end_field>>
```

```
 <<begin_field>>
   NETCDF_NAME         : wind_modulus
   NAME                : wind modulus
   UNITS               : m/s
   DESCRIPTION         : MOHID
   DIM                 : 2
   VECTOR_INTENSITY    : 1
   VECTOR_X            : wind velocity X
   VECTOR_Y            : wind velocity Y
 <<end_field>>
```

```
 <<begin_field>>
   NETCDF_NAME         : Temperature_height_above_ground
   NAME                : air temperature
   UNITS               : oC
   DESCRIPTION         : MOHID
   DIM                 : 2
   ADD_FACTOR          : -273
 <<end_field>>
```

```
 <<begin_field>>
   NETCDF_NAME         : Total_precipitation_surface_3_Hour_Accumulation
   NAME                : precipitation
   UNITS               : mm/h
   DESCRIPTION         : MOHID
   DIM                 : 2
 <<end_field>>
```

```
 <<begin_field>>
   NETCDF_NAME         : Relative_humidity_height_above_ground
   NAME                : relative humidity
   UNITS               : -
   DESCRIPTION         : MOHID
   DIM                 : 2
   MULTIPLY_FACTOR     : 0.01
 <<end_field>>
```

```
 <<begin_field>>
   NETCDF_NAME         : Downward_Short-Wave_Radiation_Flux_surface_3_Hour_Average
   NAME                : solar radiation
   UNITS               : W/m^2
   DESCRIPTION         : MOHID
   DIM                 : 2
 <<end_field>>
```

```
 <<begin_field>>
   NETCDF_NAME         : Pressure_reduced_to_MSL_msl
   NAME                : atmospheric pressure
   UNITS               : pa
   DESCRIPTION         : Malaca
   DIM                 : 2
 <<end_field>>
```

```
 <<begin_input_files>>
   !path to the input files
 <<end_input_files>>
```

```
<end_file>
```

#### MyOcean and RTOFS

STEP 1: Grib to NetCDF

In case of download grib data use the same sample than that the one prescribed in GFS. If the data downloaded is already in NETCDF than go to Step 2.

STEP 2: NetCDF to HDF5

Running tool [ConversionHDF5.exe](../tools/converttohdf5.md)

Example of the file ConvertToHDF5.dat

```
<begin_file>
```

```
    ACTION              : CONVERT NETCDF CF TO HDF5 MOHID
    HDF5_OUT            : 1
    OUTPUTFILENAME      : 20130416_rtofs_glo_3dz_nowcast_daily.hdf5
    NETCDF_OUT          : 1
    OUTPUT_NETCDF_FILE  : rtofs_glo_3dz_nowcast_daily.nc
```

```
    WINDOW_OUT : 1173 1323 1428 1643
   !is optional but usefull to generate smaller area for interpolation. In case of already extracted a small area durin the download data then change to WINDOW_OUT : 0.
```

```
  <<begin_time>>
    NETCDF_NAME : time
  <<end_time>>
```

```
  <<begin_grid>>
    NETCDF_NAME_LAT     : lat
    NETCDF_NAME_LONG    : lon
    STARTS_180W         : 0
    NETCDF_NAME_MAPPING : temperature
    MAPPING_LIMIT       : 1.2676506E29
    MAPPING_INSTANT     : 2
    NETCDF_NAME_DEPTH   : lev
    INVERT_LAYER_ORDER  : 1
    BATHYM_FROM_MAP     : 1
    BATHYM_FILENAME     : Batim_rtofs_glo_3dz_nowcast_daily.dat
  <<end_grid>>
```

```
    PROPERTIES_NUMBER   : 5
```

```
  <<begin_field>>
    NETCDF_NAME         : temperature
    NAME                : temperature
    UNITS               : ºC
    DESCRIPTION         : MOHID
    DIM                 : 3
  <<end_field>>
```

```
  <<begin_field>>
    NETCDF_NAME         : salinity
    NAME                : salinity
    UNITS               : psu
    DESCRIPTION         : MOHID
    DIM                 : 3
  <<end_field>>

```

```
  <<begin_field>>
    NETCDF_NAME         : u
    NAME                : velocity U
    UNITS               : m/s
    DESCRIPTION         : MOHID
    DIM                 : 3
  <<end_field>>
```

```
  <<begin_field>>
    NETCDF_NAME         : v
    NAME                : velocity V
    UNITS               : m/s
    DESCRIPTION         : MOHID
    DIM                 : 3
  <<end_field>>
```

```
  <<begin_field>>
    NETCDF_NAME : velocity_modulus
    NAME                : velocity modulus
    UNITS               : m/s
    DESCRIPTION         : MOHID
    DIM                 : 3
    VECTOR_INTENSITY    : 1
    VECTOR_X            : velocity U
    VECTOR_Y            : velocity V
  <<end_field>>
```

```
  <<begin_input_files>>
    !path to the input files
  <<end_input_files>>
```

```
<end_file>
```

## Interpolation

Vertical discretization of RTOFS and MyOcean model solutions correspond to the depth of the center of the cell. Since MOHID [geometry](../base/module-geometry.md) must be provided in thickness the depths must be convertes to top of the faces cell and then to thickness layers. In order to avoid this the new method of interpolating grids allow to provide the depth of the center of the cells in ConvertToHDF5.dat instead to use a thickness calculated based on the centers depth. Th both samples are provided bellow.

### GFS to MOHID

This example perform linear horizontally interpolation (e.g., 2D) between GFS grid data and MOHID.

Running tool [ConversionHDF5.exe](../tools/converttohdf5.md)

Sample of the file ConvertToHDF5.dat

```
<begin_file>
```

```
  ACTION                    : INTERPOLATE GRIDS
```

```
  TYPE_OF_INTERPOLATION     : 1
  !2D horizontal interpolation 1-bilinear and 3-Triangulation
```

```
  INTERPOLATION3D           : 0
  !No interpolation 3D in vertical
```

```
  FATHER_FILENAME           : ..\outdata\GFS_2011010100_2011011621.hdf5
  FATHER_GRID_FILENAME      : batim\GFS.dat
```

```
  OUTPUTFILENAME            : ..\outdata\Batim_Caribe_Colombia_5km.hdf5
  NEW_GRID_FILENAME         : batim\Batim_Caribe_Colombia_5km_GFS.dat
```

```
  EXTRAPOLATE_2D            : 4
  EXTRAPOLATE_LIMIT         : -10000
```

```
<end_file>
```

### RTOFS to MOHID

This example perform 2D horizontal and 3D vertical interpolation between RTOFS grid data and MOHID. In this sample is provided the depth of the centers cells in the file.dat.

Running tool [ConversionHDF5.exe](../tools/converttohdf5.md)

Example of the file ConvertToHDF5.dat

```
<begin_file>
```

```
  ACTION                    : INTERPOLATE GRIDS
```

```
  DO_NOT_BELIEVE_MAP        : 1
```

```
  TYPE_OF_INTERPOLATION     : 1
  !2D horizontal interpolation 1-bilinear and 3-Triangulation
```

```
  INTERPOLATION3D           : 1
  !3D vertical interpolation
```

```
  START                     : 2013 04 13 0 0 0
  END                       : 2013 04 18 0 0 0
```

```
  FATHER_FILENAME           : ..\work\20130410_rtofs_glo_3dz_nowcast_daily.hdf5
  FATHER_GRID_FILENAME      : ..\work\batim\Batim_RTOFS.dat
```

```
  OUTPUTFILENAME            : ..\work\Level2_RTOFS_13_18_04_2013.hdf5
  NEW_GRID_FILENAME         : ..\work\batim\Level2_.new
```

```
  EXTRAPOLATE_2D            : 4
  EXTRAPOLATE_LIMIT         : -10000
```

```
  FATHER_GEOMETRY           : ..\work\Geometry_RTOFS.dat
  !A "false" geometry must be provide, but in reality uses the depths in the file.dat
```

```
  NEW_GEOMETRY              : ..\work\Geometry_1.dat
  !MOHID geometry with thicknes layers
```

```
 <<BeginDepths>>
  0.0
  10.0
  20.0
  30.0
  50.0
  75.0
  100.0
  125.0
  150.0
  200.0
  250.0
  300.0
  400.0
  500.0
  600.0
  700.0
  800.0
  900.0
 1000.0
 1100.0
 1200.0
 1300.0
 1400.0
 1500.0
 1750.0
 2000.0
 2500.0
 3000.0
 3500.0
 4000.0
 4500.0
 5000.0
 5500.0
 <<EndDepths>>

<end_file>
```

### MyOcean to MOHID

This example perform 2D horizontal and 3D vertical interpolation between MyOcean grid data and MOHID.

Running tool [ConversionHDF5.exe](../tools/converttohdf5.md)

Example of the file ConvertToHDF5.dat

```
<begin_file>
```

```
  ACTION                    : INTERPOLATE GRIDS
```

```
  START                     : 2013 04 12 12 0 0
  END                       : 2013 04 27 12 0 0
```

```
  TYPE_OF_INTERPOLATION     : 1
  !2D horizontal interpolation 1-bilinear and 3-Triangulation
```

```
  FATHER_FILENAME           :  ..\conv\MyOcean3D.hdf5
  FATHER_GRID_FILENAME      :  ..\conv\MyOceanMaxDepth.dat
```

```
  OUTPUTFILENAME            : Level2_MyOcean_13-04_27-04_2013.hdf5
  NEW_GRID_FILENAME         : Level2_.new
```

```
  FATHER_GEOMETRY           : ..\conv\MyOceanGeometry.dat
  !The geometry corresponds to thickness layers calculated based on depth centers cells
```

```
  NEW_GEOMETRY              : geometry_1.dat
  !The geometry corresponds to thickness layers converting the depth centers cells in faces top and then in thickness
```

```
  INTERPOLATION3D           : 1
  !3D vertical interpolation
```

```
  POLI_DEGREE               : 1
  DO_NOT_BELIEVE_MAP        : 0
  EXTRAPOLATE_2D            : 4
```

```
<end_file>
```

## Boundary Conditions

### Open boundary

Open boundary conditions are imposed in [hydrodynamic](../water/modules/hydrodynamic.md) module of MOHID. These differs from Level 1 to Level n of nested aplications. Level 1 is a 2D barotropic model usualy forced only with tide (no wind in atmosphere and InterfaceWaterAir module) with a slow connection. From Level 2 to Level n of nested domains the models are 3D baroclinic (e.g., include th density gradientes efects) and the Open Boundary Conditions (OBC) are resolved by imposing a Flow Relaxation Scheme (FRS) similar to the one presented by Marchesiello et al. (2001). The FRS is applied to temperature (T), salinity (S) and velocities (U, V) (Martinsen and Engedahl, 1987) being combined with a radiation scheme from Flather (1976) for the barotropic mode.

#### Level 1

Sample of the hydrodynamic file:

```
 BAROCLINIC                : 0
```

```
 SLOWSTART                 : 86400
 !Connect the tide slowly during the time period (seconds) considered above
 !This keyword is only used in spin-up run, remove it to continued the calculus from the spin-up run
```

```
 CONTINUOUS                : 0
 !Allows to pursue calculus from a previous run (should be 0 in spin-up run)
 !Change to CONTINUOUS: 1 to continued the calculus from the spin-up run
```

```
 INITIAL_ELEVATION         : 1
 INITIAL_ELEVATION_VALUE   : 1.8
 !Impose an initial elevation value of 1.8 m(must be in accordance with mean water leve value of the area)
 !These keywords is only used in spin-up run, remove it to continued the calculus from the spin-up run
```

```
 TIDE                      : 1
 TIDEPOTENTIAL             : 1

 DATA_ASSIMILATION         : 0
 BRFORCE                   : 0
 SUBMODEL                  : 0

 WIND                      : 0
 ATM_PRESSURE              : 0
```

```
 BIHARMONIC                : 1
 BIHARMONIC_COEF           : 2e9
 !The coefficient is calculated by BIHARMONIC_COEF = Dx^3/100
```

```
 RESIDUAL                  : 1
 ENERGY                    : 1
```

```
 OUTPUT_TIME               : 21600
```

#### Level 2

Sample of the hydrodynamic file:

```
 BAROCLINIC                : 1
 BAROCLINIC_POLIDEGREE     : 1
```

```
 ADV_METHOD_H              : 4
 ADV_METHOD_V              : 4
 TVD_LIMIT_H               : 4
 TVD_LIMIT_V               : 4
 VOLUME_RELATION_MAX       : 1.3
```

```
 CONTINUOUS                : 0
 !Allows to pursue calculus from a previous run (should be 0 in spin-up run)
 !Change to CONTINUOUS: 1 to continued the calculus from the spin-up run
```

```
 TIDE                      : 0
 TIDEPOTENTIAL             : 1
 DATA_ASSIMILATION         : 1
 BRFORCE                   : 1
 ATM_PRESSURE              : 1
```

```
 IMPOSE_INVERTED_BAROMETER : 1
 !Activate the Inverted Barometer method
```

```
 RADIATION                 : 2
 LOCAL_SOLUTION            : 7
 !OBC from submodel (Level 1) + field (MyOcean or RTOFS) + gauges (pressure effect from inverted barometer effect)
```

```
 SUBMODEL                  : 1
 MISSING_NULL              : 1
 SUBMODEL_EXTRAPOLATE      : 1
```

```
 SUBMODEL_FATHER_HOT_START : 1
 !Allows to start from water level of submodel (should be 1 in spin-up run)
 !Change it to SUBMODEL_FATHER_HOT_START : 0 to continued the calculus from the spin-up run
```

```
 EXTERNAL_BAROTROPIC_2D    : 0
 !Allows to calculate barotropic velocities from velocity U and V of the 3D external solution
```

```
 RAMP                      : 1
 RAMP_PERIOD               : 259200
 !Start with baroclinic force null and only after a specific period specified (seconds) the total force is compute
 !This keyword is only used in spin-up run, remove it to continued the calculus from the spin-up run
```

```
 WIND                      : 2
 WIND_SMOOTH_PERIOD        : 259200
 !Start wind with a ramp period for wind forcing in order to allow a slowly adjustement
 !Both keywords are used only in the spin-up run
 !Remove both keywords and include WIND: 1 to continued the calculus from the spin-up run
```

```
 ATM_PRESSURE              : 1
 ATM_PERIOD                : 259200
 !Connect pressure forcing slowly for the prescribe period of time
 !ATM_PERIOD keyword is used only in te spin-up run, remove it to continued the calculus from the spin-up run
```

```
 FLATHER_COLD_PERIOD       : 432000
 !Allows the relaxation of water level solution for a period of time
 !This keyword is only used in spin-up run, remove it to continued the calculus from the spin-up run
```

```
 BIHARMONIC                : 1
 BIHARMONIC_COEF           : 2e9
 !The coefficient is calculated by BIHARMONIC_COEF = Dx^3/100
```

```
 RESIDUAL                  : 1
 ENERGY                    : 1
```

```
 OUTPUT_TIME               : 0   10800
 SURFACE_OUTPUT_TIME       : 0   900
```

#### Level 3 to Level n

Sample of the hydrodynamic file:

```
 BAROCLINIC              : 1
 BAROCLINIC_POLIDEGREE   : 1
```

```
 ADV_METHOD_H             : 4
 ADV_METHOD_V             : 4
 TVD_LIMIT_H              : 4
 TVD_LIMIT_V              : 4
 VOLUME_RELATION_MAX      : 1.3
```

```
 CONTINUOUS                : 0
!Allows to pursue calculus from a previous run (should be 0 in spin-up run)
!Change to CONTINUOUS: 1 to continued the calculus from the spin-up run
```

```
 TIDE                      : 0
 TIDEPOTENTIAL             : 1
 DATA_ASSIMILATION         : 1
 BRFORCE                   : 1
 ATM_PRESSURE              : 1
```

```
 RADIATION                 : 2
 LOCAL_SOLUTION            : 2
 !OBC from submodel
```

```
 SUBMODEL                  : 1
 MISSING_NULL              : 1
 SUBMODEL_EXTRAPOLATE      : 1
```

```
 SUBMODEL_FATHER_HOT_START : 1
 !Allows to start from water level of submodel (should be 1 in spin-up run)
 !Change it to SUBMODEL_FATHER_HOT_START : 0 to continued the calculus from the spin-up run
```

```
 RAMP                      : 1
 RAMP_PERIOD               : 259200
 !Start with baroclinic force null and only after a specific period specified (seconds) the total force is compute
 !This keyword is only used in spin-up run, remove it to continued the calculus from the spin-up run
```

```
 WIND                      : 2
 WIND_SMOOTH_PERIOD        : 259200
 !Start wind with a ramp period for wind forcing in order to allow a slowly adjustement
 !Both keywords are used only in the spin-up run
 !Remove both keywords and include WIND: 1 to continued the calculus from the spin-up run
```

```
 ATM_PRESSURE              : 1
 ATM_PERIOD                : 259200
 !Connect pressure forcing slowly for the prescribe period of time
 !ATM_PERIOD keyword is used only in te spin-up run, remove it to continued the calculus from the spin-up run
```

```
 FLATHER_COLD_PERIOD       : 432000
 !Allows the relaxation of water level solution for a period of time
 !This keyword is only used in spin-up run, remove it to continued the calculus from the spin-up run
```

```
 BIHARMONIC                : 1
 BIHARMONIC_COEF           : 2e7
 !The coefficient is calculated by BIHARMONIC_COEF = Dx^3/100
```

```
 RESIDUAL                  : 1
 ENERGY                    : 1
```

```
 OUTPUT_TIME               : 0   10800
 SURFACE_OUTPUT_TIME       : 0   900
```

### Surface

To impose surface boundary conditions (e.g., wind stress, solar radiation, temperature, pressure etc.) in [atmosphere module](../water/modules/atmosphere.md) an example of the file is provided below.

Sample of the atmosphere file.

```
 OUTPUT_TIME             : 0. 3600.
```

```
<begin_rugosity>
 INITIALIZATION_METHOD   : CONSTANT
 DEFAULTVALUE            : 0.0025
<end_rugosity>
```

```
<beginproperty>
 NAME                    : wind velocity X
 UNITS                   : m/s
 DESCRIPTION             : wind velocity X interpolated from GFS model field
 DEFAULTVALUE            : 0.
 FILE_IN_TIME            : HDF
 FILENAME                : ..\..\..\GeneralData\GFS\Abril2013\Level2_10-04-2013_25-04-2013.hdf5
 TIME_SERIE              : 0
 OUTPUT_HDF              : 1
<endproperty>
```

```
<beginproperty>
 NAME                    : wind velocity Y
 UNITS                   : m/s
 DESCRIPTION             : wind velocity Y interpolated from GFS model field
 DEFAULTVALUE            : 0.
 FILE_IN_TIME            : HDF
 FILENAME                : ..\..\..\GeneralData\GFS\Abril2013\Level2_10-04-2013_25-04-2013.hdf5
 OUTPUT_HDF              : 1
<endproperty>
```

```
<beginproperty>
 NAME                    : air temperature
 UNITS                   : ºC
 DESCRIPTION             : Temperature interpolated from GFS model field
 DEFAULTVALUE            : 15.
 FILE_IN_TIME            : HDF
 FILENAME                : ..\..\..\GeneralData\GFS\Abril2013\Level2_10-04-2013_25-04-2013.hdf5
 TIME_SERIE              : 0
 OUTPUT_HDF              : 1
<endproperty>
```

```
<beginproperty>
 NAME                    : solar radiation
 UNITS                   : W/m^2
 DESCRIPTION             : Solar radiation interpolated from GFS model field
 DEFAULTVALUE            : 0.0
 FILE_IN_TIME            : HDF
 FILENAME                : ..\..\..\GeneralData\GFS\Abril2013\Level2_10-04-2013_25-04-2013.hdf5
 TIME_SERIE              : 0
 REMAIN_CONSTANT         : 0
 OUTPUT_HDF              : 1
<endproperty>
```

```
<beginproperty>
 NAME                    : atmospheric pressure
 UNITS                   : Pa
 DESCRIPTION             : Atmospheric pressure interpolated from GFS model field
 DEFAULTVALUE            : 0.
 FILE_IN_TIME            : HDF
 FILENAME                : ..\..\..\GeneralData\GFS\Abril2013\Level2_10-04-2013_25-04-2013.hdf5
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : relative humidity
 UNITS                   : fraction
 DESCRIPTION             : Constant value
 DEFAULTVALUE            : 0.55
 TIME_SERIE              : 0
 OUTPUT_HDF              : 1
<endproperty>
```

```
<beginproperty>
 NAME                    : cloud cover
 UNITS                   : %
 DESCRIPTION             : Constant value
 DEFAULTVALUE            : 50.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 1
<endproperty>
```

## Ocean-Atmosphere Heat Fluxes

Sensible/latent heat is computed in [Module InterfaceWaterAir](../water/modules/interfacewaterair.md) based in the atmospheric parameters (e.g.,wind stress, solar radiation, temperature, pressure etc.) prescribed in the atmosphere module. An example of the file is provided below.

Sample of the InterfaceWaterAir file.

```
 OUTPUT_TIME             : 0. 3600.
```

```
<begin_rugosity>
 INITIALIZATION_METHOD   : CONSTANT
 DEFAULTVALUE            : 0.0025
 REMAIN_CONSTANT         : 0
<end_rugosity>
```

```
<beginproperty>
 NAME                    : wind shear velocity
 UNITS                   : m/s
 DESCRIPTION             : Computed wind shear velocity
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : wind stress X
 UNITS                   : N/m2
 DESCRIPTION             : Computed wind stress X
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 TIME_SERIE              : 0
 DEFAULTVALUE            : 0.
 DEFINE_CDWIND           : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : wind stress Y
 UNITS                   : N/m2
 DESCRIPTION             : Computed wind stress Y
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 TIME_SERIE              : 0
 DEFAULTVALUE            : 0.
 DEFINE_CDWIND           : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : latent heat
 UNITS                   : W/m^2
 DESCRIPTION             : Computed latent heat
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : sensible heat
 UNITS                   : W/m^2
 DESCRIPTION             : Computed sensible heat
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : surface radiation
 UNITS                   : W/m^2
 DESCRIPTION             : Computed infrared radiation
 ALBEDO                  : 0.05
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : net long wave radiation
 UNITS                   : W/m^2
 DESCRIPTION             : Computed net long wave radiation
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : downward long wave radiation
 UNITS                   : W/m^2
 DESCRIPTION             : Computed downward long wave radiation
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : downward long wave radiation
 UNITS                   : W/m^2
 DESCRIPTION             : Computed downward long wave radiation
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : upward long wave radiation
 UNITS                   : W/m^2
 DESCRIPTION             : Computed upward long wave radiation
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

```
<beginproperty>
 NAME                    : non solar flux
 UNITS                   : W/m^2
 DESCRIPTION             : Computed infrared radiation
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 0
 DEFAULTVALUE            : 0.
 TIME_SERIE              : 0
 OUTPUT_HDF              : 0
<endproperty>
```

## Initial conditions and Spin up

For the spin-up procedure, a methodology based on a slow connection of the forcing terms (baroclinic force, winds tress, atmospheric pressure) is used. This methodology consists of defining an initial condition where the initial fields of salinity and temperature are interpolated from the external 3D solution (MyOcean or RTOFS), a null velocity field is assumed, and a SSH field with null gradient is also considered. A coefficient that varies linearly between 0 and 1 along the “connection” period of 5 days is multiplied by the baroclinic force and wind stress. Because the forces are slowly connected, the velocity reference solution of the OBC also needs to be slowly connected. The nudging term in the momentum equation is multiplied by a coefficient C given by:

$C=(Elapsed Time/Connection Period)^4;     Elapsed Time<Connection Period$

$C=1                                 ;     Elapsed Time>=Connection Period$

In this way, the velocity field near the boundary also converges slowly to the reference solution. To minimize the perturbations suffered by the initial condition of salinity and temperature along the spin-up period, a relaxation period variable in time was also assumed for these properties. The idea is to assume a relaxation period that increases with time; this way, in the beginning of the run the temperature and salinity fields have a stronger nudging when the external and internal activity is more intense due to the spin-up process. In the end of the spin-up period, the nudging in the model interior (out of the FRS area) is null. For the forces, a connection coefficient was assumed with a linear evolution over 5 days. For the reference solution a quadratic evolution was imposed. For the slow connection of the forcing mechanisms, the methodology followed by Mellor (2004) for the baroclinic force was assumed. This evolution allows, in the first instants, a strong nudging across the entire domain. With time, the model tends to be free except in the flow relaxation scheme area.

Initial conditions of salinity and temperature are imposed in [WaterProperties](../water/modules/waterproperties.md) module of MOHID. These differs from Level 1 to Level n of nested aplications. Level 1 is a 2D barotropic model and does not use the calculus of density fields (no keywords in WaterProperties file). From Level 2 to Level n of nested domains the models are 3D baroclinic (e.g., include th density gradientes efects) and the initial cnditons are prescribed by imposing a densitiy 3D field.

### Level 2

Sample of the WaterProperties file:

```
 OUTPUT_TIME              : 0 10800
 SURFACE_OUTPUT_TIME      : 0 900
```

```
 ADV_METHOD_H             : 4
 ADV_METHOD_V             : 4
 TVD_LIMIT_H              : 4
 TVD_LIMIT_V              : 4
 VOLUME_RELATION_MAX      : 1.3
```

```
<beginproperty>
 NAME                     : salinity
 UNITS                    : ºC
 DESCRIPTION              : MyOcean Interpolated results
 DEFAULTVALUE             : 34.895
```

```
 OLD                      : 0
 !Allows to pursue calculus from a previous run (should be 0 in spin-up run)
 !Change to OLD: 1 to continued the calculus from the spin-up run
```

```
 TYPE_ZUV                 : z
 INITIALIZATION_METHOD    : HDF
 FILENAME                 : ..\..\..\GeneralData\MyOcean\Level2_MyOcean_13-04_27-04_2013.hdf5
 !Read the 3D salinity fields interpolated to MOHID grid from HDF5 file
```

```
 ADVECTION_DIFFUSION      : 1
 DATA_ASSIMILATION        : 1
```

```
 BOUNDARY_CONDITION       : 4
 !Assumes null gradient at the open boundary
```

```
 OUTPUT_HDF               : 1
 OUTPUT_SURFACE_HDF       : 1
<endproperty>
```

```
<beginproperty>
 NAME                     : temperature
 UNITS                    : ºC
 DESCRIPTION              : MyOcean Interpolated results
 DEFAULTVALUE             : 16.
 OLD                      : 0
 TYPE_ZUV                 : z
 INITIALIZATION_METHOD    : HDF
 FILENAME                 : ..\..\..\GeneralData\MyOcean\Level2_MyOcean_13-04_27-04_2013.hdf5
 ADVECTION_DIFFUSION      : 1
 SURFACE_FLUXES           : 1
 DATA_ASSIMILATION        : 1
 BOUNDARY_CONDITION       : 4
 OUTPUT_HDF               : 1
 OUTPUT_SURFACE_HDF       : 1
<endproperty>
```

### Level 3 to Level n

Sample of the WaterProperties file:

```
 OUTPUT_TIME              : 0 10800
 SURFACE_OUTPUT_TIME      : 0 900
```

```
 ADV_METHOD_H             : 4
 ADV_METHOD_V             : 4
 TVD_LIMIT_H              : 4
 TVD_LIMIT_V              : 4
 VOLUME_RELATION_MAX      : 1.3
```

```
<beginproperty>
 NAME                     : salinity
 UNITS                    : ºC
 DESCRIPTION              : MyOcean Interpolated results
 DEFAULTVALUE             : 34.895
```

```
 OLD                      : 0
 !Allows to pursue calculus from a previous run (should be 0 in spin-up run)
 !Change to OLD: 1 to continued the calculus from the spin-up run
```

```
 TYPE_ZUV                 : z
 ADVECTION_DIFFUSION      : 1
 DATA_ASSIMILATION        : 1
```

```
 BOUNDARY_CONDITION       : 5
 !OBC are from Submodel
```

```
 SUBMODEL                 : 1
```

```
 SUBMODEL_INI             : 1
 !Allows to initializate the fields from submodel
 !Change to SUBMODEL_INI : 0 to continued the calculus from the spin-up run or remove it
```

```
 SUBMODEL_EXTRAPOLATE     : 1
 OUTPUT_HDF               : 1
 OUTPUT_SURFACE_HDF       : 1
<endproperty>
```

```
<beginproperty>
 NAME                     : salinity
 UNITS                    : ºC
 DESCRIPTION              : MyOcean Interpolated results
 DEFAULTVALUE             : 34.895
 OLD                      : 0
 TYPE_ZUV                 : z
 ADVECTION_DIFFUSION      : 1
 DATA_ASSIMILATION        : 1
 BOUNDARY_CONDITION       : 5
 SUBMODEL                 : 1
 SUBMODEL_INI             : 1
 SUBMODEL_EXTRAPOLATE     : 1
 OUTPUT_HDF               : 1
 OUTPUT_SURFACE_HDF       : 1
<endproperty>
```

## External fields assimilation

Assimilation of a 3D external solution salinity and temperature are imposed in [Assimilation](../water/modules/assimilation.md) module of MOHID. These differs from Level 1 to Level n of nested aplications. In level that external 3D fields are assimilated must be activated the cold relax period for assimilating the 3D salinity, temperature, velocity U and V.

### Level 1

Sample of Assimilation file:

```
  OUTPUT_TIME             : 0 864000.
```

```
<beginproperty>
  NAME                    : water level
  UNITS                   : m
  DIMENSION               : 2D
  OUTPUT_HDF              : 1
 <<begin_field>>
  DEFAULTVALUE            : 0
  INITIALIZATION_METHOD   : HDF
  FILE_IN_TIME            : HDF
  FILENAME                : ..\..\..\GeneralData\MyOcean\Level2_MyOcean_13-04_27-04_2013.hdf5
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : z
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

```
<beginproperty>
  NAME                    : velocity U
  UNITS                   : m/s
  DIMENSION               : 3D
  OUTPUT_HDF              : 1
  COLD_RELAX_PERIOD       : 432000
 <<begin_field>>
  DEFAULTVALUE            : 0
  INITIALIZATION_METHOD   : HDF
  FILE_IN_TIME            : HDF
  FILENAME                : ..\..\..\GeneralData\MyOcean\Level2_MyOcean_13-04_27-04_2013.hdf5
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : u
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

```
<beginproperty>
  NAME                    : velocity V
  UNITS                   : m/s
  DIMENSION               : 3D
  OUTPUT_HDF              : 1
  COLD_RELAX_PERIOD       : 432000
 <<begin_field>>
  DEFAULTVALUE            : 0
  INITIALIZATION_METHOD   : HDF
  FILE_IN_TIME            : HDF
  FILENAME                : ..\..\..\GeneralData\MyOcean\Level2_MyOcean_13-04_27-04_2013.hdf5
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : v
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

```
<beginproperty>
  NAME                    : temperature
  UNITS                   : ºC
  DIMENSION               : 3D
  OUTPUT_HDF              : 1
  COLD_RELAX_PERIOD       : 432000
 <<begin_field>>
  DEFAULTVALUE            : 18
  INITIALIZATION_METHOD   : HDF
  FILE_IN_TIME            : HDF
  FILENAME                : ..\..\..\GeneralData\MyOcean\Level2_MyOcean_13-04_27-04_2013.hdf5
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : z
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

```
<beginproperty>
  NAME                    : salinity
  UNITS                   : ºC
  DIMENSION               : 3D
  OUTPUT_HDF              : 1
  COLD_RELAX_PERIOD       : 432000
 <<begin_field>>
  DEFAULTVALUE            : 36
  INITIALIZATION_METHOD   : HDF
  FILE_IN_TIME            : HDF
  FILENAME                : ..\..\..\GeneralData\MyOcean\Level2_MyOcean_13-04_27-04_2013.hdf5
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : z
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

### Level 2

Sample of Assimilation file:

```
  OUTPUT_TIME             : 0 864000.
```

```
<beginproperty>
  NAME                    : water level
  UNITS                   : m
  DIMENSION               : 2D
  OUTPUT_HDF              : 1
 <<begin_field>>
  DEFAULTVALUE            : 0
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : z
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

```
<beginproperty>
  NAME                    : velocity U
  UNITS                   : m/s
  DIMENSION               : 3D
  OUTPUT_HDF              : 1
  COLD_RELAX_PERIOD       : 432000
 <<begin_field>>
  DEFAULTVALUE            : 0
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : u
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

```
<beginproperty>
  NAME                    : velocity V
  UNITS                   : m/s
  DIMENSION               : 3D
  OUTPUT_HDF              : 1
  COLD_RELAX_PERIOD       : 432000
 <<begin_field>>
  DEFAULTVALUE            : 0
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : v
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

```
<beginproperty>
  NAME                    : temperature
  UNITS                   : ºC
  DIMENSION               : 3D
  OUTPUT_HDF              : 1
  COLD_RELAX_PERIOD       : 432000
 <<begin_field>>
  DEFAULTVALUE            : 0
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : z
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

```
<beginproperty>
  NAME                    : salinity
  UNITS                   : ºC
  DIMENSION               : 3D
  OUTPUT_HDF              : 1
  COLD_RELAX_PERIOD       : 432000
 <<begin_field>>
  DEFAULTVALUE            : 0
  TYPE_ZUV                : z
 <<end_field>>
 <<begin_coef>>
  DEFAULTVALUE            : 1e9
  TYPE_ZUV                : z
  FILE_IN_TIME            : NONE
  REMAIN_CONSTANT         : 1
  INITIALIZATION_METHOD   : Sponge
  SPONGE_OUT              : 1e5
 <<end_coef>>
<endproperty>
```

## Usual mistakes and known limitations

If errors occurs after MOHID project implementation with the procedements provided above then check the following:

1\. The ratio between grids can not be higher than 1/5;

2\. In the open boundaries cannot exist intertidal areas;

3\. To choose the open boundaries delimitation avoid areas with big diferences in topography;

4\. Deviate the open boundaries far away from areas with big velocities jets;

5\. If noisy is created in the open boundary of large model domain (the ones that simulate the tide only) close some cells in the boundary with land;

6\. Check if the interpolations of 3D fields were performed with double precision;
