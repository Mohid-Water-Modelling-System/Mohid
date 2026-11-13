# ConvertToHDF5

The **ConvertToHDF5** is an application which allows the making of several operations, called **actions**, involving HDF5 files: conversion of data in other formats (e.g. NETCDF) to HDF5, grid interpolation, concatenation of several files.

Running options for this application are specified by the user in a input file named [**ConvertToHDF5Action.dat**](converttohdf5.md). Several actions can be specified in the same input file, being processed sequentially by the ConvertToHDF5 application.

## Introduction

The operations involving HDF5 files performed by ConvertToHDF5, specified individually by an action, can be organized in:

- file management
- grid interpolation and
- format conversion.

These types and the respective actions are detailed in the next sections.

The input file specification for each action can be found bellow in the Input file (ConvertToHDF5Action.dat) section.

### File management

#### Glue files

This action consists in two main apspects:

- joining or glue in a single HDF5 file two or more HDF5 files having the same HDF5 data groups and referring to time periods which come in sequence. Both sets of 2D and 3D HDF5 files can be glued. The use of the OPEN POINTS information is optional.

- merge HDF files groups from different files (e.g. merge groups from Hydrodynamic and Water Properties) for HDF5 files with same period and instants.

The separation between the two processes is done trough a keyword.

**Typical use:**

Glue MOHID Water results files from several runs produced in continuous running of the model, for storage space economy reasons. Can be used to join data from other origins (e.g. results of meteorological models) as long as the HDF5 format is the one supported by MOHID Water.

Also merge hydrodynamic and water properties HDF's (or other files) so that info is recorded in one file (e.g. Opendap) or to create hydrodynamics and properties animations at the same time (e.g. velocity vectors and properties) in MOHID Post

**Data input requirements:**

HDF5 files to be glued. For the concatenation of consecutive files "Grid" and "Results" data groups should be equal in all these files. For the merging of groups bathymetry and time instants need to be the same in the HDF5 files.

**Output:**

HDF5 file with glued "Results" data. "Residual" and "Statistics" HDF5 data groups are not copied to the output file since they are time period specific (different values potentially occour in each file). General statistics can be calculated for the glued HDF5 file data using tool [HDF5Statistics](hdf5statistics.md).

**ConvertToHDF5 action:** GLUES HDF5 FILES.

**RUNING TIP:**

Could be necessary to use a single compilation because time values get changed depending on single or double compilation.

### Grid interpolation

#### Interpolate files

This action performs the conversion of one HDF5 file data existing in one 2D or 3D spatial grid to another 2D or 3D spatial grid, creating a new HDF5 file. The interpolation is performed only for the data located a time window specified by the user.

The HDF5 file containing data to be interpolated is called the **father file**.

In case of 3D interpolation the application conducts first the horizontal grid interpolation (keeping father geometry) and only after it conducts the vertical interpolation (from father geometry to new geometry).

Several types of 2D interpolation are available for use: bilinear, spline 2D and triangulation. It also allows to produce a average in cell which returns one value based on all points inside cell. This method only works if there is at least one value per cell.

For vertical interpolation (used in 3D interpolation) can be supplied several polinomial degrees for interpolation.

**Typical use:**

Obtain an HDF5 file with data for forcing or providing initial conditions for a MOHID Water model, e.g. a meteorological forcing file.

**Data input requirements:**

For 2D/3D interpolation:

- father HDF5 file;

- father horizontal data grid, in a grid data file in the format supported by MOHID;

- new horizontal data grid, in a grid data file in the format supported by MOHID;

For 3D interpolation also needed:

- father vertical geometry, in a geometry file in the format supported by MOHID;

- new vertical geometry, in a geometry file in the format supported by MOHID;

- auxiliary horizontal data grid, in a grid data file in the format supported by MOHID; this file is used for horizontal grid interpolation in 3D interpolation operations.

**Ouput:**

HDF5 file with interpolated data. In case of 3D interpolation also produced an auxiliary HDF5 file with the result of the horizontal grid interpolation, which can be inspected to check if this operation is well performed.

**ConvertToHDF5 action:** INTERPOLATE GRIDS.

**RUNING TIP:**

Could be necessary to use a double compilation because values like precipitation and relative humidity can get out of scale values due to Rounding errors. For example precipitation can be slightly smaller than zero (which is negative) and relative humidity could be slightly higher than 1.

#### Patch files

This action consists in performing an interpolation of HDF5 data between grids, as in action Interpolate files, but considering more than one HDF5 file as containing data to be interpolated to the new grid and a priority scale. The interpolation is performed only for the data located in the time window specified by the user. The present version of this action operates only on 2D data.

Each HDF5 file containing data to be interpolated is called a **father file** and has an user-attributed **priority level** to be respected in the interpolation process: for each new grid cell the ConvertToHDF5 application will look for data first on the Level 1 father file and only in the case this data is inexistent will it look for data in Level 2 file, proceeding in looking for higher level files if no data is found subsequentely.

**Typical use:**

To obtain an HDF5 file with data from several HDF5 files each containing data with different spatial resolution and only for a specific part of the new grid. This is, for instance, the case when one is preparing a best resolution meteorological HDF5 file for forcing MOHID Water from several meteorological model domains, having different spatial resolution and span, since the best resolution data is not available for all new grid cells.

**Data input requirements:**

The new horizontal data grid, in a grid data file in the format supported by MOHID, and for each father file:

- level of priority: 1 = maximum priority, priority decreases with increasing level value;

- data grid, in the form of a grid data file in the format supported by MOHID.

**Ouput:**

HDF5 file with patched data.

**ConvertToHDF5 action:** PATCH HDF5 FILES.

### Format conversion

#### Meteorological model data

Mohid does not simulate explicitly the atmosphere, but needs information about atmospheric properties in time and space. This requires that atmospheric properties are supplied to MOHID Water in supported formats. These formats can be derived from meteorological data in HDF5 format. Because the results of meteorological models are accessed in different formats conversion is required.

The formats currently convertible to HDF5 in ConvertToHDF5 include the MM5 and the ERA40. These are succintly detailed in the next sections.

##### *ERA40*

This format refers to the European Centre for Medium-Range Weather Forecasts (ECMWF) 40 years re-analysises results, acessed by site [http://data.ecmwf.int/data/d/era40\_daily/](http://data.ecmwf.int/data/d/era40_daily/). This data is available for several meteorological variables with maximum 6 hour periodicity for days in the period from 1957-09-01 to 2002-08-31.

ERA40 data files are supplied by ECMWF in a NetCDF format and with an user-costumized time window, periodicity (time step range from 6 hours to a day) and meteorological properties set. The ERA40 meteorological properties which are recognized by MOHID are presented bellow together with the correspondent MOHID name:

```
---ERA40 NAME---         ---MOHID NAME---
  sshf                     sensible heat
  slhf                     latent heat
  msl                      atmospheric pressure
  tcc                      cloud cover
  p10u                     wind velocity X
  p10v                     wind velocity Y
  p2t                      air temperature
  ewss                     wind stress X
  nsss                     wind stress Y
```

The standard ConvertToHDF5 action is to convert to HDF5 the data referring to all MOHID Water recognized property available in the ERA40 file, producing an individual HDF5 file for each property. The name of each HDF5 file generated includes the ERA40 meteorological property identificator correspondent to the data contained.

Alternatively, ConvertToHDF5 can copy to a single ASCII file the heading information concerning each meteorological variable considered in the original ERA40 file.

**Typical use:**

Obtain an HDF5 file with data suitable for being used for forcing MOHID Water models.

**Data input requirements:**

ERA40 NetCDF file.

**Output:**

One HDF5 file for each meteorological property contained in the original NetCDF file.

**ConvertToHDF5 action:** CONVERT ERA40 FORMAT.

##### *Aladin*

This format relates to Aladin meteorological model results. Some of the atmospheric property needed by MOHID Water is present in Aladin output files, enabling to run prediction simulations with MOHID Water when access to Aladin prevision files is available.

The ConvertToHDF5 action converts Aladin results files from the original format to HDF5 format, allowing the easy use of these results in the MOHID framework. Conversion is only performed for the MM5 properties and the time window specified by the user.

**Typical use:**

Produce HDF5 meteorological data usable for forcing MOHID Water models.

**Data input requirements:**

Aladin netcdf results file to convert.

**Ouput:**

An HDF5 file with Aladin results and a grid data file in MOHID format with the Aladin grid pseudo-information: a fake orography is created of 100 m depth. This last file can be used to interpolate the Aladin data from the original grid to a new grid (see Interpolate files), for instance to produce an HDF5 file suitable for forcing MOHID Water models.

**Compilation:**

Caution! The ConvertToHDF5 executable must be compiled with the [Big-Endian](../reference/big-endian-little-endian.md) option set (see compatibility in the project's settings).

**ConvertToHDF5 action:** CONVERT ALADIN FORMAT.

##### *MM5*

This format relates to the Fifth-Generation NCAR / Penn State Mesoscale Model (MM5) output files format. Almost every atmospheric property needed by MOHID Water is present in MM5 output files.

The ConvertToHDF5 action converts MM5 results files from the original format to HDF5 format, allowing the easy use of these results in the MOHID framework. Conversion is only performed for the MM5 properties and the time window specified by the user.

Besides the conversion, the application can calculate some properties not contained in the MM5 files using the available information: these are wind stress, relative humidity and mean sea level pressure.

For conversion to be completed it is required the horizontal grid information of MM5 results which is available in special TERRAIN files.

**Typical use:**

Produce HDF5 meteorological data usable to force MOHID Water models.

**Data input requirements:**

MM5 results file to convert and MM5 TERRAIN file. The TERRAIN file supplies the MM5 results grid information.

**Ouput:**

A HDF5 file with MM5 results and a grid data file in MOHID format with the MM5 grid information. This last file can be used to interpolate the MM5 data from the original grid to a new grid (see Interpolate files), for instance to produce an HDF5 file suitable for forcing MOHID Water models.

**Compilation:**

Caution! The ConvertToHDF5 executable must be compiled with the [Big-Endian](../reference/big-endian-little-endian.md) option set (see compatibility in the project's settings).

**ConvertToHDF5 action:** CONVERT MM5 FORMAT.

##### *WRF*

This format relates to the Weather Research and Forecast model (WRF) output files format. Almost every atmospheric property needed by MOHID Water is present in WRF output files.

The ConvertToHDF5 action converts WRF results files from the original format to HDF5 format, allowing the easy use of these results in the MOHID framework. Conversion is only performed for the WRF properties and the time window specified by the user.

Besides the conversion, the application can calculate some properties not contained in the WRF files using the available information: these are wind stress, relative humidity and mean sea level pressure.

**Typical use:**

Produce HDF5 meteorological data usable to force MOHID Water models.

**Data input requirements:**

WRF direct output file.

**Ouput:**

An HDF5 file with WRF results and a grid data file in MOHID format with the WRF grid information. This last file can be used to interpolate the WRF data from the original grid to a new grid (see Interpolate files), for instance to produce an HDF5 file suitable for forcing MOHID Water models.

**Compilation:**

Caution! [PROJ4](proj4.md) library is needed! See Compile\_ModuleWRFFormat

**ConvertToHDF5 action:** CONVERT WRF FORMAT.

#### Wave model data

Mohid does not simulate explicitly the wave dynamics, but needs information about wave properties in time and space. This requires that wave properties are supplied to MOHID Water in supported formats. These formats can be derived from wave data in HDF5 format. Because the results of wave models are accessed in different formats conversion is required.

The formats currently convertible to HDF5 in ConvertToHDF5 include the table format of SWAN. This is succintly detailed in the next section.

##### *SWAN*

SWAN data files are supplied in a ASCII format. The data is organized by time instant. For each time instant there is a table of data. Each column represents one property.

The standard ConvertToHDF5 action is to convert to HDF5 the data referring to gradients of wave stresses X and Y, mean amplitude, mean period and mean direction.

**Typical use:**

Obtain HDF5 SWAN data usable for forcing MOHID Water models.

**Data input requirements:**

ASCII (table type) SWAN results data files and bathymetry in MOHID grid data foramt.

**Ouput:**

One HDF5 file containing all properties contained in the recognized set of properties (wave stresses X and Y, mean amplitude, mean period and mean direction).

**ConvertToHDF5 action:** CONVERT TO AND FROM SWAN.

#### Ocean model data

Ocean model data, available in diverse formats, can be used by MOHID Water to specify boundary (open ocean boundary and surface), initial conditions or for validation. These uses require that the model data is in HDF5 format and conversion is therefore needed.

Currently the large scale ocean models formats convertible into HDF5 by ConvertToHDF5 includes MERCATOR.

##### *MERCATOR*

MERCATOR data files are supplied in a NetCDF format and with an user-costumized spatial window and periodicity. Water level and water properties (temperature and salinity) data is available in type T files, velocity component u data is available in type U files and velocity component v data is available in type V files. The type of data of a specific MERCATOR file is generally indicated in the file name.

The standard ConvertToHDF5 action is to convert to HDF5 the data referring to temperature, salinity, water level, component u of velocity and component v of velocity.

**Typical use:**

Obtain HDF5 MERCATOR data usable for forcing or validation of MOHID Water models.

**Data input requirements:**

NetCDF MERCATOR results data files and NetCDF MERCATOR grid data files. It should be provided one grid data file of each type: T, U and V. These are generally provided by the MERCATOR services together with the results files.

**Ouput:**

One HDF5 file containing all properties contained in the recognized set of properties (temperature, salinity, water level, velocity u and velocity v) and the correspondent grid data and geometry files, containing respectively the horizontal grid and the vertical discretization of the HDF5 file. The grid data and geometry files can be used afterwards to interpolate the MERCATOR data to another grid and geometry (see Interpolate files).

**ConvertToHDF5 action:** CONVERT MERCATOR FORMAT.

#### Climatological data

Climatological data can be used in MOHID Water to specify boundary (open ocean boundary and surface), initial conditions or for validation, in case more realistic data (measurements or model) data is unavailable. This data is generally supplied by producers in formats not readly usable by MOHID Water which justifies the existence of a conversion tool.

Two climatological data format conversions are implemented in ConvertToHDF5: Levitus ocean data and Hellerman Rosenstein meteorological data.

##### *Levitus*

The Levitus climatology provides results for water temperature and salinity. The ConvertToHDF5 action converts the climatological data for the properties and spatial window requested by the user. Typically, it requires 3 steps to complete the task:

- convert levitus format

- extrapolate the data to the whole levitus domain(required to avoid uncoincidental coastlines)

- interpolate with the model grid(bathymetry)

**Typical use:**

Obtain climatological data in HDF5 format to use as boundary forcing and/or initial condition specification in MOHID Water models.

**Data input requirements:**

Levitus climatological data files, one per property and per time period (e.g a month).

**Ouput:**

HDF5 file with Levitus climatological data, grid data file with the horizontal grid of the data and a geometry file with vertical discretization of the data (MOHID formats). The grid data and the geometry files can be used to interpolate the climatological data from the original grid to a new grid (see Interpolate files).

**ConvertToHDF5 action:** CONVERT LEVITUS FORMAT.

##### *Hellerman Rosenstein*

This is a meteorological climatology providing wind stress. There is a file per wind stress component. Since the data refer to surface values it is a 2D field.

The ConvertToHDF5 action converts the climatological data for the properties and spatial window provided by the user.

**Typical use:**

Obtain climatological data in HDF5 format to use as meteorological forcing in MOHID Water models.

**Data input requirements:**

Hellerman Rosenstein climatological data ASCII files, one per wind stress component.

**Ouput:**

HDF5 file with Hellerman Rosenstein climatological data and grid data file with the horizontal grid of the climatological data. This grid data file can be used to interpolate the climatological data from the original horizontal grid to a new grid (see Interpolate files).

**ConvertToHDF5 action:** CONVERT HELLERMAN ROSENSTEIN ASCII.

##### *World Ocean Atlas 2005*

The World Ocean Atlas (WOA) 2005 climatology provides results for water temperature, salinity and several water quality and biology properties.

Description, Action and Input Files are described in a separate page: [ConvertToHDF5 WOA2005](converttohdf5-woa2005.md).

#### Generic NETCDF CF FILES

Converts any netcdf file following the CF convention into Mohid HDF5 files. Check the [Input File](converttohdf5.md) and [Sample File](converttohdf5.md).

#### Radar data

##### Seasonde RADAR Data

Converts a list of ascii files from the Seasonde radar data, which measures currents over a given region of space, into an hdf5 file with the MOHID format. Useful to validate coastal and estuarine models where radars are located. ConvertToHDF5 tool must be compiled with the preprocessor variable definition \_SCANLINE\_UNSAFE set in the Mohid Base 1.

### Input file (ConvertToHDF5Action.dat)

#### General structure

```
<begin_file> (block containing instructions for running a specific action)
ACTION                    : ... (intended action)
... (action specific instructions)
<end_file>

<begin_file>
ACTION                    : ...
...
<end_file>
```

#### GLUES HDF5 FILES

```
<begin_file>
ACTION                    : GLUES HDF5 FILES

GLUE_IN_TIME              : 0/1 (0 = merge files groups at same instants, 1 = concatenate consecutive files). Default is 1

3D_FILE                   : 0/1 (0 = 2D file, 1 = 3D file)

3D_OPEN                   : 0/1 (0 = default = open points are not used  , 1 = open points are used)

TIME_GROUP                : ... (Default="Time". Other option: "SurfaceTime".)

BASE_GROUP                : ... (Default="Results". Other options: "Residual", "SurfaceResults".)

OUTPUTFILENAME            : ... (path/name of HDF5 file to be created)

(block of HDF5 data files)
<<begin_list>>
... (path/name of HDF5 file with data to be included in glue, one per line, at least two files)
...
<<end_list>>
<end_file>
```

#### INTERPOLATE GRIDS

```
<begin_file>
ACTION                    : INTERPOLATE GRIDS

TYPE_OF_INTERPOLATION     : ... (type of horizontal interpolation: 1 = Bilinear, 2 = Spline2D,
                                 3 = Triangulation, 4 = Average in Cell)

INTERPOLATION_WINDOW      : ... ... ... ... (2D spatial window to consider for interpolation:
                                             Xmin Ymin Xmax Ymax; default = all domain)

START                     : ... (start date for output file: yyyy mm dd hh mm ss)
END                       : ... (end date for output file: yyyy mm dd hh mm ss)

INTERPOLATION3D           : 0/1 (0 = 2D interpolation, 1 = 3D interpolation)

FATHER_FILENAME           : ... (path/name of input HDF5 file with data to be interpolated)
FATHER_GRID_FILENAME      : ... (path/name of input grid data file with horizontal
                                 discretization of input HDF5 file)

OUTPUTFILENAME            : ... (path/name of output HDF5 file to be created)
NEW_GRID_FILENAME         : ... (path/name of input grid data file with horizontal
                                 discretization intended for output HDF5 file)

EXTRAPOLATE_2D            : 0/1/2/3/4/5 (2D extrapolation: 0=no extrapolation, 1=medium
                                     triangulation, 2=high triangulation,
                                     3=nearest neighbour, 4=nearest cell,
                                     5=constant value)

EXTRAPOLATE_VALUE         : ... (name of the value to extrapolate to when EXTRAPOLATE_2D is
                                 set to constant value (5))

DO_NOT_BELIEVE_MAP        : 0/1 (0=consider input HDF5 file map, 1=do not consider input HDF5
                                 file map)

BASE_GROUP                : ... (name of base group of HDF5 variables containing data to be
                                 interpolated; default is "/Results")

(if INTERPOLATION3D : 1 also required:)
FATHER_GEOMETRY           : ... (path/name of file (MOHID format) with vertical discretization
                                 of input HDF5 file)
NEW_GEOMETRY              : ... (path/name of file (MOHID format) with vertical discretization
                                 intended for output HDF5 file)
POLI_DEGREE               : 1/... (degree of vertical interpolation: 1=linear, ...)

AUX_GRID_FILENAME         : ... (path/name of input grid data file with horizontal
                                 discretization intended for auxiliar output HDF5 file;
                                 default is file provided in NEW_GRID_FILENAME)

AUX_OUTPUTFILENAME        : ... (path/name of auxiliar output HDF5 file to contain result
                                 of horizontal grid interpolation)
<end_file>
```

```
Remarks:

- the file indicated in AUX_GRID_FILENAME can be different from the one indicated in
  NEW_GRID_FILENAME in terms of bathymetry, while the horizontal grid should be, commonly, the
  same: this altered bathymetry can be used to extend the water column in the original data so
  that the process of vertical interpolation is done easily;

- in case of INTERPOLATION3D : 1, ConvertToHDF5 can generate new versions of bathymetry which
  are consistent with the geometry definition (extension is '.new'); there are possibly three
  bathymetry changes referring to father grid, new grid and aux grid (the same bathymetry is
  not altered twice); although initially new and aux grid are the same they can result
  different because of bathymetry changes;

- in case the new geometry is 2D and father geometry is 3D then POLI_DEGREE : 1
  (linear interpolation) should be used;

- EXTRAPOLATE_2D : 1/2/3/4/5 should be considered if it is expected that the coast line is not
  coincidental in the father and new grids, to avoid lack of data in the interpolation
  process; extrapolation is performed for all cells even the land cells;

- in case of DO_NOT_BELIEVE_MAP : 1 the application generates a map based on cells where
  interpolation results are available; this causes that if EXTRAPOLATE_2D : 1/2/3/4/5 is used
  the AUX_GRID_FILENAME should not have land cells in order for the new map to be concurrent
  with the result of extrapolation and avoid errors generation, specially if INTERPOLATION3D :
  1 is considered.
```

#### PATCH HDF5 FILES

```
<begin_file>
ACTION                    : PATCH HDF5 FILES

TYPE_OF_INTERPOLATION     : ... (type of interpolation: 3 = Triangulation, default and only
                                 one implemented)

START                     : ... (start date for output file: yyyy mm dd hh mm ss)
END                       : ... (end date for output file: yyyy mm dd hh mm ss)

(block for each father HDF5 file, should be at least two)
<<begin_father>>
LEVEL                     : ... (integer priority level: 1 = highest, increase for lower
                                 priority)
FATHER_FILENAME           : ... (path/name of input HDF5 file with data to be interpolated)
FATHER_GRID_FILENAME      : ... (path/name of input grid data file with horizontal
                                 discretization of input HDF5 file)
N_REMOVE_FRAME            : 0   (width, in number of cells, of frame to remove)
<<end_father>>

OUTPUTFILENAME            : ... (path/name of output HDF5 file to be created)
NEW_GRID_FILENAME         : ... (path/name of input grid data file with horizontal
                                 discretization intended for output HDF5 file)
<end_file>
```

#### CONVERT ERA40 FORMAT

```
<begin_file>
ACTION                    : CONVERT ERA40 FORMAT

FILENAME                  : ... (path/name of ERA40 NetCDF file)
OUTPUTFILENAME            : ... (path/name of HDF5 file to be created)
                                (root of name for all files produced)

CONVERT_TO_ASCII          : 0/1 (1 = convert variable heading info for ASCII file; 0 = default)
CONVERT_TO_HDF5           : 0/1 (1 = convert to HDF5 file; 0 = default)
GRIDTO180                 : 0/1 (1 = convert grid from [0 360] to [-180 180], 0 = default)

XX_VARIABLE               : ... (name of longitude variable in the input file: usual name
                                 is "longitude")
YY_VARIABLE               : ... (name of longitude variable in the input file: usual name
                                 is "latitude")
TIME_VARIABLE             : ... (name of time variable in the input file: usual name is
                                 "time")
<end_file>
```

```
Remarks:

- either CONVERT_TO_ASCII : 1 or CONVERT_TO_HDF5 : 1 must be chosen for any action to be
performed by ConvertToHDF5;

- when CONVERT_TO_HDF5 : 1 an HDF5 file is produced for every variable contained in the
original ERA40 file; the name of each file is composed of the name indicated on FILENAME
concatenated with the ERA40 variable identifier;

- to the XX_VARIABLE, YY_VARIABLE and TIME_VARIABLE keywords should generally be
specified "longitude", "latitude" and "time", respectively; the option to
include as keywords was made only to make the application robust to future variable name
changes.
```

#### CONVERT ALADIN FORMAT

```
<begin_file>
ACTION                    : CONVERT ALADIN FORMAT

OUTPUTFILENAME            : aladin.hdf5
OUTPUT_GRID_FILENAME      : aladin_griddata.dat

!Put here the name of any netcdf file for grid-data generation's sake.
INPUT_GRID_FILENAME      :   D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKCLOUD_OPASYMP_19723_20088.nc

<<begin_input_files>>
(path to aladin netcdf file)\ALADIN_BULKIR_OPASYMP_19723_20088.nc
...
<<end_input_files>>

<end_file>
```

```
Remarks:

- the name of each Aladin property to convert in <<begin_input_files>>...<<end_input_files>> block must conform to the following variables

          ---ALADIN NAME---    ---MOHID NAME---
            soclotot            CloudCover_
            sohumrel            RelativeHumidity_
            sofluxir            NonSolarFlux_
            sosspres            AtmosphericPressure_
            sosolarf            SolarRadiation_
            sotemair            AirTemperature_
            sowinmod            WindModulus_
            sowaprec            Precipitation_
            sozotaux            WindStressX_
            sometauy            WindStressY_
            sowindu10           WindVelocityX_
            sowindv10           WindVelocityY_
```

#### CONVERT MM5 FORMAT

```
<begin_file>
ACTION                    : CONVERT MM5 FORMAT

FILENAME                  : ... (path/name of MM5 file)
TERRAIN_FILENAME          : ... (path/name of MM5 TERRAIN file)

OUTPUTFILENAME               : ... (path/name of HDF5 file to be created)
OUTPUT_GRID_FILENAME         : ... (path/name of grid data file with horizontal grid of MM5 data to be created)

COMPUTE_WINDSTRESS           : 0/1 (1 = compute and write wind stress field; 0 = default)
COMPUTE_WINDMODULUS          : 0/1 (1 = compute wind modulus; 0 = default)
COMPUTE_RELATIVE_HUMIDITY    : 0/1 (1 = compute and write 2-meter relative humidity field; 0 = default)
COMPUTE_RELATIVE_HUMIDITY_3D : 0/1 (1 = compute and write 3D relative humidity field; 0 = default)
COMPUTE_PRECIPITATION        : 0/1 (1 = compute and write precipitation field; 0 = default)
COMPUTE_MSLP                 : 0/1 (1 = compute and write mean sea level pressure field; 0 = default)

WRITE_XYZ                    : 0/1 (1 = write xyz center grid cells; 0 = default)
WRITE_TERRAIN                : 0/1 (1 = write MM5 TERRAIN fields; 0 = default)

START                        : ... (start date for output file: yyyy mm dd hh mm ss)
END                          : ... (end date for output file: yyyy mm dd hh mm ss)

(block of MM5 properties to convert)
<<BeginFields>>
... (name of MM5 property to convert do HDF5 format, one per line)
...
<<EndFields>>

<end_file>
```

```
Remarks:

- the name of each MM5 property to convert in <<BeginFields>>...<<EndFields>> block must
conform to the MOHID designation specified in code of ModuleGlobalData; the correspondence is
the following (see Module_InterfaceWaterAir for a more detailed explanation).

          ---MM5 NAME---    ---MOHID NAME---
            T2             air temperature
            PSTARCRS       atmospheric pressure
            U10            wind velocity X
            V10            wind velocity Y
            UST            wind shear velocity
            LHFLUX         latent heat
            SWDOWN         sensible heat
            SWDOWN         solar radiation
            LWDOWN         infrared radiation
            SWOUT          top outgoing shortwave radiation
            LWOUT          top outgoing longwave radiation
            SOIL T 1       soil temperature layer 1
            SOIL T 1       soil temperature layer 2
            SOIL T 1       soil temperature layer 3
            SOIL T 1       soil temperature layer 4
            SOIL T 1       soil temperature layer 5
            SOIL T 1       soil temperature layer 6
            Q2             2-meter mixing ratio
            TSEASFC        sea water temperature
            PBL HGT        PBL height
            PBL REGIME     PBL regime
            RAIN CON       accumulated convective precipitation        (cm)
            RAIN NON       accumulated non-convective precipitation    (cm)
            GROUND T       ground temperature
            RES TEMP       infinite reservoir slab temperature
            U              wind velocity X_3D
            V              wind velocity Y_3D
            W              wind velocity Z_3D
            T              air temperature_3D
            PP             atmospheric pressure_3D
            Q              mixing ratio_3D
            CLW            cloud water mixing ratio_3D
            RNW            rain water mixing ratio_3D
            ICE            cloud ice mixing ratio_3D
            SNOW           snow mixing ratio_3D
            RAD TEND       atmospheric radiation tendency_3D
```

#### CONVERT WRF FORMAT

```
<begin_file>
ACTION                       : CONVERT WRF FORMAT

FILENAME                     : ... (path/name of WRF file)
OUTPUTFILENAME               : ... (path/name of HDF5 file to be created)
OUTPUT_GRID_FILENAME         : ... (path/name of grid data file with horizontal grid of WRF data
                                    to be created)

COMPUTE_WINDSTRESS           : 0/1 (1 = compute and write wind stress field; 0 = default)
COMPUTE_WINDMODULUS          : 0/1 (1 = compute wind modulus; 0 = default)
COMPUTE_RELATIVE_HUMIDITY    : 0/1 (1 = compute and write 2-meter relative humidity field; 0 = default)
COMPUTE_RELATIVE_HUMIDITY_3D : 0/1 (1 = compute and write 3D relative humidity field; 0 = default)
COMPUTE_PRECIPITATION        : 0/1 (1 = compute and write precipitation field; 0 = default)
COMPUTE_MSLP_MM5             : 0/1 (1 = compute mean sea level pressure with MM5toGrads algorithm; 1 = default)
COMPUTE_MSLP_WRF             : 0/1 (1 = compute mean sea level pressure with ARWPost algorithm; 0 = default)

WRITE_XYZ                    : 0/1 (1 = write xyz center grid cells; 0 = default)
WRITE_TERRAIN                : 0/1 (1 = write MM5 TERRAIN fields; 0 = default)
OUTPUT_DT                    : real (DT interval in seconds between each output; default = 0.0)

START                        : ... (start date for output file: yyyy mm dd hh mm ss)
END                          : ... (end date for output file: yyyy mm dd hh mm ss)

(block of properties to convert)
<<BeginFields>>
... (name of property to convert do HDF5 format, one per line)
...
<<EndFields>>

<end_file>
```

```
Remarks:

- the name of each property to convert in <<BeginFields>>...<<EndFields>> block must conform to the MOHID
  designation specified in code of ModuleGlobalData; the correspondence to WRF is in function
  VariableIsToRead in ModuleWRFFormat. Conversion is not direct since WRF outputs tendencies and not
  total values.

- Caution: This action is not extensively tested.
```

#### CONVERT MERCATOR FORMAT

```
<begin_file>
ACTION                    : CONVERT MERCATOR FORMAT

READ_OPTION               : 1/2/3/4 (version of MERCATOR files)

OUTPUTFILENAME            : ... (path/name of HDF5 file to be created)
OUTPUT_GRID_FILENAME      : ... (path/name of grid data with horizontal discretization to be
                                 created)
OUTPUT_GEOMETRY_FILENAME  : ... (path/name of geometry file with vertical discretization to be
                                 created)

(if READ_OPTION : 1:)
BASE_BULLETIN             : ...
DATES_FILE                : ...
NUM_DATES                 : ...

(if READ_OPTION : 2/3:)
INPUT_GRID_FILENAME       : ... (path/name of file with horizontal discretization of water
                                 properties and water level data)
(if READ_OPTION : 2:)
INPUT_GRID_FILENAME_U     : ... (path/name of file with horizontal discretization of velocity
                                 component U data)
INPUT_GRID_FILENAME_V     : ... (path/name of file with horizontal discretization of velocity
                                 component V data)

(if READ_OPTION : 3:)
INPUT_BATHY_FILENAME      : ... (path/name of file with bathymetry)

(if READ_OPTION : 3/4:)
CALC_BAROTROPIC_VEL       : 0/1 (1 = calculate barotropic velocity, 0 = not calculate;
                                 default = 0)

(if CALC_BAROTROPIC_VEL : 1 and READ_OPTION : 3:)
INPUT_MESH_ZGRID_FILENAME : ... (path/name of file with information about layers ticknesses)

(block of MERCATOR data files)
<<begin_input_files>>
... (path/name of MERCATOR NetCDF data file, one per line, can be several)
...
<<<end_input_files>>>

<end_file>
```

#### CONVERT LEVITUS FORMAT

```
<begin_file>
ACTION                    : CONVERT LEVITUS FORMAT

OUTPUTFILENAME            : ... (path/name of HDF5 file to be created)
OUTPUT_GRID_FILENAME      : ... (path/name of grid data with horizontal discretization to be
                                 created)
OUTPUT_GEOMETRY_FILENAME  : ... (path/name of geometry file with vertical discretization to be
                                 created)

PERIODICITY               : ... (periodicity of Levitus data: "monthly"/"annual"; default is
                                 "monthly")

SPATIAL_RESOLUTION        : ... (spatial resolution (degrees) of horizontal Levitus grid)

FILL_VALUE                : ... (real value identificator for missing data; default is
                                 "-99.999900")

(definition of spatial window to be present in output HDF5 file)
LOWER_LEFT_CORNER         : ... ... (longitude and latitude (degrees) of south west corner)
UPPER_RIGHT_CORNER        : ... ... (longitude and latitude (degrees) of north east corner)

(block for each water property to be present in output HDF5 file, can be several)
<<beginfield>>
NAME                      : ... (name of property)
ANNUAL_FILE               : ... (path/name of Levitus annual file)

(block of Levitus data files)
<<<begin_input_files>>>
... (path/name of Levitus data file (e.g. a monthly data file), one per line, can be several)
...
<<<end_input_files>>>

<<endfield>>

<end_file>
```

#### CONVERT HELLERMAN ROSENSTEIN ASCII

```
<begin_file>
ACTION                    : CONVERT HELLERMAN ROSENSTEIN ASCII

OUTPUTFILENAME            : ... (path/name of HDF5 file to be created)
OUTPUT_GRID_FILENAME      : ... (path/name of grid data with horizontal discretization to be
                                 created)

PERIODICITY               : ... (periodicity of Hellerman Rosenstein data: "monthly")

SPATIAL_RESOLUTION        : ... (spatial resolution (degrees) of horizontal Hellerman
                                 Rosenstein grid: default and only allowed value is "2.")

FILL_VALUE                : ... (real value identificator for missing data; default is
                                 "-99.999900")

(definition of spatial window to be present in output HDF5 file)
LOWER_LEFT_CORNER         : ... ... (longitude and latitude (degrees) of south west corner)
UPPER_RIGHT_CORNER        : ... ... (longitude and latitude (degrees) of north east corner)

(block for each Hellerman Rosenstein data file)
<<beginfield>>
NAME                      : ... (name of property: "wind stress X"/"wind stress Y")
FILE                      : ... (path/name Hellerman Rosenstein file)
<<endfield>>

<end_file>
```

#### CONVERT GENERIC NETCDF CF

```
<begin_file>

ACTION              : CONVERT NETCDF CF TO HDF5 MOHID
HDF5_OUT            : 0/1 (1 = convert to HDF5 file; 0 = default)
OUTPUTFILENAME      : ... (path/name of HDF5 file to be created)
NETCDF_OUT          : 0/1 (1 = convert to netcdf file; 0 = default)
OUTPUT_NETCDF_FILE  : ... (path/name of netcdf file to be created)
WINDOW_OUT          : imin imax jmin jmax (indexes of cells of the spatial window the user wants to extract, if this keyword in not used the entire domain is convert)

<<begin_time>>
NETCDF_NAME         : name of the netcdf property for time, generally "time"
<<end_time>>

<<begin_grid>>
NETCDF_NAME_LAT     : ... (name of the netcdf property for latitude, generally "latitude")
NETCDF_NAME_LONG    : ... (name of the netcdf property for longitude, generally "longitude")
NETCDF_NAME_MAPPING : ... (name of the netcdf field where is located the mapping 0 water 1 land)
                          (It can also be used any field, i.e. temperature)
MAPPING_LIMIT       : real (limit for netcdf value mapping; default = 0.5, thus >0.5 water)
INVERT_LAYER_ORDER  : 0/1  (by default (0) is assumed the last layer is the surface one. This keywpord allows invert the default reading order)
BATHYM_FROM_MAP     : 0/1  (based in the mask/mapping matrix and the field depth (NETCDF_NAME_DEPTH) the tool is able to generate a bathymetry)
BATHYM_FILENAME     :  ... (File name of the bathymetry that can be generate by tool (BATHYM_FROM_MAP : 1)
NETCDF_NAME_DEPTH   : ... ! (name of the netcdf property for the cells depth, generally "depth", usually is a  1D Field)
STARTS_180W         : 0/1 ! (default 0 - only important for global solutions. Some times a solution that starts in 0º is not a good option. This allows a global solution to start at 180 W)
<<end_grid>>

PROPERTIES_NUMBER   : real (number of netcdf properties)

<<begin_field>>
NETCDF_NAME         : name of the netcdf property
NAME                : name of the HDF5 Mohid property
UNITS               : property units
DESCRIPTION         : property description
DIM                 : 2/3 Number of property dimensions
VECTOR_INTENSITY    : 0/1 If property corresponds to a vector intensity
VECTOR_X            : Decomposition of Vector Intensity X
VECTOR_Y            : Decomposition of Vector Intensity Y
BEAUFORT_SCALE      : 0/1 If property correspond to the beaufort scale
ADD_FACTOR          : To add the property a fixed value
<<end_field>>

<<begin_field>>
...
<<end_field>>

<<begin_input_files>>
path to netcdf cf file
...
<<end_input_files>>
<end_file>
```

#### CONVERT SEASONDE RADAR FORMAT

```
<begin_file>
ACTION                    : CONVERT IH RADAR FORMAT

!HDF5 path and filename that will contain the results
OUTPUTFILENAME            : IHRADAR_201205_1400.hdf5

!griddata path and filename that will contain a griddata of the proxy bathymetry used in the hdf5 output file.
OUTPUT_GRID_FILENAME      : IHRADAR_griddata.dat

!Path and filename to the grid that maps the radar data. This grid must be constructed with the GIS
!or MOHID Studio and must "fit" well where the points stand.
INPUT_GRID_FILENAME       : ../IH-Radar-GridVersion4.grd

!This is the version of the grid used by the seasonde radar system.
!If the grid version of the seasonde radar changes, then most likely a new input grid must be
!constructed with the GIS or MOHID Studio.
IH_GRID_VERSION           : 4

!This is a list with the Seasonde radar ascii data files. One file per instant.
<<begin_input_files>>
../RadarFiles/TOTL_IHOC_2012_05_01_1400.tuv.txt
<<end_input_files>>
<end_file>
```

### Samples

All sample files are named *ConvertToHDF5Action.dat*.

#### Glue several MOHID(.hdf5) files sequentially (Concatenate)

```
<begin_file>
ACTION                   : GLUES HDF5 FILES

OUTPUTFILENAME           : SurfaceHydro_OP.hdf5

<<begin_list>>
D:\Projectos\SurfaceHydrodynamic_21.hdf5
D:\Projectos\SurfaceHydrodynamic_22.hdf5
<<end_list>>
<end_file>
```

#### Glue several MOHID(.hdf5) files with same instants (Merge)

```
<begin_file>
ACTION                   : GLUES HDF5 FILES

GLUE_IN_TIME             : 0

OUTPUTFILENAME           : SurfaceHydro_OP.hdf5

<<begin_list>>
D:\Projectos\SurfaceHydrodynamic_21.hdf5
D:\Projectos\SurfaceWaterProperties_21.hdf5
<<end_list>>
<end_file>
```

#### Interpolate 2D MOHID(.hdf5) files to a new grid

```
<begin_file>
ACTION : INTERPOLATE GRIDS

TYPE_OF_INTERPOLATION    : 1
FATHER_FILENAME          : D:\Projectos\MohidRun\test\res\Lagrangian_1.hdf5
OUTPUTFILENAME           : OilSpillThickness_GridRegular.hdf5

START                    : 2006 6 21 17 22 30
END                      : 2006 6 22 17 22 0

FATHER_GRID_FILENAME     : D:\Projectos\MohidRun\GeneralData\batim\Tagus.dat_A
NEW_GRID_FILENAME        : TagusConstSpacing.dat

BASE_GROUP               : /Results/Oil/Data_2D
<end_file>
```

#### Interpolate 3D MOHID(.hdf5) files to a new grid

```
<begin_file>
ACTION : INTERPOLATE GRIDS

TYPE_OF_INTERPOLATION   : 1
FATHER_FILENAME         : D:\Projectos\MohidRun\test\res\Lagrangian_1.hdf5
OUTPUTFILENAME          : OilSpillThickness_GridRegular.hdf5

START                   : 2006 6 21 17 22 30
END                     : 2006 6 22 17 22 0

FATHER_GRID_FILENAME    : D:\Projectos\MohidRun\GeneralData\batim\Tagus.dat_A
NEW_GRID_FILENAME       : TagusConstSpacing.dat

BASE_GROUP              : /Results/Oil/Data_2D

INTERPOLATION3D         : 1
FATHER_GEOMETRY         : D:\Projectos\MohidRun\test\data\Geometry_1.dat
NEW_GEOMETRY            : TagusGeometry.dat
AUX_GRID_FILENAME       : TagusConstSpacing.dat
AUX_OUTPUTFILENAME      : Aux_GridRegular.hdf5
<end_file>
```

#### Patch several MOHID(.hdf5) files to a new grid

```
<begin_file>
ACTION : PATCH HDF5 FILES

TYPE_OF_INTERPOLATION   : 3

START                   : 2005 2 28 13 0 0
END                     : 2005 3 1 13 0 0

<<begin_father>>
LEVEL                   : 3
FATHER_FILENAME         : K:\MM5output\2005022812_2005030712\MM5OUT_D1.hdf5
FATHER_GRID_FILENAME    : K:\MM5output\2005022812_2005030712\grid1.dat
<<end_father>>

<<begin_father>>
LEVEL                   : 2
FATHER_FILENAME         : K:\MM5output\2005022812_2005030712\MM5OUT_D2.hdf5
FATHER_GRID_FILENAME    : K:\MM5output\2005022812_2005030712\grid2.dat
<<end_father>>

<<begin_father>>
LEVEL                   : 1
FATHER_FILENAME         : K:\MM5output\2005022812_2005030712\MM5OUT_D3.hdf5
FATHER_GRID_FILENAME    : K:\MM5output\2005022812_2005030712\grid3.dat
<<end_father>>

OUTPUTFILENAME          : MM5Forcing.hdf5
NEW_GRID_FILENAME       : K:\Simula\GeneralData\Batim\CostaPortuguesa.dat
<end_file>
```

#### Convert an ERA40 file to MOHID(.hdf5)

```
<begin_file>
ACTION                  : CONVERT ERA40 FORMAT

FILENAME                : D:\Aplica\ERA40\1971ERA1973.nc
OUTPUTFILENAME          : D:\Aplica\ERA40\1971ERA1973T2

CONVERT_TO_ASCII        : 0
CONVERT_TO_HDF5         : 1

XX_VARIABLE             : longitude
YY_VARIABLE             : latitude
TIME_VARIABLE           : time
<end_file>
```

#### Convert a MM5 file to MOHID(.hdf5)

```
<begin_file>
ACTION                    : CONVERT WRF FORMAT
FILENAME                  : MMOUT_DOMAIN2
TERRAIN_FILENAME          : TERRAIN_DOMAIN2
OUTPUTFILENAME            : D2.hdf5
OUTPUT_GRID_FILENAME      : D2.dat
COMPUTE_WINDSTRESS        : 1
COMPUTE_RELATIVE_HUMIDITY : 1
COMPUTE_PRECIPITATION     : 1
COMPUTE_MSLP              : 1
COMPUTE_WINDMODULUS       : 1
WRITE_XYZ                 : 0
WRITE_TERRAIN             : 0
OUTPUT_DT                 : 3600

!START                     : 2010 02 15 00 00 00
!END                       : 2010 02 15 06 00 00

<<BeginFields>>
air temperature
atmospheric pressure
downward long wave radiation
latent heat
sensible heat
sea water temperature
solar radiation
top outgoing shortwave radiation
upward long wave radiation
wind velocity X
wind velocity Y
<<EndFields>>
<end_file>
```

#### Convert a WRF file to MOHID(.hdf5)

```
<begin_file>
ACTION                    : CONVERT WRF FORMAT
FILENAME                  : wrfout_d01_2010-02-15_00.nc
OUTPUTFILENAME            : wrfout_d01.hdf5
OUTPUT_GRID_FILENAME      : wrfout_d01.dat
COMPUTE_WINDSTRESS        : 1
COMPUTE_RELATIVE_HUMIDITY : 1
COMPUTE_PRECIPITATION     : 1
COMPUTE_MSLP              : 1
COMPUTE_WINDMODULUS       : 1
WRITE_XYZ                 : 0
WRITE_TERRAIN             : 0
OUTPUT_DT                 : 3600

START                     : 2010 02 15 00 00 00
END                       : 2010 02 15 06 00 00

<<BeginFields>>
air temperature
atmospheric pressure
downward long wave radiation
latent heat
sensible heat
sea water temperature
solar radiation
top outgoing shortwave radiation
upward long wave radiation
wind velocity X
wind velocity Y
<<EndFields>>
<end_file>
```

#### Convert Mercator-Ocean(.nc) to MOHID(.hdf5)

```
<begin_file>
ACTION                   : CONVERT MERCATOR FORMAT

READ_OPTION              : 2

OUTPUTFILENAME           : Psy2v2r1v_R20060628/MercatorR20060628.hdf5
OUTPUT_GRID_FILENAME     : Psy2v2r1v_R20060628/MercatorGridR20060628.dat
OUTPUT_GEOMETRY_FILENAME : Psy2v2r1v_R20060628/MercatorGeometryR20060628.dat

INPUT_GRID_FILENAME      : GridFiles/ist_meteog-gridT.nc
INPUT_GRID_FILENAME_U    : GridFiles/ist_meteog-gridU.nc
INPUT_GRID_FILENAME_V    : GridFiles/ist_meteog-gridV.nc

<<begin_input_files>>
Psy2v2r1v_R20060628/ist_meteog-mercatorPsy2v2r1v_T_MEAN_ANA_20060621_R20060628.nc
Psy2v2r1v_R20060628/ist_meteog-mercatorPsy2v2r1v_T_MEAN_ANA_20060622_R20060628.nc
Psy2v2r1v_R20060628/ist_meteog-mercatorPsy2v2r1v_T_MEAN_ANA_20060623_R20060628.nc
<<end_input_files>>
<end_file>
```

#### Convert Levitus format to MOHID(.hdf5) and interpolate grid

##### Convert

First convert the Levitus ASCII format to a raw HDF5 format:

```
<begin_file>
ACTION                   : CONVERT LEVITUS FORMAT

OUTPUTFILENAME           : Levitus.hdf5
OUTPUT_GRID_FILENAME     : LevitusGrid.dat
OUTPUT_GEOMETRY_FILENAME : LevitusGeometry.dat

PERIODICITY              : monthly
SPATIAL_RESOLUTION       : 0.25
FILL_VALUE               : -99.9999

LOWER_LEFT_CORNER        : -16.0  31
UPPER_RIGHT_CORNER       :   1.   40

<<beginfield>>
NAME                     : salinity
ANNUAL_FILE              : DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s000hr.obj

<<<begin_input_files>>>
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s001
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s002
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s003
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s004
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s005
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s006
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s007
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s008
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s009
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s010
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s011
DataCenter\DadosBase\Ocean\Levitus\Data\Salinity\s012
<<<end_input_files>>>
<<endfield>>

<<beginfield>>
NAME                     : temperature
ANNUAL_FILE              : DataCenter\DadosBase\Ocean\Levitus\Data\Temp\t000hr.obj

<<<begin_input_files>>>
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t001
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t002
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t003
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t004
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t005
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t006
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t007
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t008
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t009
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t010
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t011
DataCenter\DadosBase\Ocean\Levitus\Data\Temperature\t012
<<<end_input_files>>>
<<endfield>>
<end_file>
```

##### Extrapolate

Then extrapolate the data (still in the raw HDF5 format):

```
<begin_file>
ACTION                   : INTERPOLATE GRIDS

TYPE_OF_INTERPOLATION    : 1

FATHER_FILENAME          : Levitus.hdf5
FATHER_GRID_FILENAME     : LevitusGrid.dat

OUTPUTFILENAME           : LeviTusAllPointsWithData.hdf5
NEW_GRID_FILENAME        : LevitusGrid.dat

START                    : -9999 1  1 0 0 0
END                      : -9999 12 1 0 0 0

INTERPOLATION3D          : 1
FATHER_GEOMETRY          : LevitusGeometry.dat
NEW_GEOMETRY             : LevitusGeometry.dat
AUX_GRID_FILENAME        : LevitusGrid.dat
AUX_OUTPUTFILENAME       : AuxLeviTusAllPointsWithData.hdf5

POLI_DEGREE              : 3
DO_NOT_BELIEVE_MAP       : 1
EXTRAPOLATE_2D           : 2
<end_file>
```

##### Interpolate

Finally, interpolate to the final grid and geometry (same as the Interpolate 3D sample):

```
<begin_file>
ACTION                   : INTERPOLATE GRIDS

TYPE_OF_INTERPOLATION    : 1
FATHER_FILENAME          : LeviTusAllPointsWithData.hdf5
OUTPUTFILENAME           : CadizMonthlyLevitus.hdf5
FATHER_GRID_FILENAME     : LevitusGrid.dat
NEW_GRID_FILENAME        : Algarve0.02SigmaSmooth_V3_CartMoreLayers.dat

START                    : -9999 1  1 0 0 0
END                      : -9999 12 1 0 0 0

INTERPOLATION3D          : 1
FATHER_GEOMETRY          : LevitusGeometry.dat
NEW_GEOMETRY             : Geometry_1.dat
AUX_OUTPUTFILENAME       : AuxCadizMonthlyLevitus.hdf5
AUX_GRID_FILENAME        : Aux12km.dat

POLI_DEGREE              : 3
DO_NOT_BELIEVE_MAP       : 1
<end_file>
```

Note that the programme may construct a new bathymetry twice. Use this bathymetry only on the AUX\_GRID\_FILENAME keyword.

#### Convert Hellerman Rosenstein ASCII format to MOHID(.hdf5)

```
<begin_file>
ACTION                   : CONVERT HELLERMAN ROSENSTEIN ASCII

OUTPUTFILENAME           : ClimatologicWindStress.hdf5
OUTPUT_GRID_FILENAME     : ClimatologicWindStressGrid.dat

PERIODICITY              : monthly
SPATIAL_RESOLUTION       : 2.
FILL_VALUE               : -99.9999

LOWER_LEFT_CORNER        : -180  -90
UPPER_RIGHT_CORNER       : 180  90

<<beginfield>>
NAME                     : wind stress X
FILE                     : D:\Aplica\Dados\Hellerman_Rosenstein\TAUXX.DAT
<<endfield>>

<<beginfield>>
NAME                     : wind stress Y
FILE                     : D:\Aplica\Dados\Hellerman_Rosenstein\TAUYY.DAT
<<endfield>>
<end_file>
```

#### Convert ALADIN(.nc) format to MOHID(.hdf5)

```
<begin_file>
ACTION                    : CONVERT ALADIN FORMAT

OUTPUTFILENAME            : aladin.hdf5
OUTPUT_GRID_FILENAME      : aladin_griddata.dat

!Put here the name of any netcdf file for grid-data generation's sake.
INPUT_GRID_FILENAME      :   D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKCLOUD_OPASYMP_19723_20088.nc

<<begin_input_files>>
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKIR_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKPRES_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKSOLAR_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKTAIR_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKWIND_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_FLUXPRE_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_STRESSU_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_STRESSV_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_U10_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_V10_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKCLOUD_OPASYMP_19723_20088.nc
D:\Aplica\BiscayAplica\FORCAGES\METEO\ALADIN_BULKHUMI_OPASYMP_19723_20088.nc
<<end_input_files>>

<end_file>
```

#### Convert generic Netcdf CF to MOHID(.hdf5)

This sample is for converting a typical GFS file

```
<begin_file>

ACTION : CONVERT NETCDF CF TO HDF5 MOHID
HDF5_OUT : 1
OUTPUTFILENAME : ..\..\data\hdf5\GFS.hdf5
NETCDF_OUT : 1
OUTPUT_NETCDF_FILE : ..\..\data\netcdf\GFS.nc
<<begin_time>>
NETCDF_NAME : time
<<end_time>>

<<begin_grid>>
NETCDF_NAME_LAT : latitude
NETCDF_NAME_LONG : longitude
NETCDF_NAME_MAPPING : TMP_10mb
MAPPING_LIMIT : -10000

<<end_grid>>

PROPERTIES_NUMBER : 6

<<begin_field>>
NETCDF_NAME : PRMSL_meansealevel
NAME : atmospheric pressure
UNITS : Pa
DESCRIPTION : GFS
DIM : 2
<<end_field>>

<<begin_field>>
NETCDF_NAME : UGRD_10mb
NAME : wind velocity X
UNITS : m/s
DESCRIPTION : MOHID
DIM : 2
<<end_field>>

<<begin_field>>
NETCDF_NAME : VGRD_10mb
NAME : wind velocity Y
UNITS : m/s
DESCRIPTION : MOHID
DIM : 2
<<end_field>>

<<begin_field>>
NETCDF_NAME              : wind_modulus
NAME                     : wind modulus
UNITS : m/s
DESCRIPTION : MOHID
DIM : 2
VECTOR_INTENSITY         : 1
VECTOR_X                 : wind velocity X
VECTOR_Y                 : wind velocity Y
<<end_field>>

<<begin_field>>
NETCDF_NAME              : wind_modulus_beaufort
NAME                     : wind modulus beaufort
UNITS : beaufort scale
DESCRIPTION : MOHID
DIM : 2
BEAUFORT_SCALE           : 1
VECTOR_X                 : wind modulus
<<end_field>>

<<begin_field>>
NETCDF_NAME : TMP_2maboveground
NAME : air temperature
UNITS : oC
DESCRIPTION : MOHID
DIM : 2
ADD_FACTOR  : -273
<<end_field>>

<<begin_input_files>>
..\..\data\netcdf\gfs_4_20111115.nc
<<end_input_files>>
<end_file>
```

#### Convert Seasonde RADAR data sample configuration file

The input grid file must be created to be the recipient of the data in the ascii files. To create this grid, it is recommended to convert manually, resorting to a good text editor, such as Notepad++, a few radar ascii files into the xyz Mohid format, then load them in the GIS or Mohid Studio, and from there create and save a regular grid that fits the points. Use that grid to convert the data into one hdf5 file.

```
<begin_file>
ACTION                    : CONVERT IH RADAR FORMAT

OUTPUTFILENAME            : IHRADAR_201205_1400.hdf5
OUTPUT_GRID_FILENAME      : IHRADAR_griddata.dat
INPUT_GRID_FILENAME       : ../IH-Radar-GridVersion4.grd
IH_GRID_VERSION           : 4

<<begin_input_files>>
../RadarFiles/TOTL_IHOC_2012_05_01_1400.tuv.txt
../RadarFiles/TOTL_IHOC_2012_05_30_1400.tuv.txt
<<end_input_files>>
<end_file>
```

Compile the ConvertToHDF5 tool with the preprocessor variable \_SCANLINE\_UNSAFE defined in Mohid Base 1.

### OceanColor modules compilation

Compiling the **ConvertToHDF5** tool with the OceanColor modules is more complicated than one might expect. A solution is proposed here for a release version using the Compaq Visual Fortran 6.6c. The difficulties rise because C code is embedded with a fortran interface and also, extra libraries such as hdf4 are required.

#### Pre-requisites

This is a list of prerequisites to successfully compile the tool:

- Compaq Visual Fortran 6.5 with patch 6.6c,
- VS .NET 2003 (Vc7 in particular),
- Hdf5 libraries (**hdf5.lib** **hdf5\_fortran.lib** **hdf5\_hl.lib**),
- Netcdf libraries (**netcdf.lib** **netcdf\_.lib**),
- Hdf4 libraries (**hd421.lib**, **hm421.lib**),
- szlib, zlib and jpeg libraries (**szlib.lib**, **zlib.lib** and **libjpeg.lib**),
- the fortran source files (**ModuleConvertModisL2.F90 ModuleConvertModisL3.F90 ModuleConvertOceanColorL2.F90**),
- the C source files and their fortran interface files (**readL2scan.c readL2Seadas.c** and **cdata.f crossp.f fgeonav.f**).

#### CVF IDE configuration

1. Configure everything as specified in [Compiling with CVF](../development/build/compiling-with-cvf.md).
2. Add the source files listed in the prerequisites above to the source files listing.
3. Go to **Tools-->Options...-->Directories**. There, add the **$DOTNET2K3/Vc7/bin** to the **Executable files'**; the **$DOTNET2K3/Vc7/include** and **$DOTNET2K3/Vc7/PlatformSDK/include** to the **Include files**; and finally, the **$DOTNET2K3/Vc7/lib**, **$DOTNET2K3/Vc7/PlatformSDK/lib** and **$DOTNET2K3/Vc7/PlatformSDK/bin** to the **Library files**.
4. Go to **Projects-->Settings-->Release-->Link-->Input**. There, add the following libraries: **netcdf.lib netcdf\_.lib hd421.lib hm421.lib libjpeg.lib**. (Make sure the hdf5 libraries as well as the szlib and zlib libraries are already mentioned).

#### Troubleshoots

**Q: I get unresolved external references during linkage, but I have all the libraries mentioned above included. What should I do?**

A: Unresolved external references can come out for two reasons:

1. you didn't specified all the libraries required or all the paths for the default libraries or,
2. [name mangling](http://en.wikipedia.org/wiki/Name_decoration) problems. Use the [dumpbin](../building/windows/dumpbin.md) utility to the libraries to checkout which language convention they are using. If that's the problem then you need to try to get new libraries with the correct naming convention.

That's it, you should now be able to build the [ConvertToHdf5](converttohdf5.md) project successfully.

**Q: I got a message saying the entry point \_NF\_PUT\_ATT\_REAL@28 could not be located in netcdf.dll**

A: copy the file netcdf.dll to the exe folder

### References

- [HDF5 Homepage](http://www.hdfgroup.org/)
- [HDF4 Homepage](http://www.hdfgroup.org/)

### Links

- [Module\_Atmosphere](../water/modules/atmosphere.md)
- [Module\_InterfaceWaterAir](../water/modules/interfacewaterair.md)
- [Coupling\_Water-Atmosphere\_User\_Manual](../water/coupling-water-atmosphere-user-manual.md)
