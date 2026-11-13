# HDF5Exporter

The application **HDF5Exporter** produces [time series](../io/formats/time-series.md) from data contained in [HDF5](hdf5.md) files in MOHID format .

Running options for this application are specified in an [input file](hdf5exporter.md) whose path is indicated in a [nomfich.dat file](../io/nomfich-dat-file.md).

## Introduction

HDF5Exporter extracts data contained in one or more HDF5 files subject to the parameter list, time window and locations specified by the user.

The time window specified by the user indicates the intended time series period. The actual time series period results from the intended period subjected to the data availability of the supplied HDF5 files.

Each time series location can consist in the grid cell that contains a point location or an area containing several grid cells. The output originates by HDF5Exporter is time series files. All the time series produced in a HDF5Exporter will have data for the same parameters and time window.

**Typical use:**

Extract time series from MOHID Water result files from several runs produced in continuous running of the model.

**Data input requirements:**

One or more HDF5 files, in MOHID format, containing data to create time series.

**Ouput:**

A time series file, in MOHID format, for each time series location defined by the user in the input file. In the input file the available option are the one described in the [Module TimeSerie](../water/modules/timeserie.md).

**Spatial definition**

The time series location can be a point (default)or an area.

**Area definition**

The area can be defined using a [polygon](../io/formats/polygons.md) file or a cell by cell mapping using a grid data file. For the area definition the additional keywords (options) are available (see below - Sample : Time Serie - Area Location) :

```
EXPORT_TYPE : 1/2 - 1 (default) point location and 2 area location
```

```
POLYGONS_FILE : **** - polygon file name. If the specify this keyword the area location will be defined as a polygon
```

```
MASK_GRID : *** - grid data file name where for each cell is given an integer value that is use to associate each cell to a specific timer serie
```

```
In the block use to define the location of each time serie the user must link each time serie with the polygons or grid data cells id's.
The user needs also to specify the correspondent layer (in the 3D cases).
```

```
<BeginTimeSerie>
...
MASK_ID : 1  ! In the case of polygons this id corresponds to the order in polygons file (1 - first polygon). In the case of grid data corresponds to the cell id.
LAYER   : KUB !  Grid layer. By default is assumed to be equal to the surface layer.
...
<EndTimeSerie>
```

## Input file

The name of the input file must be provided in the [nomfich.dat file](../io/nomfich-dat-file.md) in use.

```
(block for each HDF5 file to extract time series from; may be several, block order is
irrelevant)
<BeginHDF5File>
NAME                    : ... (path/name of HDF5 file with data to extract time series)
<EndHDF5File>

START_TIME              : ... (start time for time series: yyyy mm dd hh mm ss)
END_TIME                : ... (end time for time series: yyyy mm dd hh mm ss)

COMPUTE_RESIDUAL        : 0/1 (residual appearance in time series: 0 = appears not, 1 =
                               appears; 1 = default; residual is meanless in time series
                               produced, should always be 0)

VARIABLE_GRID           : 0/1 (use time variable grid: 0 = normal grid, 1 = variable grid;
                               0 = default)

WATERPOINTS_NAME        : ... (time independent map item name in HDF5 file)

(block for each parameter that is requested a time series; may be several)
<BeginParameter>
PROPERTY                : ... (parameter/property name)
HDF_GROUP               : ... (complete path in HDF5 file to parameter data)
<EndParameter>

(block for each time series location; may be several, one per output time series file)
<BeginTimeSerie>
NAME                    : ... (path/name of output time series file)

LOCALIZATION_I          : ... (time series location latitude cell indentification index, i)
LOCALIZATION_J          : ... (time series location longitude cell indentification index, j)
LOCALIZATION_K          : ... (time series location vertical cell indentification index, k;
                               should be 1 in 2D data)

The location of I and J could be also provided as longitude and latitude
<BeginTimeSerie>
NAME                    : ... (path/name of output time series file)

COORD_X                 : ... (time series location geographic decimated longitude)
COORD_Y                 : ... (time series location geographic decimated latitude)
LOCALIZATION_K          : ... (time series location vertical cell indentification index, k;
                               should be 1 in 2D data)

(if VARIABLE_GRID : 1:)
LATITUDE                : ... (time series location geographic decimated latitude)
LONGITUDE               : ... (time series location geographic decimated longitude)
<EndTimeSerie>
```

```
Remarks:

- option VARIABLE_GRID : 1 is usually useful when input HDF5 files contain satellite images
  data;

- the parameter/property name indicated in PROPERTY keyword should be the one used in MOHID
  framework, otherwise a warning message is issued;

- all the time series locations introduced through <BeginTimeSerie>...<EndTimeSerie> blocks
  must have valid data in the input HDF5 files for any output files to be produced: if one
  input location does not have valid data the application will produce an error and end
  without producing output files.

- The time-span of the data files must be included (but not necessarily equal)
  in the interval defined by START_TIME and END_TIME.

- The properties MUST be present in the data files, otherwise an error is returned.
```

## Sample : Time Serie - Area location

```
<BeginHDF5File>
NAME                    : K:\Simula\Nivel1CostaPort\Nivel2Tagus\res\WaterProperties_15.hdf5
<EndHDF5File>

<BeginHDF5File>
NAME                    : K:\Simula\Nivel1CostaPort\Nivel2Tagus\res\WaterProperties_16.hdf5
<EndHDF5File>

<BeginHDF5File>
NAME                    : K:\Simula\Nivel1CostaPort\Nivel2Tagus\res\WaterProperties_17.hdf5
<EndHDF5File>

START_TIME              : 2005 1 31 0 0 0
END_TIME                : 2005 3 1 0 0 0

COMPUTE_RESIDUAL        : 0
WATERPOINTS_NAME        : WaterPoints3D

EXPORT_TYPE             : 2
GRID_FILENAME           : batim.dat
POLYGONS_FILE           : poly.xy

<BeginParameter>
HDF_GROUP               : /Results/salinity
PROPERTY                : salinity
<EndParameter>

<BeginTimeSerie>
NAME                    : Monitoringbox1
COORD_X	         : -9.0
COORD_Y	         : 42.0
DEPTH_LEVEL             : 0
MASK_ID                 : 1
LAYER                   : 1
<EndTimeSerie>

<BeginTimeSerie>
NAME                    : Monitoringbox2
COORD_X	         : -9.1
COORD_Y	         : 42.1
DEPTH_LEVEL             : 0
MASK_ID                 : 2
LAYER                   : 1
<EndTimeSerie>

<BeginTimeSerie>
NAME                    : Monitoringbox3
COORD_X	         : -9.2
COORD_Y	         : 42.2
DEPTH_LEVEL             : 0
MASK_ID                 : 3
LAYER                   : 1
<EndTimeSerie>
```
