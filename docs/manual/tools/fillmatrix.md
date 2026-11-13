# FillMatrix

## Introduction

The **FillMatrix** application can create a [grid data](../io/formats/grid-data.md) file or [HDF5](hdf5.md) file from two possible origins:

1. data interpolated from time invariant points [XYZ](../io/formats/xyz-points.md) or time variant points [Time Series](../io/formats/time-series.md).
2. data not interpolated but used to fill the output grid or hdf5 according to grid id's. Also works with time invariant (single value) or time variant [Time Series](../io/formats/time-series.md) data.

FillMatrix should not be confused with [Module FillMatrix](../base/module-fillmatrix.md) which is responsible in MOHID for handling data input to 2D or/and 3D arrays.

Running options for this application are specified by the user in a input file named [FillMatrix.dat](fillmatrix.md).

## How it works?

This application can work with two different kind of data: with interpolation (from georeferenced data) and without interpolation (data to create output is associated with a property class and does not has a defined space coordinate)

### With Interpolation

In this process georeferenced data is used in an interpolation process to obtain one or a time series of 2D fields usable in the MOHID framework. The interpolation process involves spatial interpolation, in the case of a time invariant 2D field, and spatial and time interpolation, in the case of a time series of 2D fields.

The georeferenced data can be of two types:

- time invariant: XYZ data, defined in a XYZ file or directly by the user in the input file;

- time variant: time series data, defined in a time series file.

The data location defined individually by the user in the input file is called a **station**.

The output field is obtained for the desired property and can also be time invariant or time variant:

- a **time invariant field (grid data file)** is obtained from spatial interpolation and time interpolation for a time instant specified by the user;

- a **time variant field (HDF5 file)** is determined taking into account user-defined time window, a maximum time span and output times; the maximum time span is used in case of data lacking to control the process of interpolation: if a source of time variant data has data more distant in time from the desired output time than the maximum time span then this source is not considered in the interpolation.

The interpolation in space can be accomplished using triangulation or [inverse-distance weighted method](../reference/inverse-distance-weighted-method.md). In this second method the user can specify the inverse-distance weighting power and the maximum distance between the intended grid location and a data station for that station to be considered in the spatial interpolation.

The interpolation in time is linear.

**Typical use:**

Generate, from measurement networks data, initial or forcing data for MOHID Water models.

**Data input requirements:**

A [grid data](../io/formats/grid-data.md) file containing grid and mapping information for the intended output file. This is commonly a bathymetry file.

One or more files containing data to use as the interpolation basis, which can consist in:

- one [XYZ](../io/formats/xyz-points.md) file, if the intended output data field is not time variable;

- one or more [Time Series](../io/formats/time-series.md) files, each representing a station, usable both in cases of time variable or non variable fields.

Alternatively the use of data files can be ignored, being the data values for each station supplied by the user in the input file.

The coordinate system must be the same in all files.

**Ouput:**

An grid data or a HDF5 file, respectively in case a time invariant and a time variant field is generated.

### Without Interpolation

In this process data is not georeferenced and is used to obtain one grid or a time series of 2D fields usable in the MOHID framework. This method (without interpolation) is used when one property has different classes spatially distributed and each class has its own time invariant or time variant values. FillMatrix searches for the values associated with each class and fills cell values. The link between the cell classes and the fill value for the class is a unique ID.

As an example, vegetation property as LAI can be a MOHID Land model input and is defined for each vegetation type which means that different vegetation types will have different LAI values or evolution. FillMatrix, in this case, will need to run trough a grid with vegetation ID's (to know which vegetation type is growing in the cell) and fills the output cells with the correspondent values for LAI (constant value or timeserie value defined for that ID/vegetation type).

The fill values must be defined for all the property classes present and can be of two types:

- time invariant or constant values, if the property class is constant in time;

- time variant: time series data, defined in a time series file.

The fill values defined individually for each class by the user in the input file is called a **station**. (Could be named class but it would duplicate code).

The output field is obtained for the desired property and can also be time invariant or time variant:

- a **time invariant field (grid data file)** is obtained spatially according to classes ID's and time interpolation for a time instant specified by the user;

- a **time variant field (HDF5 file)** is determined taking into account user-defined time window, a maximum time span and output times; the maximum time span is used in case of data lacking to control the process of interpolation: if a source of time variant data has data more distant in time from the desired output time than the maximum time span then this source is not considered in the interpolation.

**Typical use:**

This process is typically used when the user has to input a property to the model and this property has different classes spatially distributed where each has its own values or time evolution.

**Data input requirements:**

A [grid data](../io/formats/grid-data.md) file containing grid and mapping information for the intended output file. This is commonly a bathymetry file.

A [grid data](../io/formats/grid-data.md) file containing property types or classes in a form of an ID.

One or more files containing data to use as the interpolation basis, which can consist in one or more [Time Series](../io/formats/time-series.md) files, each representing a station, usable both in cases of time variable or non variable fields.

Alternatively the use of data files can be ignored, being the data values for each station supplied by the user in the input file (constant values).

**Ouput:**

An grid data or a HDF5 file, respectively in case a time invariant and a time variant field is generated.

## Input file (FillMatrix.dat)

```
PROPERTY_NAME            : ... (name of the property whose values are being interpolated)
PROPERTY_UNITS           : ... (units of property whose values are being interpolated)

OUTPUT_FILE              : ... (path/name of output data file)

VARIABLE_IN_TIME         : 0/1 (states if fields to be created are variable in time:
                                0 = not variable (default), 1 = variable)
                               (if they are variable in time a HDF5 file is created)

GRID_DATA_FILE           : ... (path/name of input grid data file, a bathymetry/topography file)

INTERPOLATION_METHOD     : ... (spatial interpolation method: 0 = no interpolation,
                                 1 = triangulation, 2 = IWD)

(if  INTERPOLATION_METHOD = no interpolation)
FILL_ID_FILE             : ... (path/name of input data file with classes ID)

XYZ_FILE                 : ... (path/name of input XYZ file, in Mohid GIS format)

(if XYZ_FILE not provided:)
(block for each station containing data to be considered; may be several)
<begin_station>
NAME                     : ... (name of station)
(if  INTERPOLATION_METHOD /= no interpolation)
X                        : ... (XX coordinate of the station)
Y                        : ... (YY coordinate of the station)
(if  INTERPOLATION_METHOD = no interpolation)
ID                       : ... (ID associated to values; must correspond to ID's from FILL_ID_FILE)

VALUE_TYPE               : ... (type of information: "SINGLE_VALUE" (default), "TIMESERIE")
(if VALUE_TYPE : SINGLE_VALUE:)
VALUE                    : ... (property value of the station)
(if VALUE_TYPE : TIMESERIE:)
FILENAME                 : ... (path/name of station time series file)
DATA_COLUMN              : ... (column where data is stored)
<end_station>

(if VARIABLE_IN_TIME : 1:)
START                    : ... (start time: yyyy mm dd hh mm ss)
END                      : ... (end time: yyyy mm dd hh mm ss)
MAX_TIME_SPAN            : ... (maximum span in seconds for a time serie value to be
                                considered valid as data source)
OUTPUT_TIME              : ... (output time in seconds, may be provided several values)
SKIP_NULLVALUES          : 0/1 (option to do not output null values: 0 = no (default), 1 = yes)

(if VARIABLE_IN_TIME : 0:)
TIME_INSTANT             : ... (time instant in which to extract data from time series:
                                yyyy mm dd hh mm ss)

(if INTERPOLATION_METHOD : 1:)
FILL_OUTSIDE_POINTS      : 0/1 (option to perform extrapolation in triangulation: 0 = no
                                (default), 1 = yes)
WRITE_TRIANGLES          : 0/1 (option to write a polygon file with created triangles in
                                triangulation: 0 = no (default), 1 = yes)
TRIANGLES_FILE           : ... (path/name of triangles file to create)

(if INTERPOLATION_METHOD : 2:)
MAX_DISTANCE             : ... (maximum allowed distance between a data source location and
                                the output grid location)
IWD_N                    : ... (inverse-distance weighting power; 2.0 = default)
```

```
Remarks:

- output field can only be obtained for one property/parameter per FillMatrix application run;

- although not required it is advisable that the attributed property/parameter name is
consistent with the MOHID naming convention, so that the produced output file is readily used
in other MOHID applications;

- FillMatrix application assumes that the input georeferenced data is in the same coordinate
system and property units requested for the output field;

- in OUTPUT_TIME keyword (case VARIABLE_IN_TIME : 1) at least two values in seconds must be
provided: the first value indicates the time in seconds after the start time when is made
first output, the last value is the afterwards time step of the outputs; when more than two
values are provided (maximum of 240 values) the values between the first and the last are
considered the time, from the start time, of specific outputs, after which the time step value
is considered for output: e.g. OUTPUT_TIME : 0 7200 21600 3600 means outputs at the start
time, after 7200 seconds from start, after 21600 seconds from start and after every 3600
seconds;

- the default value for the IWD_N keyword is considered by N.O.A.A. 1972, National Weather
Service River Forecast System forecast procedures, TM NWS HYDRO-14, US Department of
Commerce, Washington DC. adequate for interpolation of rainfall data.
```

## Sample

### With Interpolation

```
PROPERTY_NAME            : batim
PROPERTY_UNITS           : m
VARIABLE_IN_TIME         : 0
TIME_INSTANT             : 2005 1 1 0 0 0

OUTPUT_FILE              : PolygonBatimIWN.dat

GRID_DATA_FILE           : PolygonBatim.dat
XYZ_FILE                 : vigo.xyz

INTERPOLATION_METHOD     : 2
MAX_DISTANCE             : 0.005
IWD_N                    : 2
```

### With Interpolation Example Precipitation Stations

```
PROPERTY_NAME            : precipitation
PROPERTY_UNITS           : mm

OUTPUT_FILE              : Precipitation.hdf5

VARIABLE_IN_TIME         : 1
START                    : 2001 1 1 0 0 0
END                      : 2002 1 1 0 0 0
MAX_TIME_SPAN            : 86400.
OUTPUT_TIME              : 3600.
SKIP_NULLVALUES          : 0

GRID_DATA_FILE           : Batim.dat

INTERPOLATION_METHOD     : 1

FILL_OUTSIDE_POINTS      : 1
WRITE_TRIANGLES          : 1
TRIANGLES_FILE           : Triangles.xy

<begin_station>
NAME                     : 1
X                        : -9
Y                        : 40
VALUE_TYPE               : TIMESERIE
FILENAME                 : St1.srm
DATA_COLUMN              : 2
<end_station>

<begin_station>
NAME                     : 2
X                        : -9.1
Y                        : 40.1
VALUE_TYPE               : TIMESERIE
FILENAME                 : St2.srm
DATA_COLUMN              : 2
<end_station>

<begin_station>
NAME                     : 3
X                        : -9.2
Y                        : 40.2
VALUE_TYPE               : TIMESERIE
FILENAME                 : St3.srm
DATA_COLUMN              : 2
<end_station>
```

### Without Interpolation

```
PROPERTY_NAME            : leaf area index
PROPERTY_UNITS           : m2/m2

VARIABLE_IN_TIME         : 1
START                    : 2002 10  1 0 0 0
END                      : 2004 10  1 0 0 0
OUTPUT_TIME              : 0 86400

OUTPUT_FILE              : LAI.hdf5

GRID_DATA_FILE           : MDTTrancao200mSD.dat

INTERPOLATION_METHOD     : 0
FILL_ID_FILE             : VegetationID200m.dat

<begin_station>
NAME                    : Artificial Surfaces
ID                      : 100
VALUE_TYPE              : SINGLE_VALUE
VALUE                   : 1.2
<end_station>

<begin_station>
NAME                    : Agriculture
ID                      : 200
VALUE_TYPE              : TIMESERIE
FILENAME                : ..\..\GeneralData\Vegetation\Wheat.dat
DATA_COLUMN             : 2
<end_station>

<begin_station>
NAME                    : Forest
ID                      : 300
VALUE_TYPE              : SINGLE_VALUE
VALUE                   : 5.
<end_station>
```
