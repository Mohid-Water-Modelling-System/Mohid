# HDF5Extractor

The application **HDF5Extractor** extracts from an HDF5 file data to produce one or more smaller HDF5 files.

Running options for this application are specified in an [input file](hdf5extractor.md) whose path is indicated in a [nomfich.dat file](../io/nomfich-dat-file.md).

## Introduction

The extraction operations which can be accomplished with HDF5Extractor are the following:

- extraction of data within a time window of the time span of the original HDF5 file;

- split the data in several HDF5 files (extraction by interval), each containing the data refering to a time interval;

- extraction of data using a time frequency lower than original (skipping instants);

- extraction of data within a 2D spatial window of the spatial span of the original HDF5 file;

- extraction of data from only some vertical layers when the original HDF5 file is 3D;

- extraction of data for only some results properties.

These operations can be ordered by the user in a isolated way or at the same time.

The extraction by interval is not accomplished if the requested interval is longer than the provided HDF5 time period or if the provided HDF5 file only has one time.

**Typical use:**

Operations of data compression on existing HDF5 files to remove useless data, namely to accomplish storage space economy or to prepare [glue](converttohdf5.md) operations.

**Data input requirements:**

An HDF5 file containing data to be extracted.

**Ouput:**

One HDF5 file with the extracted data, except in the case the operation of extraction of the data by time intervals is chosen.

When extraction by interval is performed it is produced one file per every time period of length equal to interval specified by the user, begining in the original HDF5 file start time, and is produced also an HDF5 file with the remainder of the spliting, so that all data from input file is extracted.

## Input file

The name of the input file must be provided in the [nomfich.dat file](../io/nomfich-dat-file.md) in use.

```
FILENAME                 : ... (path/name of HDF5 file with data to extract)

OUTPUTFILENAME           : ... (path/name of output file (HDF5) with extracted data)

START_TIME               : ... (start date for time window: yyyy mm dd hh mm ss)
END_TIME                 : ... (end date for time window: yyyy mm dd hh mm ss)

XY_WINDOW_OUTPUT         : 0/1 (0 = extract all cells; 1 = extract only cells in spatial
                                window specified in XY_WINDOW_LIMITS; 0 = default)
XY_WINDOW_LIMITS         : ... ... ... ... (window limits: ILB (i of lower-left cell), JLB
                                            (j of lower-left cell), IUB (i of upper-right
                                            cell), JUB (j of upper-right cell))

LAYERS_OUTPUT            : 0/1 (0 = extract all layers, 1 = extract only layers specified in
                                LAYERS_MIN_MAX; 0 = default)
LAYERS_MIN_MAX           : ... ... (KLB (lower layer), KUB (upper layer))

INTERVAL                 : 0/1 (0 = time extraction of period indicated in START_TIME
                                and END_TIME; 1 = time extraction of every time period
                                in file with length provided by DT_INTERVAL; 0 = default)
DT_INTERVAL              : ... (time interval for extraction, in seconds)

SKIP_INSTANTS            : 0/1 (0 = use all time instants from original HDF;
                                1 = extract time instants at a defined frequency given by
                                 DT_SKIP_INSTANTS (use lower frequency than original))
DT_SKIP_INSTANTS         : ... (time interval to extract instants, in seconds)

CONVERT_V3_TO_V4         : 0/1 (0 = do not convert file format from MOHID V3 to MOHID V4;
                                1 = convert file format from MOHID V3 to MOHID V4; 0 = default)

(block for each property/parameter which is object of extraction, may be several)
<BeginParameter>
PROPERTY                 : ... (parameter name, according with MOHID V4)
HDF_GROUP                : ... (complete path in HDF5 file to parameter data, according with
                                MOHID V4)

(if CONVERT_V3_TO_V4 : 1:)
PROPERTY_V3              : ... (parameter name in MOHID V3)
HDF_GROUP_V3             : ... (complete path in HDF5 file to parameter data in MOHID V3)
<EndParameter>
```

```
Remarks:

- in the spatial window specification i, j are the cell indentification indexes in the MOHID
framework;

- in layer specification "lower" means the layer closest to bottom, which has smaller layer
number (k), and "upper" means the layer closest to surface, which has a higher layer number
(k));

- extraction by interval option (INTERVAL : 1) produces one file per every time
period of length equal to interval specified by the user, begining in the file's start
time: e.g. if DT_INTERVAL : 21600 (6 hours interval) then a HDF5 file is produced with the
data of each 21600 seconds period and also with the remainder of the division, so that all
data from provided file is extracted;

- in case of extraction by interval the time window eventually defined by START_TIME and
END_TIME keywords is not considered;

- extraction by interval is not accomplished if interval is longer than the provided HDF5
time period (or if the provided HDF5 file only has one time);

- extraction by skipping instants is done by selecting the first HDF instant and then
selecting instants at a defined interval given by user in DT_SKIP_INSTANTS (e.g. selecting
daily results (86400s) in hourly HDF).

- exctraction by skipping instants can be used in conjunction with extaction by interval
where all the created output HDF files will have the frequency defined in DT_SKIP_INSTANTS

- only time dependent parameters can be extracted, because extraction is meaningless
otherwise: for instance "Residual" and "Statistics" group parameters cannot be extracted;

- as parameters to specify in parameter block are meant result properties (but not only those
present in the "Result" group): "Grid" time dependent variables are not to be specified but
are automatically extracted;

- at least one parameter block must be present in input file for any extraction to be
accomplished.
```

## Sample

### Example with area extraction and interval extraction

```
OUTPUTFILENAME           : compactedhydro.hdf5

FILENAME                 : Hydrodynamic_15.hdf5

START_TIME               : 2005 2 28 0 0 0
END_TIME                 : 2005 2 28 18 22 0

XY_WINDOW_OUTPUT         : 1
XY_WINDOW_LIMITS         : 1 1 70 100

LAYERS_OUTPUT            : 1
LAYERS_MIN_MAX           : 10 11

INTERVAL                 : 1
DT_INTERVAL              : 21600

CONVERT_V3_TO_V4         : 0

<BeginParameter>
PROPERTY                 : velocity modulus
HDF_GROUP                : /Results/velocity modulus
<EndParameter>

<BeginParameter>
PROPERTY                 : velocity U
HDF_GROUP                : /Results/velocity U
<EndParameter>

<BeginParameter>
PROPERTY                 : velocity V
HDF_GROUP                : /Results/velocity V
<EndParameter>

<BeginParameter>
PROPERTY                 : water level
HDF_GROUP                : /Results/water level
<EndParameter>
```

### Example with area extraction and time skipping extraction

```
OUTPUTFILENAME           : compactedhydro.hdf5

FILENAME                 : Hydrodynamic_15.hdf5

START_TIME               : 2005 2 28 0 0 0
END_TIME                 : 2005 2 28 18 22 0

XY_WINDOW_OUTPUT         : 1
XY_WINDOW_LIMITS         : 1 1 70 100

LAYERS_OUTPUT            : 1
LAYERS_MIN_MAX           : 10 11

INTERVAL                 : 0

SKIP_INSTANTS            : 1
DT_SKIP_INSTANTS         : 86400

CONVERT_V3_TO_V4         : 0

<BeginParameter>
PROPERTY                 : velocity modulus
HDF_GROUP                : /Results/velocity modulus
<EndParameter>

<BeginParameter>
PROPERTY                 : velocity U
HDF_GROUP                : /Results/velocity U
<EndParameter>

<BeginParameter>
PROPERTY                 : velocity V
HDF_GROUP                : /Results/velocity V
<EndParameter>

<BeginParameter>
PROPERTY                 : water level
HDF_GROUP                : /Results/water level
<EndParameter>
```
