# ConvertHDF5ToGridData

The application **ConvertHDF5ToGridData** converts an HDF5 file into one or more grid data files (it extracts or one instant or all the intants for the defined properties).

Running options for this application are specified in an [input file](converthdf5togriddata.md) called HDF5ToGridData.dat.

## Introduction

With this tool the user provides an HDF and specifies a given instant (integer number) or all the instants will be extracted. And the user specifies also the properties to convert (one grid data per property and per instant will be produced).

**Typical use:**

To extract grids from HDF (to allow interpolations that need the grid) or to get HDF snapshots in grid data format for specific properties.

**Data input requirements:**

An HDF5 file containing data to be converted.

**Ouput:**

One or all of the HDF instants in grid data format. One file per instant and per property.

## Input file

The input file is called HDF5ToGridData.dat.

```
HDF_FILE                 : ... (path/name of HDF5 file with data to extract)

ALL_INSTANTS             : o/1 (wether or not all instants to extraxt)

INSTANT                  : int (if ALL_INSTANTS : 0 then checks what HDF instant
                           to convert (integer number between one and maximum number
                           of instants in HDF)

(block for each property/parameter which is object of extraction, may be several)
<BeginParameter>
PROPERTY                 : ... (parameter name, according with MOHID V4)
HDF_GROUP                : ... (complete path in HDF5 file to parameter data, according with
<EndParameter>                  MOHID V4)
```

## Sample

```
HDF_FILE          : WRF_d05_1km.hdf5

ALL_INSTANTS      : 0
INSTANT           : 1

<BeginParameter>
PROPERTY                : air temperature
HDF_GROUP               : /Results/air temperature
<EndParameter>

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
