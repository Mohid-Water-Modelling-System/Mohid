# HDF5Statistics

The application **HDF5Statistics** produces time statistics in HDF5 format from data existent in HDF5 files.

Running options for this application are specified in an [input file](hdf5statistics.md) whose path is indicated in a [nomfich.dat file](../io/nomfich-dat-file.md).

## Introduction

HDF5Statistics uses [module Statistics](../water/modules/statistics.md) which calculates the following time statistics:

- **arithmetic average**;

- **geometric average** (optional);

- **sample standard deviation from arithmetic average**;

- **sample standard deviation from geometric average** (optional);

- **maximum value**;

- **minimum value**;

- **accumulated value** (optional).

- **percentage of run-time below a critical value** (optional).

Time statistics can be calculated from a sample of data defined by a spatial unit and a time period. In HDF5Statistics application the **spatial unit is each 2D or 3D cell** and the time period is defined by the user, through the indication of the limit instants of a **time window**.

Several options exist for way time statistics are calculated:

- **global statistics**: statistics are produced, for each spatial unit, considering all data contained in the user defined time window;

- **daily statistics**: statistics are produced, for each spatial unit, considering data for each calendar day contained in the user defined time window, i.e., each day has its own statistics;

- **monthly statistics**: statistics are produced, for each spatial unit, considering data for each calendar month contained in the user defined time window, i.e., each month has its own statistics;

- **specific hour statistics**: statistics are produced, for each spatial unit, considering all data contained in the user defined time window refering to an user defined hour; only one hour can be defined for each HDF5Statistics run.

**Typical use:**

Obtain statistics for data contained in a set of MOHID Water result files or from an HDF5 file resulting from [glue](converttohdf5.md).

**Data input requirements:**

One or more HDF5 files with time dependent non "Grid" data for which statistics are to be calculated and equal time independent "Grid" data. If several files are provided they should contain sequential and same parameters data for the sake of statistics consistency, although data lacking periods are supported by the application.

**Ouput:**

One HDF5 file in MOHID format containing:

- a "Grid" group with time independent data of the original files;

- a "Statistics" group containing the statistics data for every parameter requested; the statistic data for each parameter is organized first by the type of statistical time period - "Global", "Daily", "Monthly", "SpecificHour" - and within each period by the type of statistical parameter - "Accumulated", "Average", "CriticalValue", "GeomAverage", "GeomStandDev", "Maximum", "Minimum", "StandDev";

- a "Time" group containing the first and last time instants of data considered for statistical calculation; these instants are the limits of the actual time window considered for calculation, defined by the time window specified by the user and the available time instants in the input HDF5 files.

## Input file

The name of the input file must be provided in the [nomfich.dat file](../io/nomfich-dat-file.md) in use.

```
(block for each HDF5 containing data; may be several, order of blocks is irrelevant)
<BeginHDF5File>
NAME                    : ... (path/name of HDF5 file with data to calculate statistics from)
<EndHDF5File>

START_TIME              : ... (start date for time window: yyyy mm dd hh mm ss)
END_TIME                : ... (end date for time window: yyyy mm dd hh mm ss)

METHOD_STATISTIC        : ... (statistical method: 1 = Values3DStats3D, 3 = Values2DStats2D;
                               1 = default)

GLOBAL_STATISTIC        : 0/1 (1 = calculate statistics for the whole time window period;
                               0 = default)
DAILY_STATISTIC         : 0/1 (1 = calculate daily statistics; 0 = default)
MONTHLY_STATISTIC       : 0/1 (1 = calculate monthly statistics; 0 = default)
SPECIFIC_HOUR_STATISTIC : 0/1 (1 = calculate statistics for data referring to a specific hour
                               of the day; 0 = default)

SPECIFIC_HOUR           : ... (if SPECIFIC_HOUR_STATISTIC : 1, hour of day for statistics
                               calculation; 12 = default)

(block for each parameter that is object of statistics; may be several)
<BeginParameter>
PROPERTY                : ... (parameter/property name)
HDF_GROUP               : ... (complete path to parameter data in HDF5 file )
<EndParameter>

OUTPUTFILENAME          : ... (path/name of output HDF5 file with statistics)

3D_HDF5                 : 0/1 (1 = HDF5 files are 3D, 0 = HDF5 files are 2D; 0 = default)
HDF5_MAP_ITEM           : ... (map item name in HDF5 file)

GEOMETRIC_MEAN          : 0/1 (0 = do not calculate geometric mean for non negative parameters,
                               1 = calculate geometric mean for non negative parameters;
                               0 = default)
ACCUMULATED             : 0/1 (0 = do not calculate accumulated values, 1 = calculate
                               accumulated values; 0 = default)
CRITICAL                : 0/1 (1 = calculate percentage of time below a critical value; 0 =
                               default)
CRITICAL_VALUE          : ... (parameter critical value, real; default = 0.02)
```

```
Remarks:

- the time window defined by START_TIME and END_TIME keywords has to contain more than one time
instant;

- the actual time window used for statistic calculation is defined by the START_TIME and
END_TIME keywords and the data availability in the supplied input HDF5 files;

- statistical method 1 (Values3DStats3D) calculates statistics for 3D parameters for the whole
3D field (one statistic per cell and layer);

- statistical method 2 (Values2DStats2D) calculates statistics for 2D parameters for the 2D
field (one statistic per cell);

- at least one type of statistics (GLOBAL_STATISTIC, DAILY_STATISTIC, MONTHLY_STATISTIC or
SPECIFIC_HOUR_STATISTIC) has to be chosen for any statistics file to be produced;

- when performing daily/monthly statistics it is not required that day/month is complete,
day/month is defined by calendar day/month;

- the 0 hour data is accounted for statistics of the previous and next day;

- the first 0 hour data of a month (first day) is accounted for statistics in previous
and next month;

- in GEOMETRIC_MEAN option null values are considered unity because of calculation needs; this
procedure is common for the calculation of the geometric mean of coliforms but may not be
adequate for other parameters;

- ACCUMULATED option is useful for precipitation statistics;

- CRITICAL option is useful for estimating the time percentage of the run where the velocity
modulus was below a critical velocity value; however the CRITICAL option may be used for any
other parameter;

- statistics for both 2D and 3D parameters cannot be calculated in same run.
```

## Sample

```
<BeginHDF5File>
NAME                    : K:\HDF5Statistics\Lagrangian_1.hdf5
<EndHDF5File>

<BeginHDF5File>
NAME                    : K:\HDF5Statistics\Lagrangian_2.hdf5
<EndHDF5File>

<BeginHDF5File>
NAME                    : K:\HDF5Statistics\Lagrangian_3.hdf5
<EndHDF5File>

START_TIME              : 2005 5 1 0 0 0
END_TIME                : 2005 5 10 0 0 0

METHOD_STATISTIC        : 1

GLOBAL_STATISTIC        : 1
DAILY_STATISTIC         : 1
MONTHLY_STATISTIC       : 1

SPECIFIC_HOUR_STATISTIC : 1
SPECIFIC_HOUR           : 12

ACCUMULATED             : 1
GEOMETRIC_MEAN          : 1
CRITICAL                : 0
CRITICAL_VALUE          : 0.02

<BeginParameter>
HDF_GROUP               : /Results/Group_1/Data_3D/fecal coliforms
PROPERTY                : fecal coliforms
<EndParameter>

OUTPUTFILENAME          : teste_lagrangian5.hdf5
3D_HDF5                 : 1
HDF5_MAP_ITEM           : WaterPoints3D
```
