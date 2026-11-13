# Module Statistics

**Module Statistics** calculates the following time statistics:

- **arithmetic average**;

- **geometric average** (optional);

- **sample standard deviation from arithmetic average**;

- **sample standard deviation from geometric average** (optional);

- **maximum value**;

- **minimum value**;

- **accumulated value** (optional).

- **percentage of run-time below a critical value** (optional).

It is better described in the support tool [HDF5Statistics](../../tools/hdf5statistics.md), but can be used in the [MohidWater](../../reference/mohidwater.md) program.

## Input file

The statistics configuration file can be configured using the following keywords:

```
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

CRITICAL                : 0/1 (1 = calculate percentage of time below a critical value; 0 =
                               default)
CRITICAL_VALUE          : ... (parameter critical value, real; default = 0.02)

GEOMETRIC_MEAN          : 0/1 (0 = do not calculate geometric mean for non negative parameters,
                               1 = calculate geometric mean for non negative parameters;
                               0 = default)
ACCUMULATED             : 0/1 (0 = do not calculate accumulated values, 1 = calculate
                               accumulated values; 0 = default)
```

### Hydrodynamic

Assign the following keywords to the hydrodynamic configuration file:

```
STATISTICS               : 1
STATISTICS_FILE          : ..\..\GeneralData\Statistics\stats.dat

STATISTICS_2D            : 1
STATISTICS_FILE_2D       : ..\..\GeneralData\Statistics\stats2D.dat
```

Note that there are two statistics files. One for the *u*, *v* and *w* components (3D fields). And another for the *waterlevel* (2D field). Make sure you define the correct statistical method in each case.

### WaterProperties

Under each property assign the following keywords:

```
<beginproperty>
...
STATISTICS               : 1
STATISTICS_FILE          : ..\..\GeneralData\Statistics\stats.dat
...
<endproperty>
```
