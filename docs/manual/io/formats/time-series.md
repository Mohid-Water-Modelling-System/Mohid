# Time Series

## Overview

Time series files are ASCII files used to store information varying in time. These files can be easily plotted with [Mohid Time Series Editor](../../tools/gui/mohid-time-series-editor.md) and also, if associated with a location, in [MOHID GIS](../../tools/gui/mohid-gis.md). These files can be created manually to serve as input for the model and they are created by Mohid numerical programs as outputs. To learn how to define time series outputs in MOHID click here.

## File organization

### Header section

At least 2 keywords must be present for a file to be recognised as a time series file.

```
SERIE_INITIAL_DATA - definition of the series initial data given in MOHID Time format

TIME_UNITS         - time units of the data series. It can be SECONDS, MINUTES, HOURS, DAYS or MONTHS

TIME_CYCLE         - check if the time serie is a cycle one. For example in the TIME_UNITS is HOURS a
                           24 hour cycle is assumed. In this case user needs to define 24 values one for
                           each hour of the day.
```

### Data section

The data section is defined inside a [block](../keywords/blocks.md) defined by the tags **<BeginTimeSerie>** and **<EndTimeSerie>**. Inside the block data is formatted into columns separated by blank spaces. The leftmost column defines the accumulated time, in the time units defined in the header section and starting from the time series initial date, also defined in the header section. The number of columns is limited by the default number of characters allowed in one single line, which is 256. This number can be shortened or increased via optional compilation of the source code.

In order to [Mohid Time Series Editor](../../tools/gui/mohid-time-series-editor.md) or [Mohid GIS](../../tools/gui/mohid-gis.md) read the time series files and associate a certain column with the name of the property (e.g. temperature, salinity, wind speed, etc), a header line must be present on the line just above the **<BeginTimeSerie>** tag. The number of properties names in this header line must be equal to the number of columns in the data section including the accumulated time column. The names must be separated by blank spaces, thus properties names with blank spaces are not allowed. For example:

```
wind direction is wrong
```

```
wind_direction is right
```

[Mohid Water](../../water/mohid-water.md) and the other numerical programs generate this header line automatically when outputting results into time series files.

### Residual section

Mohid results written in time series format always present a residual section which consists of the average value of the property during the period of the simulation. Note, that if the output time step of the time teries is higher than the time step of the model the residual value (or average value) will be different from the average of the values written in the time series file. This section is written by default and can be switched off. See [Module TimeSerie](../../water/modules/timeserie.md) for more information.

### Sample

This sample representation of a time series file shows

```
SERIE_INITIAL_DATA : 2006 8 1 12 0 0
TIME_UNITS         : MINUTES
```

```
time temperature salinity
<BeginTimeSerie>
0  10.2  35.4
10 10.3  35.2
20 10.3  35.3
30 10.1  35.5
40 10.0  35.6
..  ..    ..
<EndTimeSerie>
```

## Related links

- [How to define time series outputs in MOHID](../../tools/gui/how-to-create-mohid-timeseries-outputs.md)
