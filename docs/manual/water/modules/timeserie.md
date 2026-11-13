# Module TimeSerie

## Overview

Module TimeSerie manages all operations involving [Time Series](../../io/formats/time-series.md) files, namely in terms of input and output.

## Time Series input

## Time Series output

Module TimeSerie handles results output in time series format. It acts as a server for other modules that which to output results in this format. This the client module needs to define the specifications with the location and other options of the output time series. This is made via a file:

```
TIME_SERIES_LOCATION     : ..\..\GeneralData\TimeSeriesFile.dat
```

or reading directly from the client modules' input data file.

The Time Series output configuration is done via the follow keywords.

```
MAX_BUFFER_SIZE  : 100000 - The maximum BufferSize is set here to 0.1Mb (for each property)

COMPUTE_RESIDUAL : 1 - Computes the time serie average value
```

```
REPLACE_PATH : ****  - the user can define a specific path for time serie file
                       instead of the default one ( ..\res).

IGNORE_ON : 0 - logical that if it is true (1) the module do not stop the program if
                  the time serie location is outside of the model grid

FIRST_OUTPUT_TIME :    - By default the first output instant is the starting time of the
                        model. However the user can specify with this keyword an
                        alternative starting output time.

DT_OUTPUT_TIME :       - Time series time step by default equal to the model time step.
```

For each time serie output the user need to define a block with the time serie location

```
<BeginTimeSerie>
NAME             :  - Name of the time serie
DEPTH_LEVEL      :  - depth of the time serie location
COORD_X          :  - x coordinate of the time serie location
COORD_Y          :  - y coordinate of the time serie location
<EndTimeSerie>
```

The time serie location can also be given specifying the grid cells index (i,j,k)

```
<BeginTimeSerie>
LOCALIZATION_I   :  - grid cell line location
LOCALIZATION_J   :  - grid cell column location
LOCALIZATION_K   :  - grid cell layer location
<EndTimeSerie>
```

The [time series](../../io/formats/time-series.md) output format is described in [http://wiki.mohid.com/wiki/index.php?title=Time\_Series](../../io/formats/time-series.md)

Example of definition of the time series output

```
DT_OUTPUT_TIME                :  60
MAX_BUFFER_SIZE               :  10000

<BeginTimeSerie>
NAME                          : !Malaysia
COORD_X                       :  100.6
COORD_Y                       :  4.2
DEPTH_LEVEL                   : 0
<EndTimeSerie>

<BeginTimeSerie>
NAME                          : !KualaTerengganu
COORD_X                       :  103.1833
COORD_Y                       :  5.26666
<EndTimeSerie>
```
