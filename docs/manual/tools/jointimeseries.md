# JoinTimeSeries

## Introduction

JoinTimeSeries is a MOHID support tool designed to join several [time series](../io/formats/time-series.md) files in a single file.

## Input file

It is needed to define a text file with the name JoinTimeSeries.dat where the input and output files names and location should be defined.

```
OUT_FILE                      : ... (path/name of output file with the joined time series data)
DT                            : ... (time interval between data, in seconds)
ADD_MODEL_DT                  : ... (time used for adjustement, in seconds)
WRITE_RESIDUALS               : 0/1 (residual appearance in time series: 0 = appears not, 1 =
                              appears)

<BeginTimeSeriesFile>
NAME                          : ... (path/name of the time series data file to be joined)
<EndTimeSeriesFile>
<BeginTimeSeriesFile>
....
<EndTimeSeriesFile>
```

## Sample

```
OUT_FILE                      : BoiaFunchalLevel2B_WQHDF.dat
DT                            : 86410
ADD_MODEL_DT                  : 10
WRITE_RESIDUALS               : 0

<BeginTimeSeriesFile>
NAME                          : .\Run1\BoiaFunchal.srw
<EndTimeSeriesFile>

<BeginTimeSeriesFile>
NAME                          : .\Run2\BoiaFunchal.srw
<EndTimeSeriesFile>
```
