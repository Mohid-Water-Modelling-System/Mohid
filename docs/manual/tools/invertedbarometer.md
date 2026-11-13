# InvertedBarometer

This application creates water level time series for tide gauge locations according with the **inverted (inverse) barometer effect**.

Running options for this application are specified in an [input file](invertedbarometer.md) whose path is indicated in a [nomfich.dat file](../io/nomfich-dat-file.md).

## Introduction

The **inverted (inverse) barometer effect** specifies that *a 1 mb decrease of atmospheric pressure below a reference pressure results in a 1.01 cm increase in water level* (e.g. Kantha, L., K. Whitmer, G. Born, 1994, "The Inverted Barometer Effect in Altimetry: A Study in the North Pacific").

This rule is generally considered to provide a rough measure of the water level change due to large scale circulation induced by water surface atmospheric pressure gradients.

The InvertedBarometer application applies this rule to obtain mean water level time series at tide gauge locations based on time and space variable atmospheric pressure data supplied in HDF5 format.

The time discretization of the time series is determined by the time discretization of the HDF5 data.

The calculation of water level is performed considering:

- the atmospheric pressure data for each time gauge location is considered the one available for the HDF5 grid cell where that tide gauge is geographically located; for this, the coordinate system in gauge files is assumed equal to coordinate system of pressure data in HDF5 files;

- the **reference tide gauge water level**, provided in tide gauge files, is considered equivalent to the **reference atmospheric pressure** assumed 101330.0 Pa in the application (following Dorandeu, J., P. Le Traon, 1999, "Effects of Global Mean Atmospheric Pressure Variations on Mean Sea Level Changes from TOPEX/Poseidon", Journal of Atmospheric and Oceanic Technology, Vol. 16, No. 9, pp. 1279-1283);

- the absolute water level at the tide gauge location in every HDF5 time instant is determined as:

```
(absolute tide gauge water level) = (reference tide gauge water level) +
(water level change correspondent to the difference between the reference and actual
atmospheric pressures)
```

The name of each time series file is automatically generated as the following concatenation:

```
(name of water level time series) =
(user defined root name specific for each gauge file) +
(name of the gauge location in gauge file) +
(number of order of appearance of the gauge location in the gauge file)
```

The application can also create new tide gauge files for the ready use of the water level time series as tide gauge mean level time series in tidal forcing of MOHID Water.

**Typical use:**

Allows the transformation of gauge files for use for tidal forcing in MOHID Water models, to take into account the time evolution of the reference level of each gauge according to large scale circulation induced by atmospheric pressure change.

**Data input requirements:**

One or more HDF5 files with surface atmospheric pressure data, with rank 2D and Pa units, for a spatial domain containing tide gauge locations. The data must be located in the HDF5 path "/Results/atmospheric pressure" (according with MOHID name convention). If several, the HDF5 files must have common grid and non-overlapping time periods.

Grid file, in MOHID format ([grid](../io/formats/grid.md) or [grid data](../io/formats/grid-data.md)), with the definition of HDF5 data grid.

One or more gauge files (tidal components files), in MOHID format, containing information of gauge locations where water level time series are to be produced.

**Ouput:**

[Time series](../io/formats/time-series.md) files of water level for each gauge location provided by the supplied gauge files. Besides water level values each time serie have a column for covered/uncovered information, as required by MOHID Water; every time data is assumed covered (valued 1.0).

In case there are gauge locations outside the provided HDF5 grid those locations' coordinates and gauge file are saved in a Noinformation.dat file in the directory provided in ROOT\_SRT in [nomfich.dat file](../io/nomfich-dat-file.md).

Optionally can also be produced new tide gauge files: a copy of the input tide gauge files with extra keywords connecting to the water level time series. These files can be readily used in MOHID Water as tidal components files with time variable gauge mean water level. For the gauge locations outside the HDF5 grid area the original gauge block is maintained in the new tide gauge file with no additional keywords.

## Input file

The name of the input file must be provided in the [nomfich.dat file](../io/nomfich-dat-file.md) in use.

```
(block for each HDF5 containing surface atmospheric pressure data; may be several,
order irrelevant)
<BeginHDF5File>
NAME                    : ... (path/name of HDF5 file containing atmospheric pressure data)
<EndHDF5File>

START_TIME              : ... (start time for water level time series: yyyy mm dd hh mm ss)
END_TIME                : ... (end time for water level time series: yyyy mm dd hh mm ss)

GRID_FILENAME           : ... (path/name of grid file with HDF5 files's grid)

MAX_BUFFER_SIZE         : ... (maximum buffer size, an integer value, for temporary time series
                               storage for each gauge location; 100000 = default)

(block for each gauge file containing locations for water level time series; may be several)
<BeginGaugeFile>
NAME                    : ... (path/name of gauge file)

TIMESERIE_NAME          : ... (intended root path/name of output time series without extension)
WRITE_MOHID_INPUT       : 0/1 (1 = create a new gauge file ready for MOHID Water use, 0 = not
                               to create; 0 = default)

(if WRITE_MOHID_INPUT : 1)
OUTPUTNAME              : ... (path/name with extension of gauge file to be created for MOHID
                               Water use)
<EndGaugeFile>
```

```
Remarks:

- time interval between pressure data values in HDF5 files may be variable: times series are
constructed maintaing these intervals within and between files;

- since no fixed DT is assumed for pressure data, no alerts are issued at InvertedBarometer
run time because of data lacking intervals between HDF5 data files;

- the path/name of time serie provided in the block for each gauge file is used as basis for
all time series paths/names originated from that gauge file; the name of each time serie is
completed in InvertedBarometer with the reference «GaugeName#», where GaugeName is the one
provided in the gauge block in gauge file and # is the order of the gauge block in gauge file.
```

## Samples

### Input gauge file

```
<begingauge>
NAME       :  Gauge
LONGITUDE  :     -6.0000     13.0000     50.3724
LATITUDE   :     43.0000     42.0000     22.8268
METRIC_X   :     -6.2307
METRIC_Y   :     43.7063
REF_LEVEL  :      2.0800
TIME_REF   :      0.0000
O1         :    0.068361    320.9841
K1         :    0.065776     74.1854
N2         :    0.265215     71.0081
M2         :    1.294223     91.5081
K2         :    0.122431    119.8621
S2         :    0.463176    124.2623
<endgauge>
<begingauge>
NAME       :  Gauge
LONGITUDE  :     -6.0000      4.0000     44.9939
LATITUDE   :     44.0000      8.0000     40.5185
METRIC_X   :     -6.0792
METRIC_Y   :     44.1446
REF_LEVEL  :      2.0800
TIME_REF   :      0.0000
O1         :    0.068554    320.9165
K1         :    0.065686     74.1463
N2         :    0.265796     71.2711
M2         :    1.296347     91.7997
K2         :    0.122623    120.1885
S2         :    0.464600    124.6422
<endgauge>
<begingauge>
NAME       :  Gauge
LONGITUDE  :     -6.0000     13.0000     50.3759
LATITUDE   :     44.0000     37.0000     53.5107
METRIC_X   :     -6.2307
METRIC_Y   :     44.6315
REF_LEVEL  :      2.0800
TIME_REF   :      0.0000
O1         :    0.068675    321.2289
K1         :    0.066201     74.5448
N2         :    0.265124     71.9428
M2         :    1.290508     92.5688
K2         :    0.122080    120.9709
S2         :    0.462477    125.5901
<endgauge>
```

### Input file

```
<BeginHDF5File>
NAME                    : 1971ERA1973msl.hdf5
<EndHDF5File>

GRID_FILENAME           : 1971ERA1973.grd

<BeginGaugeFile>
NAME                    : TidalComponents.dat
TIMESERIE_NAME          : K:\Simula\GeneralData\Tide\CP
WRITE_MOHID_INPUT       : 1
OUTPUTNAME              : K:\Simula\GeneralData\Tide\TidalComponentsWithMeanLevel.dat
<EndGaugeFile>

START_TIME              : 1972 10 1 0 0 0
END_TIME                : 1972 10 5 0 0 0
```

### Output water level time series

```
Time Serie created by InvertedBarometer.exe
Based on pressure HDF5 2D data on grid by 1971ERA1973.grd

GRID_LOCALIZATION_I     : 13
GRID_LOCALIZATION_J     : 36
SERIE_INITIAL_DATA      : 1972. 10.  1.  0.  0.  0.
TIME_UNITS              : SECONDS
LONGITUDE               :   -6.000000       13.00000       50.37240
LATITUDE                :    43.00000       42.00000       22.82680

     Seconds   YY  MM  DD  HH  MM       SS   water_level           open_point
<BeginTimeSerie>
        0.00 1972  10   1   0   0   0.0000   0.205845308304E+001   1.0
    21600.00 1972  10   1   6   0   0.0000   0.207204771042E+001   1.0
    43200.00 1972  10   1  12   0   0.0000   0.207506346703E+001   1.0
    64800.00 1972  10   1  18   0   0.0000   0.209167337418E+001   1.0
    86400.00 1972  10   2   0   0   0.0000   0.208186054230E+001   1.0
   108000.00 1972  10   2   6   0   0.0000   0.211031007767E+001   1.0
   129600.00 1972  10   2  12   0   0.0000   0.210468816757E+001   1.0
   151200.00 1972  10   2  18   0   0.0000   0.209870862961E+001   1.0
   172800.00 1972  10   3   0   0   0.0000   0.207271218300E+001   1.0
   194400.00 1972  10   3   6   0   0.0000   0.205773806572E+001   1.0
   216000.00 1972  10   3  12   0   0.0000   0.206882834435E+001   1.0
   237600.00 1972  10   3  18   0   0.0000   0.209010624886E+001   1.0
   259200.00 1972  10   4   0   0   0.0000   0.207881069183E+001   1.0
   280800.00 1972  10   4   6   0   0.0000   0.209054899216E+001   1.0
   302400.00 1972  10   4  12   0   0.0000   0.209281420708E+001   1.0
   324000.00 1972  10   4  18   0   0.0000   0.210097479820E+001   1.0
   345600.00 1972  10   5   0   0   0.0000   0.210286545753E+001   1.0
<EndTimeSerie>
```

### Output gauge file (for MOHID water use)

```
<begingauge>
NAME       :  Gauge
LONGITUDE  :     -6.0000     13.0000     50.3724
LATITUDE   :     43.0000     42.0000     22.8268
METRIC_X   :     -6.2307
METRIC_Y   :     43.7063
REF_LEVEL  :      2.0800
TIME_REF   :      0.0000
O1         :    0.068361    320.9841
K1         :    0.065776     74.1854
N2         :    0.265215     71.0081
M2         :    1.294223     91.5081
K2         :    0.122431    119.8621
S2         :    0.463176    124.2623
REF_EVOLUTION : Time Serie
TIME_SERIE_FILE : K:\Simula\GeneralData\Tide\CPGauge1.bts
REFLEVEL_COLUMN : 8
COVERED_COLUMN : 9
<endgauge>
<begingauge>
NAME       :  Gauge
LONGITUDE  :     -6.0000      4.0000     44.9939
LATITUDE   :     44.0000      8.0000     40.5185
METRIC_X   :     -6.0792
METRIC_Y   :     44.1446
REF_LEVEL  :      2.0800
TIME_REF   :      0.0000
O1         :    0.068554    320.9165
K1         :    0.065686     74.1463
N2         :    0.265796     71.2711
M2         :    1.296347     91.7997
K2         :    0.122623    120.1885
S2         :    0.464600    124.6422
REF_EVOLUTION : Time Serie
TIME_SERIE_FILE : K:\Simula\GeneralData\Tide\CPGauge2.bts
REFLEVEL_COLUMN : 8
COVERED_COLUMN : 9
<endgauge>
<begingauge>
NAME       :  Gauge
LONGITUDE  :     -6.0000     13.0000     50.3759
LATITUDE   :     44.0000     37.0000     53.5107
METRIC_X   :     -6.2307
METRIC_Y   :     44.6315
REF_LEVEL  :      2.0800
TIME_REF   :      0.0000
O1         :    0.068675    321.2289
K1         :    0.066201     74.5448
N2         :    0.265124     71.9428
M2         :    1.290508     92.5688
K2         :    0.122080    120.9709
S2         :    0.462477    125.5901
REF_EVOLUTION : Time Serie
TIME_SERIE_FILE : K:\Simula\GeneralData\Tide\CPGauge3.bts
REFLEVEL_COLUMN : 8
COVERED_COLUMN : 9
<endgauge>
```

## How to Include the Inverted Barometer Effect in MOHID Simulations

To impose the inverted barometer effect in MOHID model simulations you have regard to the **Hydrodynamic** and **Atmosphere** modules. To visualize the effect the spatial gradient of atmospheric pressure should be appreciable. It means that the effect is visible in large scale models

## Samples

### Hydrodynamic Module

```
LOCAL_SOLUTION            : 6 or 7 !6-GaugePlusSubModel and 7-AssimilaGaugeSubModel
ATM_PRESSURE              : 1 !activate the use of the property pressure
ATM_PRESSURE_TYPE         : 1 or 2 !OPTION 1 input data must be atmospheric pressure and OPTION 2 Mean Sea Level Pressure (MSLP)
ATM_RAMP                  : 1 !activate a smooth period
ATM_PERIOD                : 86400 !is the period by default in MOHID if the user not provide any one
```

```
IMPOSE_INVERTED_BAROMETER : 1 !activates the inverted barometer approximation in the open boundary
```

```
INVERTED_BAROMETER_CELLS  : 1 !allows restrain to specific cells in the open boundary where you want to impose the inverted barometer approximation
```

```
<begin_InvertBarometerCells>
 NAME                    : generic property
 DEFAULTVALUE            : 0
 TYPE_ZUV                : z
 FILE_IN_TIME            : NONE
 REMAIN_CONSTANT         : 1
 INITIALIZATION_METHOD   : ASCII_FILE
 FILENAME                : !Path and name of the file. The file is similar that that a grid data file with the specific i j cells in your domain
<end_InvertBarometerCells>
```

### Atmosphere Module

```
<beginproperty>
 NAME                    : atmospheric pressure or
 UNITS                   : Pa
 DESCRIPTION             : Meteo-IST MM5 atmospheric pressure
 DEFAULTVALUE            : 100000
 REMAIN_CONSTANT         : 0
 FILE_IN_TIME            : HDF
 FILENAME                : !path and name of the file
 OUTPUT_HDF              : 1
<endproperty>
```

```
<beginproperty>
 NAME                    : mean sea level pressure
 UNITS                   : Pa
 DESCRIPTION             : mean sea level pressure GFS
 DEFAULTVALUE            : 0.
 FILE_IN_TIME            : HDF
 FILENAME                : !Path and name of the file
 OUTPUT_HDF              : 1
 TIME_SERIE              : 1
<endproperty>
```
