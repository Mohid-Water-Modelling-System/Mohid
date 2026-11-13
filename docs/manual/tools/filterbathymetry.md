# FilterBathymetry

The application **FilterBathymetry** performs filtering of MOHID bathymetries.

Running options for this application are specified in an input file named [Filter.dat](filterbathymetry.md).

## Introduction

Two options exist for filtering operation: - calculate the new bathymetry in each cell so as to ponder the average bathymetry of cells located in a specified area centered in that cell (defined with radius valued in number of cells) and the original cell bathymetry; - calculate the new bathymetry in each cell so that slope of the bathymetry verified at that cell is no greater than a specified value.

If the first option is chosen then the user is required to supply the intended radius (number of cells) and weighting factor between the neighborhood average and the original bathymetry values. Optionally, the user can instruct the application only to filter bathymetry values which are originally between user defined minimum and maximum bathymetry values, allowing the filtering process to alter only parts of the original bathymetry.

In the case of the second option, the user should supply the maximum slope.

**Typical use:**

This application is usually used to alter a bathymetry (previously created in the MOHID framework) which, because of its irregularity, does not leads to satisfactory MOHID Water model results.

**Data input requirements:**

The bathymetry file (in MOHID grid data format) to be filtered.

**Ouput:**

Filtered bathymetry.

## Input file ([Filter.dat](../io/filter-dat.md))

```
IN_BATIM          : ... (path/name of bathymetry file to be filtered)
OUT_BATIM         : ... (path/name of output filtered bathymetry file)

SIGMA_SMOOTH      : 0/1 (1 = automatic smoothing in a way that the slope of bathymetry in each
                         grid cell is below a slope parameter, 0 = do not use this option
                        (default))

(if SIGMA_SMOOTH : 0:)
FILTER_RADIUS     : ... (number of influence cells to filter bathymetry)
FACTOR            : ... (filtering factor)
H_MIN             : ... (minimum original bathymetry value for filtering to occur; default =
                         -1x109)
H_MAX             : ... (maximum original bathymetry value for filtering to occur; default =
                         1x109)

(if SIGMA_SMOOTH : 1:)
SLOPE_PARAMETER   : ... (slope parameter value; default = 0.2)
```

```
Remarks:

- the value indicated by keyword SLOPE_PARAMETER is considered in the literature to be
  adequate for sigma vertical resolution models if not larger than 0.2;
```

## Sample

```
IN_BATIM        : U:\OndasInternas3\GeneralData\Batim\Bat_Setubal_Nazare_Level2.dat
OUT_BATIM       : U:\OndasInternas3\GeneralData\Batim\Bat_Setubal_Nazare_Level2Filtered.dat

SIGMA_SMOOTH    : 1
SLOPE_PARAMETER : 0.2

H_MIN           : 1000
H_MAX           : 5000
```
