# ConvertToHDF5 WOA2005

The World Ocean Atlas (WOA) 2005 climatology provides results for water temperature, salinity and several water quality and biology properties. The data is provided by [NOAA NODC](http://www.nodc.noaa.gov).

## Introduction

The ConvertToHDF5 action converts the climatological data for the following properties:

- Dissolved Oxygen

- Phosphate

- Nitrate

- Silicate

For temperature and salinity fields, refer to [ConverToHDF5 LEVITUS FORMAT](converttohdf5.md) because these fields have a lower resolution(1º) than Levitus(0.25º).

### Spacial Coverage

The resulting [HDF5](hdf5.md) files have a spatial resolution of 1 degree latitute and longitude. The resulting grid spans from -180º to 180º rather than the original 0º to 360º format. Thus, the lower left corner of the hdf5 grid is -180º longitude, 90ºS longitude.

### Units

The data in netcdf format provided by [NOAA](http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html) is in uM for the nutrients and ml/l for dissolved oxygen. There is an option to convert to MOHID units, namelly mg/L for nutrients and mg/l for oxygen.

### Kepler data folder

For those users working in MARETEC/HIDROMOD network, the following data folder contains the converted files, sample ConvertToHDF5 files(with WOA geometry and grid) and WOA05 documentation:

```
\\kepler\DataCenter\DadosBase\Oceanografia\World Ocean Atlas 2005\
```

## Input file (ConvertToHDF5Action.dat)

```
<begin_file>
  ACTION                      : char    -    'CONVERT WOA FORMAT' to use this module
  FILEPATH                    : char    -    Path to WOA05 original netcdf files(n0112an1.nc for e.g.)
                                                               with ending slash
  OUTPUTFILENAME              : char    -    Path to WOA05 file generated from NETCDF file
  MOHID_UNITS                  : 0/1    -    Convert uunits to mg/l
  <<begin_woafiles>>
  (netcdf file)                 char    -    WOA05 netcdf file in path defined in FILEPATH, one per line
  (netcdf file2)                char         e.g.: n0112an1.nc i0112an1.nc p0112an1.nc o0112an1.nc
   ...
  <<end_woafiles>>
  CONVERT_TO_ASCII            : 0       -    Flag to convert to ascii file: not implemented
  CONVERT_TO_HDF5             : 1       -    Flag to convert to hdf5 file
  XX_VARIABLE                 : char    -    'longitude'
  YY_VARIABLE                 : char    -    'latitude'
  ZZ_VARIABLE                 : char    -    'depth'
  TIME_VARIABLE               : char    -    'time'
<end_file>
```

```
Remarks:

- only CONVERT_TO_HDF5 : 1 is implemented
- files in block begin/end_woafiles are all concatenated into one hdf5 file
- FILEPATH must have an ending slash
```

## Sample

All sample files are named *ConvertToHDF5Action.dat*.

### Convert

```
<begin_file>
ACTION                    : CONVERT WOA FORMAT

!WOA05_o2sat.hdf5 WOA05_nut_uM.hdf5 WOA05_nut_mg_l.hdf5
OUTPUTFILENAME            : WOA05_nut_uM.hdf5
MOHID_UNITS               : 0

!filepath with ending slash!
FILEPATH                  : Z:\netcdf\monthly\

<<begin_woafiles>>
n0112an1.nc
i0112an1.nc
p0112an1.nc
<<end_woafiles>>

CONVERT_TO_ASCII          : 0
CONVERT_TO_HDF5           : 1

XX_VARIABLE               : lon
YY_VARIABLE               : lat
ZZ_VARIABLE               : depth
TIME_VARIABLE             : time
<end_file>
```

```
Remarks
The converted files(O2 Saturation and Nutrients in uM and mg/l) are in kepler folder:
\\kepler\DataCenter\DadosBase\Oceanografia\World Ocean Atlas 2005\Data\data\
Since the conversion only is done one time for the whole grid, there is no need to convert every time.
```

### Extrapolate

Extrapolate data for all regions for the interpolation(triagulation) to work fine.

```
<begin_file>
ACTION                   : INTERPOLATE GRIDS

TYPE_OF_INTERPOLATION    : 1

FATHER_FILENAME           : \\kepler\DataCenter\DadosBase\Oceanografia\World Ocean Atlas 2005\Data\data\WOA05_nut_mg_l.hdf5
FATHER_GRID_FILENAME      : ..\WOA05_Grid_3D.dat

OUTPUTFILENAME           : ..\WOA05_nut_mg_l_alldata.hdf5
NEW_GRID_FILENAME        : ..\WOA05_Grid_3D.dat

START                    : -9999 1  1 0 0 0
END                      : -9999 12 1 0 0 0

INTERPOLATION3D          : 1
FATHER_GEOMETRY          : ..\WOA05_Geometry_3D_nuts.dat
NEW_GEOMETRY             : ..\WOA05_Geometry_3D_nuts.dat
AUX_GRID_FILENAME        : ..\WOA05_Grid_3D.dat
AUX_OUTPUTFILENAME       : ..\AuxAllPointsWithData.hdf5

POLI_DEGREE              : 1
DO_NOT_BELIEVE_MAP       : 1
EXTRAPOLATE_2D           : 2
<end_file>
```

### Interpolate

Finally, interpolate to the final grid and geometry (same as the Interpolate 3D sample):

```
<begin_file>
ACTION                   : INTERPOLATE GRIDS

TYPE_OF_INTERPOLATION    : 1
FATHER_FILENAME          :  ..\WOA05_nut_mg_l_alldata.hdf5
OUTPUTFILENAME           : ..\NS_Nuts_3D.hdf5
FATHER_GRID_FILENAME     : ..\WOA05_Grid_3D.dat
NEW_GRID_FILENAME        : ..\NorthSea_submodel__.new

START                    : -9999 1  1 0 0 0
END                      : -9999 12 1 0 0 0

INTERPOLATION3D          : 1
FATHER_GEOMETRY          : ..\WOA05_Geometry_3D_nuts.dat
NEW_GEOMETRY             : ..\GeometryModel3D.dat
AUX_OUTPUTFILENAME       : ..\AuxCadizMonthlyLevitus.hdf5
AUX_GRID_FILENAME        : ..\NorthSea_submodel__.new

POLI_DEGREE              : 1
DO_NOT_BELIEVE_MAP       : 1
<end_file>
```
