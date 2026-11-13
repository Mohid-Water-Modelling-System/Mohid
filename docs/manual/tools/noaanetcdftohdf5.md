# NOAAnetcdfTOhdf5

The NOAAnetcdfTOhdf5 is an application which allows the conversion of NOAA meteorological reanalysis files in NetCDF format to HDF5 format readable by MOHID applications.

Running options for this application are specified by the user in a input file named files.dat.

## Input file (files.dat)

```
INPUT_FILENAME    : ... (path/name of NetCDF file to be converted)
OUTPUT_FILENAME   : ... (path/name of HDF5 file to be converted)

JUST_REPORT       : 0/1 (1 = obtain to file only input file header info, 0 = default)
CONVERT_TO_ASCII  : 0/1 (1 = convert data to ASCII format, 0 = default)
CONVERT_TO_HDF5   : 0/1 (1 = convert data to HDF5 format, 0 = default)

XX_VARIABLE       : ... (name of longitude variable in the input file: usual name is "lon")
YY_VARIABLE       : ... (name of latitude variable in the input file: usual name is "lat")
TIME_VARIABLE     : ... (name of time variable in the input file: usual name is "time")
```

```
Remarks
- if JUST_REPORT : 1 no HDF5 file will be produced regardless of the choice in CONVERT_TO_HDF5;
- CONVERT_TO_ASCII : 1 should not usually be attempted because ASCII file can be potentially
  very large; its main use is for the application development and debugging;
- if CONVERT_TO_HDF5 : 0 or if this keyword is inexistent in the input file no HDF5 file will
  be created;
- the names in XX_VARIABLE, YY_VARIABLE and TIME_VARIABLE are expected to be the same for
  different NOAA reanalysis NetCDF files, i.e. "lon", "lat" and "time".
```

## Samples

```
INPUT_FILENAME    : D:\Aplica\NOAAreanalysis\vwnd.sig995.2004.nc
OUTPUT_FILENAME   : D:\Aplica\NOAAWinds\Output\vwnd.sig995.2004
JUST_REPORT       : 0
CONVERT_TO_ASCII  : 0
CONVERT_TO_HDF5   : 1

XX_VARIABLE       : lon
YY_VARIABLE       : lat
TIME_VARIABLE     : time
```
