# WOAnetcdfTOhdf5

This program converts the NetCDF files of the [World Ocean Atlas 2005](http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html) (WOA05) into HDF5, it allows to create HDF5 files for Mohid using the Module LIFE (mmol/m3) or the Module WaterQuality (mg/l). The user will be prompted when the software is ran.

Running options for this application are specified by the user in a input file named PARAMETER\_READ.dat.

## Input file (PARAMETER\_READ.dat)

```
ACTION                   : CONVERT LEVITUS FORMAT
OUTPUT_FILENAME          : ... (path/name of HDF5 file to be converted)
OUTPUT_GRID_FILENAME     : ... (path/name of grid data with horizontal discretization to be
                                created)
OUTPUT_GEOMETRY_FILENAME : ... (path/name of geometry file with vertical discretization to be
                                created)
<<begin_input_files>>
... (path/name of the WOA05 NetCDF data files, one per line, can be several)
...
<<end_input_files>>
```

## Samples

```
<begin_file>
ACTION                   : CONVERT LEVITUS FORMAT

OUTPUTFILENAME           : WOA05_climatology.hdf5
OUTPUT_GRID_FILENAME     : WOA05_Grid1.dat
OUTPUT_GEOMETRY_FILENAME : WOA05_Geometry1.dat

<<begin_input_files>>
t0112an1.nc
s0112an1.nc
n0112an1.nc
i0112an1.nc
p0112an1.nc
o0112an1.nc
<<end_input_files>>
<end_file>
```
