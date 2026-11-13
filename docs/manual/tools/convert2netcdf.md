# Convert2netcdf

## Overview

**Convert2netcdf** is a tool to convert MOHID results files written in [HDF5](hdf5.md) format into [netcdf](../archive/linux/netcdf.md). This tool makes use of [ModuleNETCDF](../water/modules/modulenetcdf.md) (under development and to be included in [MOHID Base 1](../water/README.md#mohid-base-1#mohid-base-1)). Optionally, you can extract only specified datasets.

## Input data file

The input data file must be placed in the working folder and must be named *Convert2netcdf.dat*.

```
HDF_FILE             : [char]    !Path to the HDF5 file to be converted
HDF_SIZE_GROUP       : [char]    !HDF group containing the HDF dataset from where
                                 !the dimensions will be read
HDF_SIZE_DATASET     : [char]    !Name of the HDF dataset from where the dimensions
                                 !will be read
HDF_TIME_VAR         : [char]    !Name of the HDF group containing the
                                 !time values
HDF_READ_LATLON      : int       !Type of projection (0 - metric, 1 - lat/lon)
HDF_READ_SIGMA       : int       !Type of vertical coordinate (0 - cartesian, 1 - sigma)

IMPOSE_MASK          : bool      !In case we want to specify the mask variable.
HDF_MASK             : [char]    !If we chose to specify a mask variable, then we define here
                                 !the variable name
HDF_MASK_IS_3D       : bool      !Is the mask variable 3D? yes or no.
RESULTS_ARE_2D       : bool      !Sometimes the mask variable is 3D, but the results

CONVERT_EVERYTHING   : int       !0 - Don't convert everything, 1 - Convert everything
VGROUP_TO_CONVERT    : [char]    !Path to dataset(s) to converto. (e.g. /Results/water level)

NETCDF_FILE          : [char]    !netcdf file to be created
NETCDF_TITLE         : [char]    !netcdf file title
NETCDF_CONVENTION    : [char]    !netcdf naming convention
NETCDF_VERSION       : [char]    !netcdf library version
NETCDF_HISTORY       : [char]    !netcdf file history
NETCDF_SOURCE        : [char]    !netcdf file source
NETCDF_INSTITUTION   : [char]    !netcdf institution
NETCDF_REFERENCES    : [char]    !netcdf references
NETCDF_DATE          : int       !current year
```

## Example

```
HDF_FILE             : WaterProperties_1.hdf5
HDF_SIZE_GROUP       : /Grid
HDF_SIZE_DATASET     : WaterPoints3D
HDF_TIME_VAR         : Time
HDF_READ_LATLON      : 1
CONVERT_EVERYTHING   : 0
VGROUP_TO_CONVERT    : /Results/salinity

NETCDF_FILE          : WaterProperties_1.nc
NETCDF_TITLE         : WaterProperties file converted to netcdf
NETCDF_CONVENTION    : CF-1.0
NETCDF_VERSION       : 3.6.1
NETCDF_HISTORY       : Converted from hdf5
NETCDF_SOURCE        : File created by MOHID-Convert2ncdf
NETCDF_INSTITUTION   : Instituto Superior Técnico (IST)
NETCDF_REFERENCES    : http://www.mohid.com/
NETCDF_DATE          : 2006
```

```
CONVERT_EVERYTHING            : 0
```

```
<begin_groups>
/Results/temperature
/Results/salinity
<end_groups>
```
