# ModuleNETCDF

## Overview

ModuleNETCDF is an interface module between the [netcdf](../../archive/linux/netcdf.md) libraries and [Mohid Framework](../../development/programming/README.md#mohid-framework#mohid-framework), and it was developed in order to read and write netcdf files in a managed and straighforward way.

## Using ModuleNETCDF to write a netcdf file

Here is the steps to be taken in order to write a new netcdf file.

### Ask for access code to create a new file

```
call GetNCDFFileAccess (NCDF_CREATE, NCDF_READ, NCDF_READWRITE)
```

### Create a new instance of ModuleNETCDF

```
call subroutine ConstructNETCDF (ObjNCDFID, FileName, Access, STAT)
```

### Write netcdf file header (CF-1.0)

```
call NETCDFWriteHeader(NCDFID, Title, Convention, Version, History, &
                       iDate, Source, Institution, References, STAT)
```

### Write temporal dimension

```
call NETCDFWriteTime(NCDFID, InitialDate, nInstants, Times, STAT)
```

### Write spatial dimension

```
call NETCDFSetDimensions (NCDFID, IUB, JUB, KUB, STAT)
```

```
call NETCDFWriteLatLon(NCDFID, Lat, Lon, Lat_Stag, Lon_Stag, GeoCoordinates, STAT)
```

### Write data

```
call NETCDFWriteData(NCDFID, Name, LongName, StandardName, Units, &
                     ValidMin, ValidMax, MinValue, MaxValue,      &
                     OutputNumber, ARRAY*, STAT)
```

```
ARRAY can be 2D or 3D and integer or real (single or double precision)
```
