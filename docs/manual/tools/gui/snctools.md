# Snctools

**Snctools** are [Matlab](../../archive/linux/matlab.md) scripts that allow to manipulate local and remote([OpenDAP](../../archive/linux/opendap.md)) [netcdf](../../archive/linux/netcdf.md) files. It requires [mexnc](../../archive/linux/matlab.md) to be installed prior. Proceed [here](../../archive/linux/matlab.md) for installation instructions.

## Commands

### nc\_dump

```
>>> help nc_dump
 NC_DUMP:  a Matlab counterpart to the NetCDF utility 'ncdump'.
     NC_DUMP(NCFILE) prints metadata about the netCDF file NCFILE.
     NC_DUMP(NCFILE,VARNAME) prints metadata about just the one netCDF variable
     named VARNAME.
```

### nc\_varget

```
>>> help nc_varget
 NC_VARGET:  Retrieve data from a netCDF variable.

 DATA = NC_VARGET(NCFILE,VARNAME) retrieves all the data from the
 variable VARNAME in the netCDF file NCFILE.

 DATA = NC_VARGET(NCFILE,VARNAME,START,COUNT) retrieves the contiguous
 portion of the variable specified by the index vectors START and
 COUNT.  Remember that SNCTOOLS indexing is zero-based, not
 one-based.  Specifying a -1 in COUNT means to retrieve everything
 along that dimension from the START coordinate.

 DATA = NC_VARGET(NCFILE,VARNAME,START,COUNT,STRIDE) retrieves
 a non-contiguous portion of the dataset.  The amount of
 skipping along each dimension is given through the STRIDE vector.

 NCFILE can also be an OPeNDAP URL if the proper SNCTOOLS backend is
 installed.  See the README for details.

 NC_VARGET tries to be intelligent about retrieving the data.
 Since most general matlab operations are done in double precision,
 retrieved numeric data will be cast to double precision, while
 character data remains just character data.

 Singleton dimensions are removed from the output data.

 A '_FillValue' attribute is honored by flagging those datums as NaN.
 A 'missing_value' attribute is honored by flagging those datums as
 NaN.  The exception to this is for NC_CHAR variables, as mixing
 character data and NaN doesn't really seem to work in matlab.

 If the named NetCDF variable has valid scale_factor and add_offset
 attributes, then the data is scaled accordingly.

 EXAMPLE:
 #1.  In this case, the variable in question has rank 2, and has size
      500x700.  We want to retrieve starting at row 300, column 250.
      We want 100 contiguous rows, 200 contiguous columns.

      vardata = nc_varget ( file, variable_name, [300 250], [100 200] );
```

## Examples

### dump

```
>> nc_dump('http://data.mohid.com/opendap/nph-dods/mercator-ist/Portugal/20080409_Portugal_WaterProperties.nc')
netcdf  http://data.mohid.com/opendap/nph-dods/mercator-ist/Portugal/20080409_Portugal_WaterProperties.nc {

dimensions:
	time = UNLIMITED ; (8 currently)
	lat = 177 ;
	lon = 117 ;
	lat_staggered = 178 ;
	lon_staggered = 118 ;
	depth = 42 ;

variables:
	double time(time), shape = [8]
		time:long_name = "time"
		time:standard_name = "time"
		time:units = "seconds since 2008-04-04 12:00:00"
	float lat(lat), shape = [177]
...
```

### read dimension variable

```
>> lon=nc_varget('http://data.mohid.com/opendap/nph-dods/mercator-ist/Portugal/20080409_Portugal_WaterProperties.nc', 'lon')

lon =

  -12.5700
  -12.5100
  -12.4500
  -12.3900
  -12.3300
  -12.2700
...
```

### read bathymetry variable

```
>> bathym=nc_varget('http://data.mohid.com/opendap/nph-dods/mercator-ist/Portugal/20080409_Portugal_WaterProperties.nc', 'bathymetry', [0 0], [116 176], [1 1]);
```

## External links

- [mexnc and snctools homepage](http://mexcdf.sourceforge.net/)
- [Snctools manual](http://mexcdf.sourceforge.net/tutorial/ch02.html)
