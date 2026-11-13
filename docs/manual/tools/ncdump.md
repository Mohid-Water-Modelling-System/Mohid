# Ncdump

**Ncdump** is a little [netcdf](../archive/linux/netcdf.md) utility that allows to dump in ASCII the data contained in the netcdf file.

## Example

ASCII dumps the datafile header only:

```
> ncdump -h ALADIN_BULKCLOUD_OPASYMP_19723_20088.nc
```

```
netcdf ALADIN_BULKCLOUD_OPASYMP_19723_20088 {
dimensions:
       time_counter = UNLIMITED ; // (2928 currently)
       lat = 245 ;
       lon = 251 ;
variables:
       float soclotot(time_counter, lat, lon) ;
               soclotot:long_name = "Couverture nuageuse" ;
               soclotot:unit = "%" ;
               soclotot:missing_value = 0.f ;
       double time_counter(time_counter) ;
               time_counter:units = "seconds since 1950-01-01 00:00:00" ;
               time_counter:calendar = "gregorian" ;
               time_counter:title = "Time" ;
               time_counter:long_name = "Time axis" ;
               time_counter:time_origin = "1950-01-01 00:00:00" ;
       float lat(lat) ;
               lat:units = "degrees_north" ;
               lat:valid_min = 0.f ;
               lat:valid_max = 0.f ;
               lat:long_name = "Latitude" ;
               lat:nav_model = "Default grid" ;
       float lon(lon) ;
               lon:units = "degrees_east" ;
               lon:valid_min = 0.f ;
               lon:valid_max = 0.f ;
               lon:long_name = "Longitude" ;
               lon:nav_model = "Default grid" ;
}
```
