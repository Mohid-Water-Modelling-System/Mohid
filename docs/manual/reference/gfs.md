# GFS

This article is about **GFS** data conversion, data can be accessed from [http://nomads.ncep.noaa.gov/](http://nomads.ncep.noaa.gov/).

The Global Forecasting System (GFS) files are packed under the GRIB protocol, following the World Meteorological Organization (WMO) protocol and is the primary format at the operational centers (ECMWF,NCEP). There are two different types of GRIB files: **GRIB1** and **GRIB2.** GRIB2 is the most recent format, because of that historical simulation files (e.g. reanalysis) are usually found in GRIB1 while the forecasts and recent analysis files are in GRIB2 format.

Reanalysis data, in the GRIB1 format, can be collected from the NCEP FNL Operational Model Global Tropospheric Analyses [http://dss.ucar.edu/datasets/ds083.2/](http://dss.ucar.edu/datasets/ds083.2/), from July 1999. To collect the data it is needed to become a registered user.

The GFS forecasts, in the GRIB2 format, can be downloaded from [ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/](ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/). This ftp has a different folder for each simulation. The name indicates the begining of the simulation (**gfs.yYYYYMMDDHH**). For example, the folder named gfs.2009081406 indicated results produced at 6 a.m. (UTC) of the 14th of August 2009. Inside each folder, there can be found different files, but our interest is mainly in the files named with the following format **gfs.tHHz.pgrb2bfHZ** that indicate that are gfs files, produced at the simulation hour HH and that correspond to the HZ forecast hour (after the beginning of the simulation). Results are produced every 3 hours and the size of each file is around 40 MB.

More information on GFS files can be found on [http://www.nco.ncep.noaa.gov/pmb/products/gfs/](http://www.nco.ncep.noaa.gov/pmb/products/gfs/) .

## Ungrib

The GRIB1 files can be converted into ASCII or binary files using the software **wgrib.exe** while GRIB2 files can be converted into ASCII, binary or NETCDF using the software **wgrib2.exe**. Once the files are converted into ASCII text files they would be converted into HDF5 files through a Matlab code.

The wgrib and wgrib2 software can be downloaded at the NOAA´s Climate Prediction Center website [http://www.cpc.noaa.gov/products/wesley/](http://www.cpc.noaa.gov/products/wesley/)

### wgrib.exe

Before using the wgrib.exe (in this case the Windows version) it would be needed the additional software grep, that can be downloaded from:

[http://gnuwin32.sourceforge.net/packages/grep.htm](http://gnuwin32.sourceforge.net/packages/grep.htm)

Once installed it is needed to include a path to the variables in order to be recognised the grep command. The path is set in Control panel -> System -> Advanced -> Environmental variables -> Path then adding a line at the end with the location of the software i.e. C:\\Program Files (x86)\\GnuWin32\\bin

To obtain the variables included in the Grib file type:

```
wgrib -s  grib_file
```

To write each of them into a separate ascii file:

```
wgrib -s  grib_file | grep "variable_name" | wgrib -s -i -text grib_file -o output_file.dat
```

To obtain the grid

```
wgrib -V -o /dev/null  -d 1 grib_file
```

### wgrib2.exe

The wgrib2.exe is a command-line program and presents a wide list of options. The best option is to convert the grib2 file into a netcdf file using [ConvertToHDF5 Generic\_NETCDF\_CF\_FILES](../tools/converttohdf5.md)

Also can be obtained an ascii file and in order to reduce the computing space, the extracted domain can be limited through the option -small\_grib. A typical command-line using this option would be:

```
wgrib2 20090101/gfs.t00z.pgrb2f00 -small_grib -10:11 30:45 Smallergrib00.grib
```

where -10:11 is the longitude interval, 30:45 is the latitude interval and Smallergrib00.grib correspond to the new grib file.

It has been created an script to extract the variables needed for modelling purposes, that can be found in \\\\einstein\\www\\ftp.mohid.com\\operational\\meteoIST\\GFSGRIB2

Once executed, the outputs are .dat files each containing a different property with the following correspondences:

- HR2m - Relative Humidity 2 m above ground 

$$
%
$$


- LAND - Land Cover (1=land, 0=sea) 

$$
Proportion
$$


- MSLP - Pressure Reduced to Mean Sea Level 

$$
Pa
$$


- TCDC - Total Cloud Cover 

$$
%
$$


- TMP2m - Temperature 2 m above ground 

$$
K
$$


- U - U-Component of Wind 10 m above ground 

$$
m/s
$$


- V - V-Component of Wind 10 m above ground 

$$
m/s
$$



## ASCII to HDF5 format

Once obtained, the ASCII results outcoming from the wgrib2.exe have to be converted into HDF5 to serve as input for the Mohid atmosphere module. This operation is performed through a Matlab code available in:

\\\\einstein\\www\\ftp.mohid.com\\operational\\meteoIST\\GFSGRIB2\\Matlab

The Matlab version for using this code have to be the R2008b version, or previous, due to the uncompatibility found with new versions (i.e. R2009a) with the HDF5 libraries version. The new Matlab versions work with the 1.8.1 HDF5 libraries version while the previous Matlab versions work with 1.6.5 HDF5 libraries. All the Mohid modelling system is compatible with 1.6.5 HDF5 libraries.
