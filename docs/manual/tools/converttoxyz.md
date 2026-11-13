# ConvertToXYZ

**ConvertToXYZ** is a small tool that extracts baseline three-dimensional data into [XYZ points files](../io/formats/xyz-points.md). It is useful, for instance, to extract baseline data to create bathymetries (or topography in land projects).

## Quick Start

1. Map a network drive to the data file (ETOPO, GEBCO, NOAA etc...).
2. Create and edit a *ConvertToXYZ.dat* file. Save it in the same path of the executable.
3. Run the executable.

## Options file

***WARNING:** Only use integer coordinate lat/lon numbers*.

```
FILE_TYPE: is the parameter that defines the type of input data.
It accepts values from 1 to 6 that indicates the following possible origins of data:
1 NASA SRTM (Shuttle Radar Topography Mission)
2 ETOPO 2 (every 2')
3 ETOPO 5 (every 5')
4 Shoreline (Mapgen)
5 GEBCO  (every 30" or 1')
6 ASCII  (ArcView.....)
7 NASA SRTM30 (every 30")
INPUT_FILENAME: input data file,
OUTPUT_FILENAME: where the results will be written,
TOP, BOTTOM, LEFT, RIGHT: Latitude and longitude minimum and maximum values (integers only!)
WRITE_AS_BATHYMETRY: indicates if equal to 1 that the output is written in bathymetry format
WRITE_AS_POLYGONS: indicates that would produce outcoming files with polygons format *.xy when 1.
```

### ETOPO2

Here's what the options *ConvertToXYZ.dat* file looks like:

```
<begin_window>
FILE_TYPE                   : 2
INPUT_FILENAME              : Z:\etopo2.dat
OUTPUT_FILENAME             : WestIberiaET2.xyz
TOP                         : 46
BOTTOM                      : 33
LEFT                        : -14
RIGHT                       : -5
MAXIMUM                     : 1
MINIMUM                     : -9999
WRITE_AS_BATHYMETRY         : 1
CONVERT_TO_HDF5             : 0
<end_window>
```

### NOAA Shoreline/Coastline

The file ConvertToXYZ.dat to convert a coastline downloaded from the [NOAA Marine Geology & Geophysics Shoreline/Coastline Data](http://www.ngdc.noaa.gov/mgg/shorelines/shorelines.html) would look like that:

```
<begin_window>
FILE_TYPE                   : 4
INPUT_FILENAME              : Aysen-lowres.dat
OUTPUT_FOLDER               : C:\Aplica\Ecomanage\Data\NGDC_NOAA_shorelines\Software\
OUTPUT_NAME                 : Aysen-lowres
WRITE_AS_POLYGONS           : 1
MAX_BUFFER_SIZE             : 250000
<end_window>
```

### GEBCO

For [GEBCO](http://www.gebco.net/) bathymetric downloaded data the file would look like:

```
<begin_window>
FILE_TYPE                    : 5
INPUT_FILENAME               : 60W30S_40W10S.cdf
OUTPUT_FILENAME              : 60W30S_40W10S.xyz
MAXIMUM                      : 1
MINIMUM                      : -9999
WRITE_AS_BATHYMETRY          : 1
CONVERT_TO_HDF5              : 0
<end_window>
```

### ASCII Grid

Use this option if you have and ASCII grid. ASCII grids can be produced from raster datasets for example in ArcGIS ("Raster to ASCII").

```
<begin_window>
FILE_TYPE                   : 6
INPUT_FILENAME              : DTM.asc
OUTPUT_FILENAME             : DTM.xyz
<end_window>
```

Optionally it is possible now to reduce resolution of the xyz. This is quite handy when dealing with very detailed DTM's that are not openable in GUI because of the file size. This can be done by defining the step to pick points (2 will select half of the points and so on).

```
<begin_window>
FILE_TYPE                   : 6
INPUT_FILENAME              : DTM.asc
OUTPUT_FILENAME             : DTM.xyz
REDUCE_RESOLUTION           : 1
REDUCE_RESOLUTION_STEP      : 2
<end_window>
```

Other KEYWORDS that can be useful:

IGNORE\_NEGATIVE\_VALUES : 1

ignores negative values and these are not written to output file.

ADD\_FACTOR  : xx

adds the xx value to the original z (e.g. units conversion). it can have negative values

MULTIPLY\_FACTOR  : xx

multiplies the xx value to the original z (e.g. units conversion). it can have negative values

### SRTM30

*Caution: make sure the ORIGIN\_X and ORIGIN\_Y correspond to the BOTTOM-LEFT origin of the srtm30 original source file! [Here's](http://www.dgadv.com/srtm30/) a map of srtm30 tiles.*

```
<begin_window>
FILE_TYPE                   : 7
TOP                         : 46
BOTTOM                      : 33
LEFT                        : -14
RIGHT                       : -5
MAXIMUM                     : 1
MINIMUM                     : -9999
INPUT_FILE                  : \\kepler\DataCenter1\DataCenter\DadosBase\Oceanografia\Batimetria\SRTM30\Data\w020n40.Bathmetry.srtm
XYZ_OUTPUT_FILENAME         : WestIberiaSRTM.xyz
ORIGIN_X                    : -20
ORIGIN_Y                    : -10
<end_window>
```

### NASA Topography

1) If you have NASA binary files (.hgt):

For land topography look on \\\\kepler\\DataCenter1\\DataCenter\\DadosBase\\Terra\\Topografia\\NASA

In Eurasia folder you will find data for Europe. Copy the .zip files you need (or all) to your Hard Disk and unzip them.

Run the ConvertToXYZ with a structure adapted to your needs:

```
<begin_window>
FILE_TYPE             : 1
LEFT                  : -9.28
RIGHT                 : -9.056
TOP                   : 38.99
BOTTOM                : 38.75
INPUT_FOLDER          : D:\Dados\NASA\
XYZ_OUTPUT_FILENAME   : D:\Aplica\Trancao\Trancao.xyz
SAMPLE_POINTS         : 1
CONVERT_TO_UTM        : 1
<end_window>
```

2) If you have NASA grids (e.g. SRTM from site):

Convert DEM grid to ascii (e.g. using ArcGIS "Raster to ASCII")

Convert ascii to .xyz using this tool (see example above).

## See also

- [Digital Terrain Creator](digital-terrain-creator.md)
