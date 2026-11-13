# Bathymetry

## Overview

A Mohid Water bathymetry file is a [Grid Data](grid-data.md) file with information on the horizontal grid (geographical positioning, grid size, grid rotation, etc) and the bathymetry values to be used by the model. Bathymetry values are generally positive and intertidal areas normally have negative values (depending on the referential being used - normally is hydrographic zero) and land grid cells have a default value of -99.

## Generating a new bathymetry file

A bathymetry file can be created using the [Digital Terrain Creator](../../tools/digital-terrain-creator.md) which has a graphical user interface included in [MohidGIS](../../tools/gui/mohid-gis.md). The bathymetry file can be created via Mohid GIS, or from executing Digital Terrain Creator the command line. It can also be constructed manually or using another software as long as it fulfills the requirements established by the [Grid Data](grid-data.md) format.

### Required data

To create a new bathymetric grid data it is required to have:

- if necessary, a coastline of the study area (in MOHID [polygon](polygons.md) format)
- bathymetric data in x, y and z coordinates(in MOHID [xyz](xyz-points.md) format)
- a MOHID [grid](grid.md) file

### Data sources

Below are listed some places where one can find data available to construct a bathymetry file.

#### Coast lines

[NOAA - NGDC-WDC MGG, Boulder-Shoreline-Coastline Databases](http://rimmer.ngdc.noaa.gov/coast/) - In this website coastlines with different resolution can be obtained for any place in the World with a maximum resolution of 1:250.000. Note: this data has to be processed in order to fulfill the [Polygon](polygons.md) format.

#### Bathymetric data

[General Bathymetric Chart of the Oceans (GEBCO)](http://www.bodc.ac.uk/projects/international/gebco/) - In this website data can be downloaded, once registered, bathymetries of any place of the World with a 1´x1´ resolution. However, it has to be considered, that these data are not very reliable for coastal areas. In 2009, they released an updated version of their 30 minute GEBCO - grid which also includes satellite data. A description of both products can be found here [

$$
1
$$

](http://www.gebco.net/data_and_products/gridded_bathymetry_data/). Note: this data has to be processed in order to fulfill the [XYZ points](xyz-points.md) format.

[ETOPO](http://www.ngdc.noaa.gov/mgg/global/relief/) - In these website, a complete database with 2'x2' or 5'x5' regular grid resolution can be downloaded. Note: this data has to be processed in order to fulfill the [XYZ points](xyz-points.md) format.

### Converting data to MOHID format

The files with raw information must be converted to MOHID format so that they can be recognised by Mohid GIS and Digital Terrain Creator. MOHID has a simple support tool able to convert data from recognised sources (mentioned above) called [ConvertToXYZ](../../tools/converttoxyz.md).
