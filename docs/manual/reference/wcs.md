# WCS

**WCS**, Web Coverage Service, is a standard that allows interoperability of georeferenced gridded data across web services by defining an xml schema of how the served datasets metadata should be presented.

The [Thredds](../archive/linux/thredds.md) server, besides gathering an opendap server and a catalog, it also allows an experimental wcs compliant metadata server.

## The purpose of yet another standard

1. [Netcdf](../archive/linux/netcdf.md) is a binary file format that allows to store gridded data.
2. The CF convention standard defines unambiguously how variable names and attributes should be set and organized within a netcdf dataset that is CF-compliant. Ex: The time dimension name should be set to *time* (and not *t*, nor *n*) and the sea water velocity modulus must have a *long\_name* attribute set to *sea water speed* (and not *velocity modulus*, nor *speed modulus*).
3. [OpenDAP](../archive/linux/opendap.md) is a netcdf files server, that reads and displays the netcdf metadata and the size of dimensions and variables contained in each netcdf file. Thus the client or the user, provided it already knows the name and address of the netcdf dataset, can **see** what the netcdf file contains **before** dowloading it. Furthermore, it allows to download only a subset of the netcdf gridded data in ASCII or netcdf format. Ex: it allows to download only the Atlantic ocean basin from a worldwide netcdf gridded dataset.
4. [Thredds](../archive/linux/thredds.md) is a netcdf files catalog server. It means that it allows the client or the user to **browse** through an organized list of netcdf files avaible to download via opendap or http. The hierarchical organization of the datasets is fully parametrized and controlled by the service provider.
5. **WCS** is an xml standard that defines unambiguously how georeferenced variables inside large datasets should be defined. It's similar to [kml](../archive/linux/kml.md), but more extensive. The interest in adopting such an xml standard, is to allow machines to search for data available within a geospatial range (mintime, maxtime, bottom-left corner, top-right corner). The [Thredds](../archive/linux/thredds.md) server has an experimental **WCS** implementation, that means that it maps it's own xml catalogs containing geospatial coverage metadata into the **WCS** standard.

## External references

- [Thredds WCS implementation](http://www.unidata.ucar.edu/projects/THREDDS/tech/reference/WCS.html)
- [WCS OGC official page](http://www.opengeospatial.org/standards/wcs)
