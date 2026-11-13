# OpenDAP

## Available software

### libdap

Developer libraries and tools. Provides the command-line **getdap** tool that allows to download opendap files.

### libnc-dap

Provides substitute netcdf libraries, yet opendap*\-able*.

### Hyrax

Is the opendap server. Requires handlers for specific file format such as netcdf, hdf4, hdf5, etc...

### Matlab loaddap

Is a matlab client which allows to access and read netcdf files from opendap urls. Requires the libdap libraries.

## Definitions

Below lie the different concepts associated with **OpenDAP**.

### [DAP](README.md#dap#dap)

is an [acronym](http://en.wikipedia.org/wiki/Acronym) for **D**ata **A**ccess **P**rotocol. It's a [protocol](http://en.wikipedia.org/wiki/Protocol_%28computing%29) for accessing and extracting information from large dataset files through the internet.

### [DAP](README.md#dap#dap) [server](http://en.wikipedia.org/wiki/Server_%28computing%29)

is a suite of [CGI](http://en.wikipedia.org/wiki/Common_Gateway_Interface) [scripts](http://en.wikipedia.org/wiki/Scripting_language) and [executables](http://en.wikipedia.org/wiki/Executable) which serve large dataset files. An http [daemon](http://en.wikipedia.org/wiki/Daemon_%28computer_software%29) [a.k.a](http://en.wikipedia.org/wiki/Aka) [httpd](http://en.wikipedia.org/wiki/Httpd), such as [Apache](http://www.apache.org/), working on a given port (typically port 80) is required.

### [DAP](README.md#dap#dap) handler

is a program that deals with a specific file type (ex: [netcdf](http://www.unidata.ucar.edu/software/netcdf/), [HDF](http://hdf.ncsa.uiuc.edu/), [JGOFS](http://www1.whoi.edu/), FreeForm etc ...). It requires a [DAP](README.md#dap#dap) [server](http://en.wikipedia.org/wiki/Server_%28computing%29). An [HDF5](http://hdf.ncsa.uiuc.edu/) handler is on the way (as of April 2006).

### [DAP](README.md#dap#dap) [client](http://en.wikipedia.org/wiki/Client_%28computing%29)

is a program that can communicate with a [DAP](README.md#dap#dap) [server](http://en.wikipedia.org/wiki/Server_%28computing%29). Reknowned clients are [ncBrowse](http://www.epic.noaa.gov/java/ncBrowse/) , [Ferret](http://ferret.wrc.noaa.gov/Ferret/), [matlab](http://www.opendap.org/download/matlab_gui.html). There's even a command line [client](http://en.wikipedia.org/wiki/Client_%28computing%29) *geturl* that comes with the [DODS](README.md#dods#dods) package.

### [DAP](README.md#dap#dap) interface

is an solution for users to access, browse and retrieve data from a data repository. For example, [LAS](las.md) is a DAP interface.

### **OpenDAP**

is an open source project for a network data access protocol. **OpenDAP** is also known as [DODS](README.md#dods#dods), **D**istributed **O**cean **D**atasets **S**ystem. It provides a package encapsulating a [DAP](README.md#dap#dap) server with file handlers such as [netcdf](http://www.unidata.ucar.edu/software/netcdf/) or [HDF](http://hdf.ncsa.uiuc.edu/). The official website is found at [http://www.opendap.org](http://www.opendap.org).

### [DODS](README.md#dods#dods) [URL](http://en.wikipedia.org/wiki/URL)

A [DODS](README.md#dods#dods) [URL](http://en.wikipedia.org/wiki/URL) is the http address that allows to browse and retrieve data from a [DAP](README.md#dap#dap) [server](http://en.wikipedia.org/wiki/Server_%28computing%29). A [DODS](README.md#dods#dods) [URL](http://en.wikipedia.org/wiki/URL) can append a query string, indicating the variable(s), the time and the space window. It can also query for the metadata of the file.

## End user's guide

This is a local user's guide. These are the official [user's guide](http://www.opendap.org/user/guide-html/guide.html) and [quick start guide](http://www.opendap.org/user/quick-html/quick.html).

### Quick start

#### From [ncBrowse](http://www.epic.noaa.gov/java/ncBrowse/)

1. Download and install the latest version of the [Java platform](http://java.sun.com/j2se/) from [Sun](http://java.sun.com/).
2. Download and install the latest version of [ncBrowse](http://www.epic.noaa.gov/java/ncBrowse/).
3. Run [ncBrowse](http://www.epic.noaa.gov/java/ncBrowse/).
4. Go to *file --> OpenDAP*, write down the [URL](http://en.wikipedia.org/wiki/Url) to the desired file [http://data.mohid.com/opendap/nph-dods/mercator-ist/Portugal/20061213\_Portugal\_WaterProperties.nc](http://data.mohid.com/opendap/nph-dods/mercator-ist/Portugal/20061213_Portugal_WaterProperties.nc), then hit *Accept*.
5. You can now browse at will the content of the remote file.

#### From a web-browser

1. From a web-browser, compose the address where the [DAP](README.md#dap#dap) files are being served. ex: [http://data.mohid.com/opendap/nph-dods/LASoutput.nc.html](http://data.mohid.com/opendap/nph-dods/LASoutput.nc.html) or [http://www.cdc.noaa.gov/cgi-bin/nph-nc/Datasets/](http://www.cdc.noaa.gov/cgi-bin/nph-nc/Datasets/). You can browse the directories to see the available files. [netcdf](http://www.unidata.ucar.edu/software/netcdf/) files are most likely to be found.
2. Select a [netcdf](http://www.unidata.ucar.edu/software/netcdf/) file by clicking on it or simply append it's full address in the address bar followed by *.html* suffix. ex: [http://data.mohid.com/opendap/nph-dods/mercator-ist/Portugal/20061213\_Portugal\_WaterProperties.nc.html](http://data.mohid.com/opendap/nph-dods/mercator-ist/Portugal/20061213_Portugal_WaterProperties.nc.html) or [http://www.cdc.noaa.gov/cgi-bin/nph-nc/Datasets/coads/1degree/global/std/air.mean.nc.html](http://www.cdc.noaa.gov/cgi-bin/nph-nc/Datasets/coads/1degree/global/std/air.mean.nc.html).
3. Select the desired variables by clicking on their select boxes. Then write their subdomain limits. For example, putting in latitude *

$$
12:2:24
$$

* selects the subdomain ranging from 12ºN to 24ºN with a 2º step. Only integers are allowed.
4. Get the data by clicking in *Get Ascii* or *Binary Data Object*. Respectively you will obtain the data in an [ASCII](http://en.wikipedia.org/wiki/Ascii) file or in [netcdf](http://www.unidata.ucar.edu/software/netcdf/) format (rewrite the extension *.nc.dods* to *.nc*).

### Quick ref

Here's a list of useful [suffixes](http://en.wikipedia.org/wiki/Suffix) for extracting metadata from netcdf files:

- *.help* generates a web based help document.
- *.html* generates a web page containing information and action buttons over the netcdf file.
- *.info* returns the netcdf file metadata.
- *.dds* returns information about the variables.
- *.das* returns information about the variable attributes.
- *.ver* returns the OpenDAP server version.
- *.asc* returns the netcdf file in [ASCII](http://en.wikipedia.org/wiki/Ascii) format.

Here's a quick reference of [suffixes](http://en.wikipedia.org/wiki/Suffix) for extracting data from netcdf files:

- *.asc?time* returns the variable time in ascii format.
- *.asc?time

$$
0:2:6
$$

* returns the variable time in ascii format but within a subset of stride 2.
- *.dods?time

$$
0:2:6
$$

* returns the same variable as above but in [binary](http://en.wikipedia.org/wiki/Binary_file) format (i.e. netcdf).
- *.dods?time

$$
2:3:8
$$

,sst

$$
2:3:8
$$



$$
43:47
$$



$$
215:234
$$

,ssh

$$
2:3:8
$$



$$
43:47
$$



$$
215:234
$$

* returns the time, the SST and the SSH variables for a given domain.
- *.dods?time

$$
2:3:8
$$

,sst

$$
2:3:8
$$



$$
43:47
$$



$$
215:234
$$

&pressure>50.0&pressure<100.0* returns the time and the SST variables for a given domain with the additional constraint that the variable pressure must range between 50.0 and 100.0 units.

## Netcdf libraries

The Opendap team maintains alternative netcdf libraries to be linked against any netcdf tool. They provide the ability to the netcdf libraries to read opendap urls besides local netcdf files.

## External references

- [OpenDap home](http://www.opendap.org/)
- [NCO homepage](http://nco.sourceforge.net/)
