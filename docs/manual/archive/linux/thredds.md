# Thredds

**Thredds** (Thematic Realtime Environmental Distributed Data Services) is a tomcat servlet that serves xml catalogs of [openDAP](opendap.md) datasets.

## Setup

A working [tomcat](tomcat.md) server is a prerequisite. Simply **copy the thredds.war file** obtained from the thredds homepage **into the** [tomcat](tomcat.md) **webapps directory**.

### data.mohid.com

- **How do I restart the [tomcat](tomcat.md) and thredds server in data.mohid.com?**

```
> cd /home/guillaume/Software/las7/tomcat/tomcat/bin
> ./Tomcat5.sh stop
> ./Tomcat5.sh start
```

- **How do I know if tomcat is really running?**:

```
> netstat -tlnap | grep 8080
```

If it's running then a port 8080 must appear in the results of the command.

### VmMrtSrv06

- **192.168.20.136** when it is resident on the LAN.
- To configure the tomcat to start during boot:

```
> sudo vim /etc/rc.d/rc.local
rc.local> /home/user/Download/apache-tomcat-6.0.18/bin/startup.sh
```

#### Catalogs architecture

We propose the following architecture as a good practice for maintaining catalogs.

## Remote web-based administration

Use the following link

- [Remote administration](http://opendap.mohid.com:8080/thredds/admin/debug)

Username and password are the usual maretec administrator login and password.

## Directory tree structure

### $TOMCAT\_HOME/content/thredds

- **threddsConfig.xml**: configuration xml file of the Thredds server.
- **catalog.xml**: catalog file containing links to all the subcatalogs containing the netcdf data.
    -   **pcomsCatalog.xml**: catalog file pointed by the main catalog.

#### threddsConfig.xml

Simply edit the file and enable the netcdf subset service.

```
#threddsConfig.xml
 <NetcdfSubsetService>
   <allow>true</allow>
   <dir>/data/tmp/thredds/ncSubsetCache/</dir>
   <scour>15 min</scour>
   <maxAge>30 min</maxAge>
 </NetcdfSubsetService>
```

#### catalog.xml

```
#catalog.xml
<?xml version="1.0" encoding="UTF-8"?>
<catalog name="THREDDS Server Default Catalog : You must change this to fit your server!"
       xmlns="http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0"
       xmlns:xlink="http://www.w3.org/1999/xlink">
  <catalogRef xlink:title="PCOMS Catalog" xlink:href="/opendap/pcoms/pcomsCatalog.xml" name=""/>
</catalog>
```

#### pcomsCatalog.xml

```
#pcomsCatalog.xml
<?xml version="1.0" encoding="UTF-8"?>
<catalog xmlns="http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0" xmlns:xlink="http://www.w3.org/1999/xlink"
  name="Unidata THREDDS-IDD NetCDF-OpenDAP Server" version="1.0.1">

 <service name="latest" serviceType="Resolver" base="" />
 <service name="all" serviceType="Compound" base="">
   <service name="ncdods" serviceType="OPENDAP" base="/thredds/dodsC/" />
   <service name="HTTPServer" serviceType="HTTPServer" base="/thredds/fileServer/" />
   <service name="ncss" serviceType="NetcdfSubset" base="/thredds/ncss/grid/" />
 </service>

 <dataset name="PCOMS Model Data">
   <metadata inherited="true">
     <serviceName>all</serviceName>
     <authority>maretec.ist.utl.pt</authority>
     <dataType>Grid</dataType>
     <dataFormat>NetCDF</dataFormat>
     <documentation type="rights">Freely available</documentation>
     <documentation xlink:href="http://www.mohid.com/operational/" xlink:title="Portuguese Coast Operational Modeling System"></documentation>
     <creator>
       <name vocabulary="DIF">Maretec/IST</name>
       <contact url="http://www.mohid.com/" email="http://www.mohid.com/contacts.htm" />
     </creator>
     <publisher>
       <name vocabulary="DIF">Maretec/IST</name>
       <contact url="http://www.mohid.com/" email="http://www.mohid.com/contacts.htm" />
     </publisher>
     <timeCoverage>
       <start></start>
       <end>present</end>
     </timeCoverage>
   </metadata>

   <datasetScan name="PCOMS Data" ID="PcomsDatasetScan"
                path="pcoms" location="/opendap/pcoms/"
                harvest="true">
     <metadata inherited="true">
       <documentation type="summary">The PCOMS model consists of a downscaling of the Mercator-Ocean PSY3V1 north-atlantic solution for the Portuguese coast, where the fes2004 tide solution is added, and the surface forcing is performed using the atmospheric forecasts from the 9 km resolution MM5 model run at http://meteo.ist.utl.pt. The model is run in full baroclinic mode with an horizontal resolution of 6 km and with 43 layers in the vertical axis with finer resolution of 5 m near the surface.</documentation>
       <geospatialCoverage zpositive="down">
         <northsouth>
           <start>34.4</start>
           <size>10.6</size>
           <units>degrees_north</units>
         </northsouth>
         <eastwest>
           <start>-12.6</start>
           <size>7.1</size>
           <units>degrees_east</units>
         </eastwest>
         <updown>
           <start>0</start>
           <size>5500</size>
           <units>m</units>
         </updown>
       </geospatialCoverage>
       <variables vocabulary="CF-1.0" />
       <variables vocabulary="">
         <variable name="temperature" vocabulary_name="temperature" units="degC">temperature</variable>
         <variable name="salinity" vocabulary_name="salinity" units="1e-3">salinity</variable>
         <variable name="ssh" vocabulary_name="ssh" units="m">sea water level</variable>
         <variable name="u" vocabulary_name="u" units="m/s">east-west current velocity</variable>
         <variable name="v" vocabulary_name="v" units="m/s">north-south current velocity</variable>
         <variable name="vm" vocabulary_name="vm" units="m/s">current velocity modulus</variable>
       </variables>
     </metadata>

     <filter>
       <include wildcard="PCOMSV1_*.nc" />
     </filter>
     <addID/>
     <sort>
       <lexigraphicByName increasing="false"/>
     </sort>
     <addLatest/>
     <addDatasetSize/>
     <addTimeCoverage datasetNameMatchPattern="PCOMSV1_([0-9]{4})([0-9]{2})([0-9]{2}).*.nc$"
                      startTimeSubstitutionPattern="$1-$2-$3T00:00:00"
                      duration="3 days" />
   </datasetScan>
 </dataset>
</catalog>
```

## Building thredds urls

### Dataset

- Viewing a ${catalog} file tree named ${catalog}.xml:

```
http://opendap.mohid.com:8080/thredds/${catalog}.html
```

- Viewing a specific ${dataset} thredds metadata:

```
http://opendap.mohid.com:8080/thredds/${catalog}.html?dataset=${dataset.id}
```

- Viewing a specific ${dataset}:

```
http://opendap.mohid.com:8080/${dataset.service.base}/${dataset.urlpath}
```

### DatasetScan

- Viewing a particular ${datasetscan} tree from any ${catalog} in the thredds server:

```
http://opendap.mohid.com:8080/thredds/catalog/${datasetscan.path}/catalog.html
```

- Viewing a specific ${datasetscan.filteredfilename} thredds metadata:

```
http://opendap.mohid.com:8080/thredds/catalog/${datasetscan.path}/catalog.html?dataset=${datasetscan.id}/${datasetscan.filteredfilename}.nc.html
```

- Viewing a specific ${datasetscan.filteredfilename} dataset:

```
http://opendap.mohid.com:8080/${datasetscan.service.base}/${datasetscan.path}/${datasetscan.filteredfilename}.nc.html
```

## A sample catalog file

`[XML]`

```
<?xml version="1.0" encoding="UTF-8"?>

<catalog name="MOHID THREDDS Server Catalog"
       xmlns="http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0"
       xmlns:xlink="http://www.w3.org/1999/xlink">

 <service name="MOHID DODS" serviceType="OpenDAP" base="/thredds/dodsC/" />
 <datasetRoot path="test" location="content/testdata/"/>

 <dataset name="Test Single Dataset" ID="testDataset" serviceName="MOHID DODS"
          urlPath="test/testData.nc"/>

 <datasetScan name="Test all files in a directory" ID="testDatasetScan"
              path="testAll" location="content/testdata">

   <metadata inherited="true">
     <serviceName>MOHID DODS</serviceName>
   </metadata>

   <filter>
     <include wildcard="*.nc"/>
   </filter>

 </datasetScan>

 <catalogRef xlink:title="Test Enhanced Catalog" xlink:href="enhancedCatalog.xml" name=""/>

</catalog>
```

## Troubleshoot

### Tomcat works, Thredds doesn't

**Q: the Thredds won't work in [opendap.mohid.com](http://opendap.mohid.com:8080/thredds), but the [tomcat](http://opendap.mohid.com:8080) seems to be working. What should I do?**

A: There's a possibility that the old LAS tomcat is running (booted at startup) and doesn't lets the new thredds to load properly. Checkout the [LAS](las.md) wiki entry and see how it goes to deactivate the old LAS. Then stop and starp the Thredds as shown in the sections above.

### Thredds cache

**Q: The Thredds aggregated timeserie won't correct the displayed values of the time variable. What can I do?**

A: You need to clean your Thredds cache. You can check-out them here:

```
/home/user/Download/apache-tomcat-6.0.18/content/thredds/aggcache
/home/user/Download/apache-tomcat-6.0.18/content/thredds/cacheAged
/home/user/Download/apache-tomcat-6.0.18/content/thredds/cache
```

Remove any suspicious file related with your aggregated timeserie.

### Controling the time variable in Thredds NcML aggregated datasets

**Q: The thredds aggregated timeserie shows all the values indexed with time. But I only want a subset of the time-indexed values to be aggregated in order to avoid overlaps. How can I do that?**

A: Add the attribute "ncoords" in the xml tag "netcdf". As an example, you could change the following line

```
<netcdf location="/opendap/pcoms/20090323_PCOMSV1.nc"/>
```

with this line

```
<netcdf ncoords="8" location="/opendap/pcoms/20090323_PCOMSV1.nc"/>
```

That way, only the first 8 time-indexed values of every variables will be concatenated.

## See also

- [Tomcat](tomcat.md)
- [OpenDAP](opendap.md)
- [LAS](las.md)
- [WCS](../../reference/wcs.md)

## External references

- [Getting started](http://www.unidata.ucar.edu/projects/THREDDS/tech/tutorial/GettingStarted.html)
- [MOHID Thredds data server](http://opendap.mohid.com:8080/thredds/catalog.html)
- [Thredds WCS server](http://www.unidata.ucar.edu/projects/THREDDS/tech/reference/WCS.html)
