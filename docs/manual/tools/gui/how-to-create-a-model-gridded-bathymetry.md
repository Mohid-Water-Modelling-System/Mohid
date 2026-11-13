# How to create a model gridded bathymetry?

How to create a model gridded bathymetry?

## Step 1 - Extract Ocean bathymetry data from ETOPO 2'

### Google Earth

[Make sure you have Google Earth correctly configured!!](../../archive/web/google-earth.md)

- Obtain the ETOPO 2' extraction domain limits by reading the Max and Min values of Longitude and Latitude for your domain from Google Earth.



### Etopo data file

- The internet address where to download the data is [http://www.ngdc.noaa.gov/mgg/global/global.html](http://www.ngdc.noaa.gov/mgg/global/global.html);
- Download the ETOPO2V2 raw binary data to use with the [ConvertToXYZ](../converttoxyz.md) program OR, if possible, rather use the intranet address which is *\\\\kepler\\DataCenter1\\DataCenter\\DadosBase\\Oceanografia\\Batimetria\\Etopo\\Etopo2\\Data*



A sample configuration file **ConvertToXYZ.dat** is:

```
<begin_window>
FILE_TYPE                   : 2
INPUT_FILENAME              : \\kepler\DataCenter1\DataCenter\DadosBase\Oceanografia\Batimetria\Etopo\Etopo2\Data\etopo2.dat
OUTPUT_FILENAME             : SLouisCoast_ETOPO2.xyz
TOP                         : 17
BOTTOM                      : 14
LEFT                        : -18
RIGHT                       : -14
MAXIMUM                     : 1
MINIMUM                     : -9999
WRITE_AS_BATHYMETRY         : 1
CONVERT_TO_HDF5             : 0
<end_window>
```





- Sample of XYZ data extracted from ETOPO 2'. How the bathymetry XYZ file should look like.



### MOHID GIS

- Open new project in Mohid GIS



- Give name and location to new project, then save



- Add new data item



- Select data item and open



- Right-click data, zoom select, right-click, properties



- Click on settings to set the colar palette and color range



- Choose min depth (negative if above sea level), max depth and color palette



- IMPORTANT: Click *Save All* from the *Data Items* menu. Otherwise you risk to loose all your work!



## Step 2 - Create a coast line

### Create coast line in Google Earth

- Create a polygon shaped like the line of the coast



- Right-click on polygon, then "save as"



- Save the file in Kml format (NOT in Kmz format)



### Run Kml2mohid

- Get the [Kml2Mohid](../../development/programming/kml2mohid.md) program, edit a batch file, and run the program



- Successfully run [Kml2Mohid](../../development/programming/kml2mohid.md) program.



- Sample polygon (.xy) file



### Import coastline polygon in Mohid

- Click *Data Items* menu, then select *Add Data Items...*



- Select folder, select polygon file type, select file and click open



- Right-click polygon data, then select *Properties*



- Click on settings



- Select *Constant* and choose a color



- Very important: *Save All* from the *Data Items* menu



- General look of ETOPO points data layer and coastline polygon layer



## Step 3 - Create the model domain regular grid in MOHID GIS

### Define the grid

- Choose bottom-left origin for the model domain



- Choose model domain grid resolution with [Google Earth](../../archive/web/google-earth.md)



- Click on *Data Item* menu, then select *New Data Item*



- Select *Simple Grid*, then click on *Browse*



- Give name, select *Grid file* as file type, click on *Save*



- Click *Ok*



- Right-click the Grid layer, then select *Properties* from the menu



- Define origin, define resolution, define grid size, choose *Simple Geographic* coordinates type (and NOT *Geographic* as seen in the figure. CAUTION!)



- Very important: *Save All* from the *Data Items* menu



- General look of a bathymetry points data from ETOPO, a coastline and land polygon, and a model domain grid



### Create Digital Terrain from Mohid GIS

- Click on menu *Tools*, select *Create Digital Terrain ...*



- Select grid layer, data points layers, land polygon layers. Then click *Browse*



- Choose *Options* menu, select *Triangulation*, select *All points*, click on *Run*



- Successful *Done* message



## Admire the result

You can now proceed to create a model simulation after you [create the tidal gauges](how-to-generate-tide-for-mohid.md)



## Next steps

- [How to generate tide for Mohid?](how-to-generate-tide-for-mohid.md)

## Previous steps

- [Software to run Mohid](../../troubleshooting/software-to-run-mohid.md)

## About the screen capturing tool

The screen capturing tool used for this tutorial was the [Jing Project](http://www.jingproject.com/).
