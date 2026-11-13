# How to create discharges in MOHID

In this tutorial we will continue to work on the simulation we built on the [previous tutorial](how-to-configure-a-2d-model-forced-with-tide-with-mohid.md) and we will show how to add an instantaneous discharge and a gradual discharge in the model, by setting up and configuring the [Module Discharges](../../water/modules/discharges.md).

## Step 1 - Identify the discharge location

Right now we're going to identify exactly what are the discharges locations. We want to simulate an accidental instantaneous discharge of a passive tracer and a continuous discharge as a function of time nearby a dam.

### In Google Earth

[Make sure you have Google Earth correctly configured!!](../../archive/web/google-earth.md)

- Select the placemark icon, give the discharge a name, and give it a position

Repeat the procedure, for the point nearby the dam.

- Save the location by right-clicking on them, and selecting *save as...*

- Select a folder, choose a filename, select *kml* filetype, click on *save* button

### In Mohid GIS

#### Convert from kml to xyz

- Get the [kml2mohid.exe](http://www.mohid.com/hydrogroup/data_sources/kml2mohid.exe) program and save it next to the kml file, in the same folder, then type from the command line shell:

```
> kml2mohid.exe discharges.kml
```

#### Add data as a points GIS layer, identify corresponding cell grids

- Open the Mohid GIS, and open the GIS project we left before:

- Add new data item from the *Data Items menu*

- Drag'n drop the xyz points data layer to the top, then right-click it and select *Properties*

- Set size to 6, then select *Pick Color*

- Pick a bright color, then click on *Ok*, twice

- Check the grid layer, and make sure that it is on top of the grid data layer, so it is displayed in the view area of the MOHID GIS

- Select the *area information* icon, then select a small area containing the discharge point. Then read the I and J grid coordinates for later usage. Repeat this process for every discharge point.

## Step 2 - Create a time-serie simulating the discharge

In order to create a time-serie that simulates a continuous discharge (to simulate, say, a river), you must open a new text file and format the data accordingly with this sample

```
TIME_UNITS  :  DAYS
SERIE_INITIAL_DATA  :   2008    10   27   0   0   0

Time  Flow     T     S
 days  m3/s     C   psu
    1     2     3     4
<BeginTimeSerie>
0  1520.00  11.4  0.0
1  1440.00  11.3  0.0
2  1340.00  11.2  0.0
3  1290.00  11.2  0.0
4  1270.00  11.1  0.0
5  1200.00  11.1  0.0
6  1160.00  11.0  0.0
7  1190.00  10.9  0.0
8  1320.00  10.9  0.0
<EndTimeSerie>
```

Then, once you have created a time-serie with the flow discharge, you simply save it in your [MohidWater](../../reference/mohidwater.md) project *GeneralData* folder.

Notes: the time containing column should be composed of integers. Thus, adjust the TIME\_UNITS parameter accordingly. The flow column is mandatory, whereas the T (temperature) and S (salinity) columns in the example are optional.

## Step 3 - Configure the Discharges module and run MohidWater

### Configure the [Module Discharges](../../water/modules/discharges.md)

- Open the project from the [last tutorial](how-to-configure-a-2d-model-forced-with-tide-with-mohid.md)

- Right-click the *run\_1* and select *Properties*

- Select in the *associated modules* the *module discharges*

- Edit the [Module Discharges](../../water/modules/discharges.md) according to the sample below, then save it. Notice we're adding two discharges: a continuous discharge, and an accidental discharge. Notice also that we're adding a single property, salinity. Properties are optional; there could be several as well as there could be none.

```
<begindischarge>
NAME                    : Dam discharge
DESCRIPTION             : Daily timeserie for the year 2008
I_CELL                  : 280
J_CELL                  : 222
DISCHARGE_UNIFORM       : 1
DATA_BASE_FILE          : D:\Aplica\Senegal-Test\GeneralData\Discharges\river_discharge_sample2.srw
FLOW_COLUMN             : 2
TIME_SERIE_COLUMN       : 1

<<beginproperty>>
NAME                   : salinity
UNITS                  : psu
DESCRIPTION            : freshwater accidental discharge
CONSTANT_CONC          : 1
DEFAULTVALUE           : 0.01
<<endproperty>>

<enddischarge>

<begindischarge>
NAME                   : Accidental discharge
DESCRIPTION            : Daily timeserie for the year 2008
I_CELL                 : 131
J_CELL                 : 165
DISCHARGE_UNIFORM      : 1
DATA_BASE_FILE         : D:\Aplica\Senegal-Test\GeneralData\Discharges\river_discharge_sample2.srw
FLOW_COLUMN            : 2
TIME_SERIE_COLUMN      : 1

<<beginproperty>>
NAME                   : salinity
UNITS                  : psu
DESCRIPTION            : freshwater accidental discharge
CONSTANT_CONC          : 1
DEFAULTVALUE           : 0.01
<<endproperty>>

<enddischarge>
```

- Edit the [Module Hydrodynamic](../../water/modules/hydrodynamic.md), add the following line, then save it

```
WATER_DISCHARGES        : 1
```

- If there's a property in the discharge (in occurence, there's salinity), then edit the [Module WaterProperties](../../water/modules/waterproperties.md), and add the following lines, then save it. The single most relevant pair parameter/value for discharges in the module WaterProperties is **DISCHARGES : 1**.

*Notes: Relatively to the project from the [last tutorial](how-to-configure-a-2d-model-forced-with-tide-with-mohid.md), we're adding the salinity property as a passive tracer. It will be advected with the same advection-diffusion method used for momentum and it will also yield a time serie for the same locations as the hydrodynamic's and it will also results every three hours in another HDF5 file. The initial model salinity will be homogeneous, with its default value.*

```
OUTPUT_TIME             : 0 10800.

ADV_METHOD_H            : 4
ADV_METHOD_V            : 4
TVD_LIMIT_H             : 4
TVD_LIMIT_V             : 4

TIME_SERIE              : 1
TIME_SERIE_LOCATION     : D:\Users\Soussou\GeneralData\TimeSeries\StLouisTimeSeries.dat

<beginproperty>
NAME                        : salinity
UNITS                       : ºC
DESCRIPTION                 : salinity by default is 36
DEFAULTVALUE                : 36

ADVECTION_DIFFUSION         : 1
OUTPUT_HDF                  : 1
DISCHARGES                  : 1

TIME_SERIE                  : 1
<endproperty>
```

### Run Mohid and inspect the results

If all is correctly configured, and all the input files are double-checked, we can hopefully start running our model!

#### Run Mohid

- Select *Tools* from menu, then *Launch Mohid...*

- Select the *Run\_1*, then select *Run Mohid*, then select *Launch after Ok*, then click on *Ok*

- Supervise the screen-outputs to check for any error, or to simply read the estimated end time of simulation

#### Inspect the results

Once the model has finished running, (or while the model is running, after it has created a few HDF5 outputs), you will want to switch in analysis mode and graphically inspect the results to see if all went well (or is going well). We will proceed to visualize the instantaneous velocity as a vector field, and the instantaneous water level as a colormap field. The proceedings are exactly the same for any other variables.

- Click on HDF file, then click on Grid HDF group

- Map the GridY, GridX, VerticalZ and Mapping with the Group items Latitude, Longitude, VerticalZ and OpenPoints, respectively.

- Map the time HDF Group to the Time selection button

- Map the Velocity U, Velocity V and Water level HDF results group's items to the Vector X, Vector Y and Color selection buttons, respectively

- From the *Action menu*, select *Run OpenGL*

- Select *settings*

- Select *Color*

- Define the minimum water level (m) and the maximum water level (m). Then click on *Ok*.

- Select *Vectors*

- Define the vector *Scale* to 0.2 and the vector stepping along X and along Y to 6. Then click on *Ok*.

If all went well, you should see results with patterns and values as expected for the domain.

- Here is the result of discharging water in a river

- Now, in order to render the results with the salinity map, you need to close the current post-processor rendering

- delete all the color layer instants of waterlevel,

- open the waterproperties hdf5 file,

- load the salinity property as a color layer,

- re-run in OpenGL mode

- Edit the color layer settings to a minimum of 0 and a maximum of 36.5 (as salinity in this case varies between these extrema).

- Once you apply with the settings, you can now render each frame and watch the freshwater discharge flow along your simulation. Thus, here is the result of an accidental freshwater discharge

## Next steps

- How to use the Blumberg-Kantha radiation method
- [How to couple atmosphere forcing](../../water/coupling-water-atmosphere-user-manual.md)

## Previous steps

- [Software to run Mohid](../../troubleshooting/software-to-run-mohid.md)
- [How to create a model gridded bathymetry?](how-to-create-a-model-gridded-bathymetry.md)
- [How to generate tide for Mohid?](how-to-generate-tide-for-mohid.md)
- [How to configure a MOHID 2D simulation with tide?](how-to-configure-a-2d-model-forced-with-tide-with-mohid.md)
- [How to create MOHID time series outputs?](how-to-create-mohid-timeseries-outputs.md)

## About the screen capturing tool

The screen capturing tool used for this tutorial was the [Jing Project](http://www.jingproject.com/).
