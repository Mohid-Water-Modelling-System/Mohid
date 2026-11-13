# How to generate tide for Mohid?

Once you [created a gridded bathymetry](how-to-create-a-model-gridded-bathymetry.md), before you run the model simulation, you may want to generate the tidal gauges that will generate the tide at the open-boundaries during the simulation.

## Step 1 - Create the tidal gauges points in the GIS

### Create a new points layer in GIS

- Make sure that you already have created your bathymetric gridded data file and that it's currently open in your GIS project. Now, simply add a new data item from the menu:



- Select *XYZ* type, then click on *Browse*

- Choose a name for the tidal gauge points, then hit save:

- Hit the Ok button

### Edit the points layer

- Type *4* in the size text-box, and then select settings:

- Select option *Constant*, then click on Pick Color (a new window pops up), then select a bright color and hit Ok on the new window (it will make the new window close), then click on Ok, and the new configuration will be enabled.

- Close and comply with changes by hitting the Ok button

### Add the tidal gauges to the points layer

- Select the *TideGaugesTest.xyz* points layer, then select the add a point command:

- Hit Ok ...

- Click to add points all along the open boundaries, roughly 10 cells wide apart:

- Keep adding those points alternatively on each side of the open-boundary

- Final look:

### Save the tidal gauges points layer

- From the *Data Items* menu, select the *Save All* option

## Step 2 - Generate the tidal gauges from FES2004

In order to generate tidal gauges from FES2004, you will need to obtain the data file and the program that performs the extraction into MOHID tidal-gauge format. [Here](../fes2004.md) you will find details about that.

- Get the [mohid-tide](http://maretec.mohid.com/PublicData/products/Software/Mohid-Tide_v0.4.zip) package and extract it somewhere on your disk;

- Transfer the xyz file to the mohid-tide folder. *WARNING: Remove any blank lines from the xyz file!*

- Edit the run-tide.bat accordingly to the following example:

```
>mohid-tide ..\fes2004data\tide.nc TideGaugesTest.xyz output.dat 0 2.08
```

- Run the batch file. This will create the tidal gauge file.

- Verify if your tidal gauge file is correctly formed. Edit it and inspect.

### That's it!

If all went well, then you have a valid gridded bathymetry for the domain and you have tidal gauges extracted from the [Fes2004](../fes2004.md) solution that will generate tide to the model. You're ok to start building a project using [Mohid GUI](mohid-gui.md).

## Next steps

- Read the [Mohid GUI](mohid-gui.md) to create a new simulation project then,
- [How to configure a 2D model forced with tide with MOHID](how-to-configure-a-2d-model-forced-with-tide-with-mohid.md)

## About the screen capturing tool

The screen capturing tool used for this tutorial was the [Jing Project](http://www.jingproject.com/).
