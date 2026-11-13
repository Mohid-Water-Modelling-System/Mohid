# MatLab Tools

This is a colletion of simple scripts to enable the analysis of MOHID's results in [MatLab](../../archive/linux/matlab.md). Many of these are useful for code reference, for implementing your own scripts.

These tools are available to [MOHID](../../overview/mohid.md) users and developers in the [SourceSafe](../../building/windows/sourcesafe.md) database "OtherTools/MatLab".

## Scripts availiable

Basic scripts for **opening and reading** the contents of MOHID HDF5 files, there are the following scripts:

\* mh5grid.m : reads and outputs the Grid group from a file
\* mh5read.m : reads and outputs a specific dataset instant(or several) from a file
\* mh5info.m : reads and outputs some information on the file, as well as the number of instants.

For plotting MOHID HDF5 files, usually saving as png file:

\* mh5ts\_plot.m : plots a timeserie of a property in model cell(i,j,z)
\* can plot a station timeserie and the model's for comparision(using waterbase\_ts.mat)
\* can plot several properties in one plot using morewaterpropz
\* mh5xy\_plot.m : plots an x-y graph at layer z for several instants within a file
\* can do a gif of the images generated
\* mh5zt\_plot.m : plots a z-t graph, creating a vertical profile versus time graph.

Converting a [netcdf](../../archive/linux/netcdf.md) file to a MOHID HDF5 file:

\* nc2mh5\_convert.m

Model's timeserie analysis, for e.g., box files:

\* ts\_BoxAnalyser.m : plots each box's properties versus time
\* ts\_PhytoLimFactors.m : plots each box's phytoplankton limitation
\* ts\_BoxFluxAnalyser.m : for reference, plots the fluxes between boxes. Can use a box image to plot arrows.

Extra scripts for plotting(may be usefull for code reference):

```
 * mh5strat_analyser.m : plots(x-y) the difference between two layers.
        Usefull to check oxygen depletion or stratification.
 * mh5temp_comp.m : compares the temperature field with the assimilation temperature
```
