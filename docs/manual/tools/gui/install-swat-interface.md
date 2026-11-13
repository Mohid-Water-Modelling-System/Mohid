# Install SWAT interface

## Introduction

SWAT interface runs as an extension of ArcView or [MapWindow](http://www.mapwindow.org/). This means that to run SWAT interface you need to install any of the softwares.

## ArcView Installation Steps

SWAT interface also needs Arcview Spatial Analyst functions, for the Watershed Delineation and Landuse and Soil Definition components. You can read more details on SWAT interface on [AVSWAT](http://www.brc.tamus.edu/swat/avswat.html).

- Install Arcview 3.1 or higher (\\\\Kepler\\Software\\AVSWAT-X\\Step\_1 **or get it from** ArcView CD)
- Install Arcview extensions: 3D Analyst, Network Analyst, Spatial Analyst (\\\\Kepler\\Software\\AVSWAT-X\\Step\_2 **or get it from** ArcView CD)
- Install AvSWAT-X (\\\\Kepler\\Software\\AVSWAT-X\\Step\_3 **or get it from** the site [AVSWAT](http://www.brc.tamus.edu/swat/avswat.html))
    -   Unzip the avswatx.zip
    -   Double click the **Start.apr** project file or use ArcView to open it (**do NOT run the instalation file AVSWATX.msi**).
    -   The opening project starts the AVSWATX Install Wizard driving through the remainder of the installation process.
    -   When prompt for the target of this application, select preferable the "Only for me option"
    -   Once the installation is complete, open a new-empty ArcView session (do **NOT** use a working project with Tables or Views in it), select the menu File, menu item Extensions, select “AVSWATX Extendable”, and press the OK button. **Note: Avoid setting the extension as default (i.e. do NOT check the Make Default Option).**
    -   The main interface dialog will be displayed. Continue creating a new project (New Project ) or opening a previous one (Open Project ) identified by the file extension .avsx.

## MapWindow Installation Steps

[MapWindow](http://www.mapwindow.org/) is an open source Geographical Information System that is free to download. The SWAT version for this GIS is called MWSWAT and can be downloaded from the [Waterbase](http://www.waterbase.org/) website along with other sofware.

- Install the 32-bit MapWindow version (\\\\kepler\\Software\\MapWindow or **get it from** the site [MapWindow](http://www.mapwindow.org/))
- Install MWSWAT (\\\\kepler\\Software\\MWSWAT or **get it from** the site [Waterbase](http://www.waterbase.org/))
- Install the SWAT Editor 2005 v.2.1.5 (\\\\kepler\\Software\\MWSWAT, also available from the [SWAT website](http://swatmodel.tamu.edu/)).

From the [Waterbase](http://www.waterbase.org/) site also can be downloaded two more pieces of software:

- SWATPlot is a tool to extract output from SWAT runs. Together with its companion tool SWATGraph you can quickly visualise your data.
- MWAGNPS is an interface to the AGricultural Non-Point Source pollution model (AGNPS). It is a distributed model that simulates agricultural watersheds for a single storm event assuming uniform precipitation patterns. It performs estimates of runoff quality with primary emphasis on nutrients, pesticides and sediments.

```
IMPORTANT ISSUES

Before using, the next files have to be replaced with the ones in \\kepler\Software\MWSWAT\Fix.zip

- MapWinGeoProc.dll located in the MapWindow installation folder (i.e. C:\Program Files\MapWindow)
- CreateHRU.dll     located in the subfolder ..\Plugins\MWSWAT
- MWSWAT.dll        located in the subfolder ..\Plugins\MWSWAT

For 64 Bit machines, the SWAT editor needs to be forced to run as a 32 Bit application. In order to
achieve it, copy in the folder where the swatEditor.exe is located the program corflags.exe and run
the following line:

corflags swatEditor.exe /32BIT+

The corflags.exe file it is included in the Microsoft® .NET Framework and the version I used was
v3.5.21022.8. This will solve the problem of compatibility with the database OLEDB.4.0.
```
