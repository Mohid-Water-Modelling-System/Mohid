# Operational Modelling

## Overview

This wiki intends to explain how operational models are set in the Mohid Maretec framework.

Before intending to put any model in operational, the model has to be strongly tested.

## Prerequisites

Mohid and Matlab have to be installed. Mohid version and dll needed in case of being 32 or 64 bits.

## Folders Needed

In the folder where the model is going to be set up, it is needed to define the following generic folders:

- Backup
- GeneralData
- Logs
- Plots
- Work (only if atmospheric interpolation is needed)

Also a folder for the first level model will be defined, and it will contain the following subfolders:

- data
- exe
- figures
- res
- a folder with the name of the submodel in case it exist, this folder must replicate this same estructure

## GeneralData folder

This folder must contains the following subfolders:

- Bathymetry: Here it should be contained the bathymetries used by the model
- BoundaryConditions: Folders for the boundary conditions should be defined here (Atmosphere, Hydrodynamics, etc.). In them, the temporary files for each run should be copied.
- Exe: In this folder should be located the MohidWater.exe file, this is the only name that is accepted for the Mohid executable file.

Also should contain any other subfolders needed for input/output definitions as discharges, timeseries or boxes files

## Input Data

For each of the different modelling levels, a set of typical input data for each of the modules have to be selected (model.dat, hydrodynamic.dat, waterproperties.dat, turbulence.dat, etc.) also a set of hdf5 final files (.fin) for the hydrodynamic, GOTM and waterproperties in case they exist.

## Atmospheric forcing

If the model have interpolated atmospheric forcing data for the particular gridded data, a location should be provided. In the case that it need to be calculated, a folder named work should be included with the configuration files and programs. Atmospherical forcing conditions can be provided by different sources: GFS, MM5 and WRF.

### MM5

Two consecutive programs would be executed here, a first program GetWindForecast.exe that will go to find the files corresponding for the selected dates and that automatically will glue the files using a second software ConvertToHDF5Action.exe. The input file for the first software is described below.

```
GetMM5Forecast.dat
MM5_PATH                     : \\ftpserver\ftp.mohid.com\LocalUser\meteoIST\mm5_6h\
WORK_PATH                    : D:\Aplica\PCOMS\Work\MM5\
OUTPUT_FOLDER                : D:\Aplica\PCOMS\Work\MM5\History\
BACKUP_FOLDER                : \\datacenter\alexandria4\modelos\Mohid_IST\PCOMS\Portugal\Meteo\
MOHID_GRID                   : Portugal.dat
MAIL_SERVER                  : mail.ist.utl.pt
MAIL_SENDER_NAME             : Portuguese Coast Operational Modelling System
MAIL_SENDER_ADDRESS          :
CHECK_NUMBER_OF_FILES        : 1
NUMBER_OF_FILES              : 13
<begin_out_email_list>
<end_out_email_list>
```

Once the first software has performed the find and glue task, the ConvertToHDF5Action.exe will generate a input file for the desired period of time and will interpolate for the selected grid.

```
ConvertToHDF5Action.dat
<begin_file>
ACTION                        : PATCH HDF5 FILES
TYPE_OF_INTERPOLATION         :  3
START                         : 2011 4 19 0 0 0
END                           : 2011 4 22 0 0 0
<<begin_father>>
FATHER_FILENAME               : MM5_D2_Portugal.hdf5
FATHER_GRID_FILENAME          : D2.dat
LEVEL                         :  2
<<end_father>>
<<begin_father>>
FATHER_FILENAME               : MM5_D3_Portugal.hdf5
FATHER_GRID_FILENAME          : D3.dat
LEVEL                         :  3
<<end_father>>
OUTPUTFILENAME                : MM5_Portugal_2011-04-19_2011-04-22.hdf5
NEW_GRID_FILENAME             : Portugal.dat
<<BeginFields>>
air temperature
downward long wave radiation
wind velocity X
wind velocity Y
wind modulus
relative humidity
solar radiation
mean sea level pressure
<<EndFields>>
<end_file>
```

### WRF

Two consecutive programs would be executed here, a first program GetMeteoForecast.exe that will go to find the files corresponding for the selected dates and that automatically will glue the files using a second software ConvertToHDF5Action.exe. The input file for the first software is described below.

```
GetMeteoForecast.dat
METEO_MODEL                              : WRF
METEO_MODEL_PATH                         : \\DATACENTER\alexandria3\modelos\Meteo_IST\wrf_6h
WORK_PATH                                : D:\Aplica\OPTagus3D\Work\WRF\
OUTPUT_FOLDER                            : D:\Aplica\OPTagus3D\Work\WRF\History\
BACKUP_FOLDER                            : \\datacenter\alexandria4\modelos\Mohid_IST\OPTagus3D\Tagus\Meteo\
MOHID_GRID                               : D:\Aplica\OPTagus3D\GeneralData\Bathymetry\TagusRibeiras.dat
MOHID_GRID_LABEL                         : Tagus3D
MAIL_SERVER                              : mail.ist.utl.pt
MAIL_SENDER_NAME                         : OPTagus3D
MAIL_SENDER_ADDRESS                      :
!So por ser hindcast
START                                    : 2011 2 6 0 0 0
END                                      : 2011 2 7 0 0 0
!NUMBER_OF_DAYS                          : 1
!CHECK_NUMBER_OF_FILES                   : 1
!NUMBER_OF_FILES                         : 5
<begin_domain_list>
WRF_3km Lisbon 3
<end_domain_list>
<begin_out_email_list>
<end_out_email_list>
```

Once the first software has performed the find and glue task, the ConvertToHDF5Action.exe will generate a input file for the desired period of time and will interpolate for the selected grid.

```
ConvertToHDF5Action.dat
<begin_file>
ACTION                        : INTERPOLATE GRIDS
TYPE_OF_INTERPOLATION         :  3
START                         : 2011 2 6 0 0 0
END                           : 2011 2 7 0 0 0
<<begin_father>>
FATHER_FILENAME               : WRF_WRF_3km_Lisbon.hdf5
FATHER_GRID_FILENAME          : WRF_3km.dat
LEVEL                         : 3
<end_father>>
OUTPUTFILENAME                : WRF_Tagus3D_2011-02-06_2011-02-07.hdf5
NEW_GRID_FILENAME             : D:\Aplica\OPTagus3D\GeneralData\Bathymetry\TagusRibeiras.dat
<<BeginFields>>
air temperature
downward long wave radiation
wind velocity X
wind velocity Y
wind modulus
relative humidity
solar radiation
mean sea level pressure
<<EndFields>>
<end_file>
```

## Running Mohid

In each of the exe folder a nomfich.dat should be placed and modified with the needed paths.

```
nomfich.dat
IN_BATIM                      : ..\..\GeneralData\Bathymetry\Portugal_Window.dat
ROOT                          : D:\Aplica\OPTagus3D\Portugal_Window\res\
ROOT_SRT                      : D:\Aplica\OPTagus3D\Portugal_Window\res\Run1\
IN_MODEL                      : D:\Aplica\OPTagus3D\Portugal_Window\data\Model_1.dat
SURF_DAT                      : D:\Aplica\OPTagus3D\Portugal_Window\data\Atmosphere_1.dat
SURF_HDF                      : D:\Aplica\OPTagus3D\Portugal_Window\res\Atmosphere_1.hdf
DOMAIN                        : D:\Aplica\OPTagus3D\Portugal_Window\data\Geometry_1.dat
AIRW_DAT                      : D:\Aplica\OPTagus3D\Portugal_Window\data\InterfaceWaterAir_1.dat
AIRW_HDF                      : D:\Aplica\OPTagus3D\Portugal_Window\res\InterfaceWaterAir_1.hdf
BOT_DAT                       : D:\Aplica\OPTagus3D\Portugal_Window\data\InterfaceSedimentWater_1.dat
BOT_HDF                       : D:\Aplica\OPTagus3D\Portugal_Window\res\InterfaceSedimentWater_1.hdf
IN_DAD3D                      : D:\Aplica\OPTagus3D\Portugal_Window\data\Hydrodynamic_1.dat
OUT_DESF                      : D:\Aplica\OPTagus3D\Portugal_Window\res\Hydrodynamic_1.hdf
IN_TURB                       : D:\Aplica\OPTagus3D\Portugal_Window\data\Turbulence_1.dat
TURB_HDF                      : D:\Aplica\OPTagus3D\Portugal_Window\res\Turbulence_1.hdf
DISPQUAL                      : D:\Aplica\OPTagus3D\Portugal_Window\data\WaterProperties_1.dat
EUL_HDF                       : D:\Aplica\OPTagus3D\Portugal_Window\res\WaterProperties_1.hdf
OUTWATCH                      : D:\Aplica\OPTagus3D\Portugal_Window\res\Outwatch_1.txt
DT_LOG                        : D:\Aplica\OPTagus3D\Portugal_Window\res\DTLog_1.txt
```

and also a Tree.dat file with the locations of the respective submodels to be run, in case there is a problem and only one of them is to be ran, it is enough to include an empty line between the different domains.

```
Tree.dat
Automatic Generated Tree File
by MOHID GUI
+D:\Aplica\OPTagus3d\Portugal_Window\exe
++D:\Aplica\OPTagus3d\Portugal_Window\3DTagusWQ\exe
```

## Hindcasting

In this case it will be used the executable AutoTagusROFI.exe and the configuration file is TagusROFI-Auto.dat, in this file should be indicated the location of the specific configuration input files and the minimum and maximum time that the model should run, to stop the operational procedure in seconds.

```
TagusROFI-Auto.dat
START                      : 2011 01 08 0 0 0
END                        : 2011 03 01 0 0 0
TAGUSROFI_INPUT_FILE       : \\turbo\aplica\OPSAndre\TagusROFI.dat
TAGUSROFIPLOTS_INPUT_FILE  : \\turbo\aplica\OPSAndre\Plots\MatlabPlots.dat
TAGUSROFI_METEO_INPUT_FILE : \\turbo\aplica\OPSAndre\Work\WRF\GetMeteoForecast.dat
TAGUSROFI_BATCH            : \\turbo\aplica\OPSAndre\TagusROFI.bat
PLOTS_ONLY                 : 0
MIN_TIME                   : 200
MAX_TIME                   : 30000
TagusROFI.dat
MAIN_PATH                                : D:\Aplica\OPTagus3D\
MAIL_SERVER                              : mail.ist.utl.pt
MAIL_SENDER_NAME                         : OPTagus3D
MAIL_SENDER_ADDRESS                      :
START                                    : 2011 2 12 0 0 0
END                                      : 2011 2 13 0 0 0
OPENMP                                   : 1
OPENMP_NUM_THREADS                       : 6
<begin_out_email_list>
List of e-mail addresses to be contacted in the format name <email address>
<end_out_email_list>
<begin_model>
NAME                                     : Portugal_Window
PATH                                     : Portugal_Window\
MODEL_HAS_OBC                            : 1
OBC_WORK_PATH                            :
OBC_SIMULATED_DAYS                       : 1 (By default 3)
MODEL_HAS_METEO                          : 0
GRID_FILE                                :
MODEL_RUN_ID                             : 1
DT                                       : 180
BACKUP_PATH                              : D:\Aplica\OPTagus3D\BackUp\Portugal_Window\
STORAGE_PATH                             : \\DATACENTER\alexandria4\modelos\Mohid_IST\OPTagus3D\Portugal_Window\
MODEL_HAS_WATERPROPERTIES                : 1
MODEL_HAS_SURFACE_HDF                    : 0
MODEL_HAS_GOTM                           : 0
MODEL_HAS_OUTPUT_WINDOW                  : 0
MODEL_HAS_SOLUTION_FROM_FILE             : 1
<end_model>
<begin_model>
NAME                                     : 3DTagusWQ
PATH                                     : Portugal_Window\3DTagusWQ\
MODEL_HAS_OBC                            : 0
MODEL_HAS_METEO                          : 1
METEO_WORK_PATH                          : D:\Aplica\OPTagus3D\Work\WRF\History\
GRID_FILE                                : Portugal_Window\GeneralData\Bathymetry\TagusRibeiras.dat
MODEL_RUN_ID                             : 1
DT                                       : 12
BACKUP_PATH                              : D:\Aplica\OPTagus3D\BackUp\3DTagusWQ\
STORAGE_PATH                             :
MODEL_HAS_WATERPROPERTIES                : 1
MODEL_HAS_INTERFACESEDIMENTWATER         : 1
MODEL_HAS_SURFACE_HDF                    : 1
MODEL_HAS_GOTM                           : 1
MODEL_HAS_OUTPUT_WINDOW                  : 0
<end_model>
```

## Forecasting

In this case the executable to be used will be TagusROFI-OP.exe

## Running Manually

The batch file to be used will be MOHID\_run - manual.bat and should be edited according to the project, as an example:

```
cd "D:\Aplica\OPSAndre\Portugal_Window\exe"
D:
"D:\Aplica\OPSAndre\GeneralData\Exe\MohidWater.exe"
pause
```

## Backup and Storage

A folder containing the name of each submodels should be created in the backup folder, i.e.

```
D:\Aplica\OPSAndre\BackUp\Portugal_Window\
```

and also a location in the dataserver should be created for the storage, i.e.

```
\\datacenter\alexandria4\modelos\Mohid_IST\OPSAndre\Portugal_Window\
```

## Plots

It is needed to have matlab and m\_map previously installed. The matlab can be found in the computer and the m\_map can be downloaded from [http://www.eos.ubc.ca/~rich/map.html](http://www.eos.ubc.ca/~rich/map.html) and once unzipped it should be placed in the toolbox folder where the Matlab has been installed.

```
To install M_Map as a toolbox accessible to all users this, unpack the files into
$MATLAB/toolbox/m_map,
add that directory to the list defined in $MATLAB/toolbox/local/pathdef.m,
matlabroot,'\toolbox\m_map;', ...
and update the cache file using
rehash toolboxcache
You should restart before this take effect.
```

In order to plot automatically the mohid results, it will be needed the software MatlabPlot.exe that is configured by the file MatlabPlots.dat and executed by the file make\_plots.bat.

```
make_plots.bat
matlab -nodesktop -nosplash -r mat_main_plot -minimize -wait
```

```
MatlabPlots.dat
WORK_PATH                                : D:\aplica\OPTagus3D_2010\Plots\
MAIL_SERVER                              : mail.ist.utl.pt
MAIL_SENDER_NAME                         : OPTagus3D_2010 - Plots
MAIL_SENDER_ADDRESS                      :
START                                    : 2011 2 14 0 0 0
END                                      : 2011 2 15 0 0 0
<begin_out_email_list>
name <email address>
<end_out_email_list>
SUCCESS_FILE                             : D:\aplica\OPTagus3D_2010\MOHID_SUCCESS.dat
MOHID_PATH                               : D:\aplica\OPTagus3D_2010\BackUp\3DTagusWQ\Results_HDF\
PUBLISH_PATH                             : D:\aplica\OPTagus3D_2010\Plots\Publish\
PLOT_TEMP_SALT_VEL                       : 1
PLOT_BIO                                 : 1
BIO_MODEL                                : waterquality
ADD_DAYS                                 : 1
```

Also some matlab files have to be edited and modified for each application,

```
mat_setup_plots.m
mat_setup_dir = 'D:\Aplica\OPTagus3D_2010\Plots';
```

```
mat_config_coastline.m
coastline_file        = 'tagus_coastlineh.mat';
save_coastline        = false;
coastline_res         = 'i';        %'l'=low, 'i'=intermediate, 'h'=high, 'f'=full;
```

To set the limits and the colour scale of each property for instance,

```
mat_config_contour.m
sst_contourlimits     = [12 26];
sst_contour_scale     = [12 13 14 15 16 17 18 19 20 21 22 23 24 25 26];
```

## Validation

to validate with OSTIA SST images first it should be installed the mexcdf package for matlab that allows operations with netcdf.

- [mexnc R2010b and beyond](http://mexcdf.sourceforge.net/downloads/)

matlabroot,'\\toolbox\\mexcdf\\mexnc;', ... matlabroot,'\\toolbox\\mexcdf\\snctools;', ...

### Sections

```
Variable List
1 =  temperature'
2 =  salinity'
3 =  density'
4 =  Density - SigmaT'
_______________________________'
Life Model'
_______________________________'
7 =  diatoms chlorophyll'
8 =  autotrophic flagellates chlorophyll'
9 =  nitrate'
10 = silicate acid'
11 = inorganic phosphorus'
12 = oxygen'
13 = diatoms carbon'
14 = autotrophic flagellates carbon'
15 = semi-labile dissolved organic carbon'
16 = labile dissolved organic carbon'
17 = particulate organic carbon'
18 = microzooplankton carbon'
19 = heterotrophic bacteria carbon'
20 = biogenic silica'
26 = short wave solar radiation'
27 = clorofila_total (diatoms + flagellates)'
28 = carbono_total (diatoms + flagellates)'
_______________________________'
Hydrodynamic'
_______________________________'
22 = velocity U'
23 = velocity V'
24 = velocity W'
25 = velocity modulus'
_______________________________'
Energy'
_______________________________'
29 = power density'
_______________________________'
Water Quality Model'
_______________________________'
5 =  cohesive sediment'
6 =  phytoplankton'
12 = oxygen'
26 = short wave solar radiation'
30 = nitrate'
31 = nitrite'
32 = ammonia'
33 = dissolved non-refractory organic nitrogen'
34 = dissolved refractory organic nitrogen'
35 = particulate organic nitrogen'
36 = dissolved refractory organic phosphorus'
37 = particulate organic phosphorus'
38 = zooplankton'
_______________________________'
ATMOSPHERE Model'
_______________________________'
39 = wind velocity X'
40 = wind velocity Y'
41 = solar radiation'
42 = precipitation'
43 = atmospheric pressure'
44 = air temperature'
45 = relative humidity'
46 = cloud cover'
47 = mean sea level pressure'
48 = wind modulus'
49 = wind direction'
50 = specific humidity'
_______________________________'
Water air interface'
_______________________________'
51 = latent heat'
52 = sensible heat'
53 = evaporation'
54 = net long wave radiation'
55 = upward long wave radiation'
56 = downward long wave radiation'
57 = short wave solar radiation'
58 = long wave solar radiation'
59 = short wave solar radiation extinction'
60 = long wave solar radiation extinction'
61 = oxygen flux'
62 = wind shear velocity'
63 = surface radiation'
64 = wind stress X'
65 = wind stress Y'
66 = surface water flux'
67 = non solar flux'
69 = albedo'
```

## Publishing online

## Operational Models

PCOMS running in mar computer

BIOPCOMS WQ running in pcoms-bio computer

PCOMS-BIO2 Life running in pcoms-bio computer

OPTagus3D running in Tagus computer

OPTagus3D\_2010 running in Turbo computer

## To be done

- Include a delay time of 6 hours to test again

- To do independent forecast pcoms and pcoms-bio

- Harmonize tools
