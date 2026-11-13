# Mohid Swan

In this wiki it would be commented how to implement the Swan model into the Mohid modelling software.

## Software Needed

In addition to the Mohid model it will be needed the following software to achive the goal.

- [Swan Model](http://swanmodel.sourceforge.net/)
- [Matlab](http://www.themathworks.com/)
- [SIGMA - A SWAN Graphical Interface in MATLAB](http://users.coastal.ufl.edu/~alex/downloads.html) (Optional)

For Maretequers, the software (with configuration examples) and user manuals can be found internally here

```
\\neptuno\Software\Swan
```

To get an example of the whole data pre-processing, look internally, here

```
\\neptuno\Software\MOHID\Mohid Swan
```

## Bathymetry

To create a SWAN bathymetry, one may convert a MOHID bathymetry (grid data file) using the matlab script found in the *\\\\neptuno\\Software\\MOHID\\Mohid Swan\\MohidGriddataToSWAN* folder. In Matlab, select the script folder as the working directory then type at the matlab prompt:

```
> MohidBatim2SWAN.m
```

Then follow the GUI instructions. It should create a SWAN bathymetry file (.sty).

The bathymetry should be visualised using notepad or a software called sigma that needs to be localed in the root directory c:\\

Sigma -> Task File -> Show so far -> Linha CGWID

## Boundary Conditions

The Swan model would need of wind, tide and deep-Ocean swell as boundary conditions.

### Wind

Wind vectors should also be provided from other model results (i.e. MOHID, MM5 or WRF) and then converted to the SWAN format using [ConvertToHDF5](converttohdf5.md) using the option "ConvertHDF5ToSWANorMOHID" that will provide a hourly file with values for x and y directions, be aware of the file size.

One way to extract winds to force the SWAN, is to use the wind velocity present in the Atmosphere\_x.hdf5 file, coming from a 2D MOHID simulation pre-made for the study region.

An alternative way is to use MM5 HDF5 Mohid results and then [glue](converttohdf5.md) and [interpolate](converttohdf5.md) them to the MOHID model(s) grid(s) (which is(are) the same as SWAN's). To perform standard [glue](converttohdf5.md) and [interpolate](converttohdf5.md) operations with [ConvertToHDF5](converttohdf5.md), please refer to the [ConvertToHDF5](converttohdf5.md) wiki ...

#### ConvertToHDF5 configuration

Here's the [ConvertToHdf5](converttohdf5.md) configuration file *ConvertToHDF5Action.dat* for the action *CONVERT FROM HDF5 TO SWAN OR MOHID* to use (notice that there are two beginfile/endfile blocks, one for X direction, and another for the Y direction):

```
<begin_file>

ACTION                         : CONVERT FROM HDF5 TO SWAN OR MOHID

INPUTFILENAME                  : Atmosphere_1.hdf5
INPUT_GRID_FILENAME            : ../../Batim/MOHID/BatimSO.dat

<<beginproperty>>
wind velocity X
<<endproperty>>

START                          : 2011 04 04  0 0 0
END                            : 2011 04 07  0 0 0
OUTPUT_TIME                    : 0  3600
OUTPUT_OPTION                  : 1
<end_file>

<begin_file>

ACTION                         : CONVERT FROM HDF5 TO SWAN OR MOHID

INPUTFILENAME                  : Atmosphere_1.hdf5
INPUT_GRID_FILENAME            : ../../Batim/MOHID/BatimSO.dat

<<beginproperty>>
wind velocity Y
<<endproperty>>

START                          : 2011 04 04  0 0 0
END                            : 2011 04 07  0 0 0
OUTPUT_TIME                    : 0  3600
OUTPUT_OPTION                  : 1
<end_file>
```

Note that this configuration will produce hourly outputs. This will generate one .txt file per time instant.

#### Criating a wind.ini SWAN input file

The SWAN input file is a *wind.ini* file that contains the listing of the wind .txt files created in the subsection above. Example:

```
Wind/windX_20110404000000.txt
Wind/windX_20110404010000.txt
Wind/windX_20110404020000.txt
Wind/windX_20110404030000.txt
Wind/windX_20110404040000.txt
Wind/windX_20110404050000.txt
Wind/windY_20110404000000.txt
Wind/windY_20110404010000.txt
Wind/windY_20110404020000.txt
Wind/windY_20110404030000.txt
Wind/windY_20110404040000.txt
Wind/windY_20110404050000.txt
```

Notice how there are first the X-direction files and then, after, come the Y-direction files.

#### Troubleshoot

**Q: The conversion won't work or the SWAN model won't run because of the winds.**

A: The Atmosphere\_x.hdf5 of the MOHID results needs to be edited and changed (using [HDFView](README.md#hdfview#hdfview)). 1 - The *MappingPoints2D* data needs to be renamed to *WaterPoints* in the *GRID/* group. 2 - The first instant of *Time* is could be wrong; it could display an offset of one dt. The time instant data should be edited and the offset should be removed. Example: Once, the first time instant in the HDF5 atmosphere\_x file would read *2006 10 15 17 00 20* when it should read instead *2006 10 15 17 00 00*.

**Q: The SWAN won't run because it can't find a wind text file, whose name in the log is a few characters short.**

A: The path and filenames to the wind text files found in the *wind.ini* file exceed the maximum allowed number of characters in a single line. Try reducing the number of characters, by renaming the *wind.ini* file and/or by renaming its files.

### Tide

The tide sould be obtained from other models (i.e. PCOMS), or calculated (i.e. [Tideprev](tide-preview.md)) and interpolated to the user domain grid.

For the interpolation method between the source grid and the user domain grid, the user is recommended to use the [ConvertToHdf5](converttohdf5.md) tool with the *INTERPOLATE* action. [Please refer to this wiki](converttohdf5.md).

Alternatively, the user domain model could be pre-run with tide only and outputs be given with the same time-period and same frequency as the wind's. Thus the user would only need to use the *water level* results from the hydrodynamic\_x.hdf5 file.

Then, the tide should be converted into SWAN input file for tide, using the *CONVERT FROM HDF5 TO SWAN OR MOHID* action method of the [ConvertToHDF5](converttohdf5.md) tool. Here's the input configuration to copy/paste/edit in the *ConvertToHdf5Action.dat* file:

```
<begin_file>

ACTION                         : CONVERT FROM HDF5 TO SWAN OR MOHID

INPUTFILENAME                  : Hydrodynamic_x.hdf5
INPUT_GRID_FILENAME            : Batim/MOHID/BatimSO.dat

<<beginproperty>>
water level
<<endproperty>>

START                          : 2011 04 04  0 0 0
END                            : 2011 04 05  0 0 0
OUTPUT_TIME                    : 0  3600
OUTPUT_OPTION                  : 1
<end_file>
```

#### Criating a WaterLevel.ini SWAN input file

The SWAN input file is a WaterLevel.ini file that contains the listing of the waterlevel .txt files created in the subsection above. Example:

```
Wlevel/WL_20110404000000.txt
Wlevel/WL_20110404010000.txt
Wlevel/WL_20110404020000.txt
Wlevel/WL_20110404030000.txt
Wlevel/WL_20110404040000.txt
Wlevel/WL_20110404050000.txt
```

### Swell

The swell can be read from a global wave model, such as [WaveWatch III](http://polar.ncep.noaa.gov/waves/wavewatch/).

Alternatively, a constant swell can be imposed at each of the lateral open boundaries of the user domain (North, South, East or West sides). This is defined and examplified in the SWAN configuration file sample below in this wiki article.

The North position must be defined (ex: 90º) and the orientation is direct (or anti-clockwise).

## SWAN Model Input

Please refer to sample input file in the subsection below for a reference input file. This section will only describe the parts that the user NEEDS to edit. For more information, please refer to the [SWAN user manual](http://citg.tudelft.nl/index.php?id=17787&L=1).

```
Sect 3
CGRID is the origin, starting from the lower left corner
0   -> Rotation (0 no inclination)
1.3 -> Grid step x in degrees
2   -> Grid step y in degrees
130 -> cells -1 in X direction
200 -> cells -1 in Y direction
CIR -> Whole circle of directions
36  -> Number of directions
0.04 -> Minimum frequency
0.5 -> Maximum frequency
50 -> Number of frequencies -1 (in total then 51)
```

```
Sect 4
Equal to Sect 3 but define the grid resolution instead of the number of cells, for each direction, (ex: dx = 0.01 = dy)
wlevel.in
wind.in
```

```
Sect 5
WavewatchIII
BOUN (boundary) SIDE W CONST PAS 5 10 210 30
                                 H  T Dir Dispersion
```

```
Sect 6 Physics
Do not change
```

```
Sect 8 Time series
To obtain a window it should be edited
NGRID Indicating the lower left cell for the new grid, the size of the grid in x and y and the number of cells
SUBG Indicates the subgrid and in should be provided the lower left corner of the father grid
```

```
Sect 10
It makes the domain 4 times and when reach the 98% probability then stops.
Defines the initial and final date
It should be defined a day as the initial date
```

#### Excerpt from word manual

On the contrary of STWAVE there was not introduced any modifications to the original code and so the files format does not follows the “MOHID standard”. In this case it will be necessary to provide the bathymetry in SWAN format. Besides the bathymetric data file it is necessary to provide a file with the computational conditions usually named INPUT. The basic philosophy still lies in the definition ok keywords o specify the model options. For learning more about the different computational options please make use of the user manual. Here, for reference, it is presented an example of an INPUT file.

INPUT

$-------- Section 1. Start

PROJECT 'Quarteira' '01' Title and run ID SET 2 90 0.05 200 1 9.81 1 0 5 water level, north, depmin, max messages, warnings, gravity, output based n energy, nautical, power of high frequency tail

MODE TWOD mode twodimensional

$-------- Section 3. Model Description

CGRID REG 185000 -8000 0 29900 22450 299 449 CIR 36 0.05 1 50 Computational grid regular, Xorig, Yorig, X axis direction, domain length in X direction, domain length in Y direction, grid cells in X direction, grid cells in Y direction, spectral directions cover the full circle, number of meshes in -space (in this case it means that we are considering 10º divisons; 360/36), lower frequency, higher frequency, number of frequencies

$-------- Section 4. Input

INPGRID BOT REG 185000 -8000 0 299 449 100 50 EXCEPTION -99 READ BOT 1 'd/450x300\_Swan\_Ref.bty' 1 0 FREE WIND 0.1 0 Grid coordinates: referred to the bottom, regular, Xorig, Yorig, X axis direction, grid cells in X direction, grid cells in Y direction, DX, DY, Ignore points -99 Bathymetry file: data sign, filename, read mode, header lines, free format reading Wind characteristics: Velocity, direction

$-------- Section 5. Boundary/Initial Conditions

BOUN SHAPE JON 1 PEAK DSPR DEGREES BOUN SIDE E CONST PAR 7.5 16 45 30 BOUN SIDE S CONST PAR 7.5 16 45 30 BOUN SIDE W CONST PAR 7.5 16 45 30 Boundary conditions: Jonswap spectrum, gamma parameter, peak period, directional distribution, degrees Wave characteristics in the boundaries (E,S,W). spectral parameters (PAR), Hs, Tp, Dir, directional spreading.

$-------- Section 6. Physics

GEN3 KOM 2.36e-5 3.02e-3 QUAD 2 AGROW 0.0015 SETUP 0 BREAKING CONSTANT 1.0 0.73 FRICTION JON 0.067 TRIAD

Komen exponential growth, rate of whitecapping dissipation, wave steepness, Quadruplets, fully implicit computation, Cavaleri and Malanotte grow terms, Cavaleri and Malanotte coefficient Setup computation activated. 0 denotes that no constant is added. Breaking constant, alpha and gamma coefficients. Frction: Jonswap type, Jonswap coefficient Triad wave-wave interaction activated

$-------- Section 9. Output BLOCK 'COMPGRID' NOHEAD 'o/S' LAY 3 SETUP 1 BLOCK 'COMPGRID' NOHEAD 'o/H' LAY 3 HS 1 BLOCK 'COMPGRID' NOHEAD 'o/T' LAY 3 PER 1 BLOCK 'COMPGRID' NOHEAD 'o/F' LAY 3 FOR 1 BLOCK 'COMPGRID' NOHEAD 'o/D' LAY 3 DIR 1 BLOCK 'COMPGRID' NOHEAD 'o/W' LAY 3 WIND 1

These are the spatial distribution of the variables that will be written in the computational grid format. Noheaders will be added. Layout type 3, $ TABLE 'COMPGRID' NOHEAD 'o/data.tbl' XP YP DEP HS DIR SETUP FOR WIND TABLE 'COMPGRID' NOHEAD 'o/dataS.tbl' XP YP SETUP TABLE 'COMPGRID' NOHEAD 'o/dataF.tbl' XP YP FOR

Writes on a table the selected variables. A program named READSWANRES was set up to write these variables in a “MOHID standard format”

$-------- Section 10. Run NUM ACCUR 0.02 0.02 0.02 98 30 COMPUTE STOP

Numerics: criterion for terminating the iterative procedure: fraction of the wave height change, fraction of the average significant wave height, fraction the mean wave period, fraction of the points to met the conditions, maximum number of iterations for stationary computations

### Sample input file

This sample input file configures the SWAN to run forced with tide, wind and a constant swell of 5m and 210º in direction (with a 30º of dispersion) coming in from the South and West boundaries of the domain. It runs for 5 hours and demands hourly results in table format.

```
$-------- Section 1. Start
PROJECT  'SAndre'  '1'
SET 0.0 90 0.05 200 1 9.81 1025 0 5 NAUTICAL
MODE NONSTAT TWOD
COORD SPHE CCM
$
$-------- Section 2. General
$
$-------- Section 3. Model Description
CGRID REG -9.22 37.637 0 0.445 0.615 89 123 CIR 36 0.04 0.5 50
$
$-------- Section 4. Input
INPGRID BOTTOM REG -9.22 37.637 0 89 123 0.005 0.005  EXCEPTION -99.
READ BOT 1 'batim/Swan/BatimSO.bty' 1 0 FREE
$
INPGRID WLEVEL REG -9.22 37.637 0 89 123 0.005 0.005  EXCEPTION -999 &
 NONSTAT 20110404.000000 1 HR 20110404.050000
READ WLEV 1 SERIES 'WLevel/WLevel.Ini' 1 0 FREE
$
INPGRID WX REG -9.22 37.637 0 89 123 0.005 0.005  EXCEPTION -999 &
 NONSTAT 20110404.000000 1 HR 20110404.050000
$
INPGRID WY REG -9.22 37.637 0 89 123 0.005 0.005  EXCEPTION -999 &
 NONSTAT 20110404.000000 1 HR 20110404.050000
READ WIND 1 SERIES 'Wind/Wind.ini' 1 0 FREE
$
$
$-------- Section 5. Boundary/Initial Conditions
$BOUNDNEST3 WWIII 'WW3_Conditions/NBCWW3.nbc' OPEN
$
BOUN SHAPE JON 3.3 PEAK DSPR DEGREES
BOUN SIDE W CONST PAR 5 10 210 30
BOUN SIDE S CONST PAR 5 10 210 30
$
$-------- Section 6. Physics
GEN3 JANSSEN 4.5 0.5 AGROW 0.003
BREAKING CONSTANT
FRICTION JONSWAP
TRIAD
DIFFRACtion
$
$-------- Section 7. Numerics
$
$-------- Section 8. Output Sets
$
TABLE 'COMPGRID' NOHEAD 'Results/20110404.000000.tbl' DEP HS HSWELL PER RTP DIR PDIR FOR WIND  &
 WATLEV OUTPUT 20110404.000000 1 HR
$
 TABLE 'COMPGRID' NOHEAD 'Results/Time.dat' TIME &
OUTPUT 20110404.000000 1 HR
$
$-------- Section 10. Run
NUM ACCUR 0.02 0.02 0.02 98 30
COMPUTE NONSTAT 20110404.000000 20 MI 20110404.050000
STOP
```

## Running Swan

Open a SWAN shell (terminal) and execute the following commands:

```
> cd [path to my INPUT_new.nrw file]
> swanrun INPUT_new
```

## SWAN Model Output

Use the *ReadSwanRes.exe* to convert the .tbl file(s) into MOHID HDF5 format file. Edit the ReadSwanResults.dat input file accordingly from the following sample configuration:

```
../batim/MOHID/BatimSO.dat
20110404.000000.tbl
6
```

where the second line contains the SWAN output in table format and where the third and last line contains the number of time instants generated by the SWAN (we have 6 hourly instants in this case).

Then, open a dos terminal (dos shell) and type:

```
> ReadSwanRes
```

This will execute the program and generate the results.

The ReadSwanRes conversion will yield one hdf5 file per time instant (6 files, in the example above). Thus a final [glue](converttohdf5.md) action with the [ConvertToHdf5](converttohdf5.md) tool will allow to generate a final Waves.hdf5 to which configure MOHID [Module Waves](../water/modules/waves.md). Please refer to the [ConvertToHDF5](converttohdf5.md) wiki for samples on how to perform a [glue](converttohdf5.md) action.

## Configuring MOHID to run with the module waves

You can now go check this [wiki](../water/modules/waves.md).

## See also

- [Module Waves](../water/modules/waves.md)
- [SWAN Home page](http://citg.tudelft.nl/index.php?id=17787&L=1)
