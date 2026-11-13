# PROJ4

This page explains what is [PROJ4](http://proj.maptools.org/) library and how to use it from the command line, in Windows and in Linux.

[PROJ4](http://proj.maptools.org/) is an open-source Cartographic Projections library originally written by Gerald Evenden then of the USGS. This library has almost all projection types and it is very usefull to convert coordinates for a series of points.

It is written in C but a Fortran wrapper was developed as part of the [GLIMMER](http://www.cbcb.umd.edu/software/glimmer/) project. More information can be found in [libfproj4](libfproj4.md)

## Installation

#### Windows

- Download [http://download.osgeo.org/proj/proj446\_win32\_bin.zip](http://download.osgeo.org/proj/proj446_win32_bin.zip)
- Extract and copy "proj" folder to C:\\
- Right-click "My Computer" and select "Properties"
- Select "Advanced" tab
- Select "Environment Variables"
- In the "System Variables" box select "Path" and click "Edit".
- At the end of the line, add "C:\\proj\\bin" and click "OK".
- In the "System Variables" box, click "New".
- "Variable Name" is PROJ\_LIB and "Variable Value" is "C:\\proj\\nad".
- Click OK to close windows.

Now test it:

- Windows/Run/cmd
- cd .. until prompt is C:\\>
- write "proj -V" to see the version of the library.

#### Linux

- Go to your Package Manager or download source code from [http://trac.osgeo.org/proj/](http://trac.osgeo.org/proj/) and follow the instructions.

## Usage

The Proj4 is well documented. Here are some commands that I found useful.

#### Basic

##### degrees to meters (forward projection)

```
proj -f "%.2f" +proj=utm +zone=29 +ellps=WGS84
-7.1 38.2                                         (input)
666368.11       4229712.04                        (output)
```

##### meters to degrees (inverse projection)

```
proj -I -f "%.2f" +proj=utm +zone=29 +ellps=WGS84
666368.11 4229712.04                                 (input)
-7.10   38.20                                        (output)
```

#### Portuguese Military Coordinates (Hayford-Gauss Military System)

- Projection: Gauss-Kruger
- Ellipsoid: Hayford (internacional 1924): a=6378.388Km, b=6356.912 Km
- Datum : Lisbon (Castelo S.Jorge: 38º42’43.631’’N, 9º07’54.862’'W, 111.229m)
- Central Point: 39º40’N, 8º07’54.862’’W
- Origin: -200Km horizontally and -300 Km vertically

```
proj -f "%.2f" +init=epsg:20790
```

or

```
proj -f "%.2f" +init=esri:102164
```

or

```
proj -f "%.2f" +proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1.000000
               +x_0=200000 +y_0=300000 +ellps=intl +pm=lisbon +units=m +no_defs  no_defs
```

#### Universal Transverse Mercator (zone 29N)

```
proj -f "%.2f" +init=esri:32629
```

or

```
proj -f "%.2f" +proj=utm +zone=29 +ellps=WGS84 +datum=WGS84 +units=m
```

#### List of projections

List of projections in PROJ4: proj -lp

#### List of ellipsoids

List of ellipsoids in PROJ4: proj -le

More popular:

- -   International 1909 (Hayford) intl a=6378388.0m b=6356912m rf=297.
    -   New International 1967 new\_intl a=6378157.5m b=6356772m
    -   WGS 84 WGS84 a=6378137.0m b=6357002m rf=298.257223563
    -   Normal Sphere (r=6370997) sphere a=6370997.0m b=6370997m

b = a\*(1-f) = a\*(1-1/rf)

#### Portuguese Coordinate Systems in PROJ4

You might want to get the libraries containing the coordinate transformations with Portuguese datum pre-configured [here](http://www.fc.up.pt/pessoas/jagoncal/coordenadas/parte2.htm). Install the libs in the proj\_lib folder, inside the [FWTools](http://fwtools.maptools.org/) folder. By the way, you may want to [install FWTools first](http://fwtools.maptools.org/).

Here's the command line on how to convert from military portuguese to wgs84:

```
 cs2cs +init=pt:dlxhgmil +to +init=pt:etrs89 < coords_input.txt
```

## Internal Links

- [libfproj4](libfproj4.md)
- [ConvertToHDF5](converttohdf5.md)

## External Links

- [PROJ4](http://proj.maptools.org/)
- [Sistemas Referencia Portugueses](http://www.isa.utl.pt/dm/sigdr/sigdr04-05/9_SistRefPortugueses.ppt)
- [PROJ4\_HandlingCoordinateSystems.pdf](http://dl.maptools.org/dl/omsug/osgis2004/PROJ4_HandlingCoordinateSystems.pdf)
- [Projections and Projected Coordinate System Formulas](http://www.posc.org/Epicentre.2_2/DataModel/ExamplesofUsage/eu_cs34.html)
