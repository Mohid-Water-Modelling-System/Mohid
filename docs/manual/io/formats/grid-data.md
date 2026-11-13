# Grid Data

## Overview

Grid Data items relate to 2D or 3D information referenced to a grid. Grid data is stored in an ASCII text file and has multiple and flexible formatting options. In order to view a Grid Data file in [Mohid GIS](../../tools/gui/mohid-gis.md) the file must contain [grid](grid.md) information. In Mohid, the grid information is not used except for [bathymetry file](bathymetry.md) to [Mohid Water](../../water/mohid-water.md), or the input topography file to [Mohid Land](../../land/mohid-land.md) where the simulation domain is defined. Grid Data items can be used to initialize 2D or 3D fields (e.g. temperature, salinity) via [Module FillMatrix](../../base/module-fillmatrix.md).

## Header Section

The header section contains some specifications on the type of organization of the information.

```
FILL_VALUE         : -9.9e15
TYPE_ZUV           : z
```

FILL\_VALUE specifies a value for points with no data.

TYPE\_ZUV specifies if the data is located in the center of the grid cell or referenced to the faces of the grid cell. It can be **z** for centered data, **u** for east-west faces or **v** for north-south faces.

## Data Section

### Grid Data 2D

The data section contains the values of the Grid Data item. The data is stored inside a block defined by the following keywords: **<BeginGridData2D>** and **<EndGridData2D>**. Note that these keywords are case sensitive.

#### Filling data for all grid cells

You can specify the values by giving a list of all values of all the grid cells. Each value must be stored in a single line. The values are read from the lower left corner to the upper right corner of the grid that is values reading is done by row-column order. The first grid data values are read in the following order: (ILB, JLB), (ILB, JLB+1), … (ILB, JUB), (ILB+1, JLB), …, (IUB, JUB-1), (IUB, JUB). Note that the index “i” refers to the Y-axis and the “j” index to the X-Axis. In this case, if you have a (IUB x JUB) grid then the list must have (IUB x JUB) values.

```
<BeginGridData2D>
12.5
13.1
10.2
...
...
<EndGridData2D>
```

#### Filling data for some grid cells only

If you don't want to define all the values in the grid but only in a subset of cells in the grid (for example, if you want to simulate a bump in a flat bathymetry. Everywhere is 1000 m depth except in the bump area, so you want to define only the bump area); then you can also use ***"i j value(i,j,:)"***. Values not specified in the list will be given the DEFAULTVALUE.

```
<BeginGridData2D>
2   1   12.8
2   2   13.4
5   1   11.3
..  ..  ..
<EndGridData2D>
```

If the array is 3D and you wish to specify that in the vertical direction all the values are constant then use ***i j value(i,j,:)***. In the following example, if you have 3 layers, then grid cells (1,1,1), (1,1,2) and (1,1,3) will all be given the value of 10.2. Values not specified in the list will be given the DEFAULTVALUE.

```
<BeginGridData2D>
2   1   12.8
2   2   13.4
5   1   11.3
..  ..  ..
<EndGridData2D>
```

### Grid Data 3D

The data section contains the values of the Grid Data item. The data is stored inside a block defined by the following keywords: <BeginGridData3D> and <EndGridData3D>. Note that these keywords are case sensitive.

Between these keywords, for each grid cell, a value must be supplied (.e.g. for the [bathymetry file](bathymetry.md) write the bathymetry value).

If the array is 3D and you which to specify that in the vertical direction the values are constant then use *"i j value(i,j,:)"*. In the following example, if you have 3 layers, then grid cells (1,1,1), (1,1,2) and (1,1,3) will all be given the value of 10.2. Values not specified in the list will be given the DEFAULTVALUE.

```
<BeginGridData3D>
    1   1   10.2
    1   2   12.3
    2   5   10.7
    ..  ..  ..
<EndGridData3D>
```

If the array is 3D and you wish to specify to each grid cell a specific value the use *i j k value(i,j,k)*. Values not specified in the list will be given the DEFAULTVALUE.

```
<BeginGridData3D>
    1   1   1   10.2
    1   1   2   10.1
    1   1   3   9.9
    1   2   1   12.3
    ..  ..  ..  ..
<EndGridData3D>
```

You can also specify the values by giving a list of all values of all the grid cells. In this case, if you have a (IUB x JUB x KUB) grid then the list must have (IUB x JUB x KUB) values.

```
 <BeginGridData3D>
    10.2
    10.1
    9.9
    12.3
    12.1
    ...
 <EndGridData3D>
```
