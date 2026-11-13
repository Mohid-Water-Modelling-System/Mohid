# XYZ points

## Overview

XYZ points are use do define discrete geo-referenced information. This information is stored associating X and Y coordinates to an attribute. Examples of applications of these items are: bathymetric or topographic surveys, water quality analysis results, etc.

## XYZ Points file

XYZ Points items are stored in an ASCII text file. The default extension for XYZ Points files in [MohidGIS](../../tools/gui/mohid-gis.md) is **\*.xyz**. This file is organized in [blocks](../keywords/blocks.md), each containing the coordinates for a certain attribute. One block is defined by the following keywords: **<begin\_xyz>** and **<end\_xyz>**. Note that these keywords are case sensitive. One file can contain an infinite number of blocks. Inside each block there at least 3 columns. The left column represents the XX coordinates of the polygon vertices and the right one to the YY coordinates of the polygon vertices. The third column contains the attribute value. A fourth column can be optionally included and it is handled as a character string which can also be used as a legend for each point.

## Example

### Three columns

```
<begin_xyz>
-9.138328  38.737320  292.
<end_xyz>
```

### Four columns

```
<begin_xyz>
-9.138328  38.737320  292.  MARETEC - IST (Lisbon)
<end_xyz>
```
