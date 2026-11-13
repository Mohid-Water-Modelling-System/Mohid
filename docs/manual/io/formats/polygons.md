# Polygons

## Overview

Polygon items relate to areas defined by a collection of vertices. Polygon items can be stored as a collection of polygons. This type of items can be used to define types of zones such as differentiate land and sea, soil use, etc.

Polygons can be manually drawn on top of a Background Image or a Grid Data item, which usually represent geographical information. This allows visual distinguishing of areas in loaded maps.

## Polygons file

Polygons items are stored in an ASCII text file. The default extension for Polygons files in [MohidGIS](../../tools/gui/mohid-gis.md) is **\*.xy**. This file is organized in [blocks](../keywords/blocks.md), each containing the vertices for a certain polygon. One block is defined by the following keywords: **<beginpolygon>** and **<endpolygon>**. Note that these [keywords](../keywords/keywords.md) are case sensitive. One file can contain an infinite number of blocks. Inside each block there are 2 columns: the left relates to the XX coordinates of the polygon vertices and the right one to the YY coordinates of the polygon vertices.

A third column can also be added in order to aid in the construction of curvilinear grids. In this column a value of -1, 0 or 1 can be defined indicating the direction into which the polygon corner turns to (clockwise direction). Learn more in [Generating curvilinear grids](../../tools/generating-curvilinear-grids.md).

## Sample file

```
<beginpolygon>
-9.5190   38.8148   0
-9.4545   38.5699   0
-9.2612   38.4098   0
-9.0372   38.4686   0
-9.0550   38.7375   0
-9.1291   38.8921   0
-9.2644   38.8792   0
-9.4384   38.8438   0
-9.5254   38.8373   0
-9.5190   38.8148   0
<endpolygon>
```
