# Lines

## Overview

Lines items are defined by the starting point and ending point. Lines items can be stored as a collection of lines. This type of items can be used to define, for example, reaches or a coast line.

## Lines file

Lines items are stored in an ASCII text file. The default extension for Lines files in [MohidGIS](../../tools/gui/mohid-gis.md) is **\*.lin**. This file is organized in [blocks](../keywords/blocks.md), each containing the vertices for a determined polygon. One block is defined by the following [keywords](../keywords/keywords.md): **<begin\_line>** and **<end\_line>**. Note that these keywords are case sensitive. One file can contain an infinite number of blocks. Inside each block there are 2 columns and 2 rows: the first row is relative to the line starting point. The second row is relative to the line ending point. The left column relates to the XX coordinates of a point and the right one to the YY coordinates.

## Sample file

```
<begin_line>
-9.5190   38.8148
-9.4545   38.5699
-9.2612   38.4098
-9.0372   38.4686
-9.0550   38.7375
-9.1291   38.8921
-9.2644   38.8792
-9.4384   38.8438
-9.5254   38.8373
-9.5190   38.8148
<end_line>
```
