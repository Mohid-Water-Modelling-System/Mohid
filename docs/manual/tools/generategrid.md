# GenerateGrid

## Overview

GenerateGrid is a [MOHID support tool](mohid-support-tools.md) used to create regular grids with a non-constant grid size. It reads an input data file named *GridGenerator.dat*.

## Input data file (GridGenerator.dat)

Below is a sample input data file.

```
OUTPUT_FILE            : Sample.grd          !Output file name
ORIGIN_X               : -9                  !Coordinate X of the lower left corner
ORIGIN_Y               : 39                  !Coordinate Y of the lower left corner
GRID_ANGLE             : 0                   !Rotation angle of the grid
COORD_TIP              : 4                   !Coordinate type

<begin_grid_xx>
NUMBER_OF_SEGMENTS     : 1                   !Number of segments in the XX direction
NODES                  : 0.    2             !Coordinates starting in 0. of the XX nodes
RESOLUTION             : 0.1   0.1           !Resolution at each of the defined XX nodes
<end_grid_xx>

<begin_grid_yy>
NUMBER_OF_SEGMENTS     : 2                   !Number of segments in the YY direction
NODES                  : 0.     1.5    3     !Coordinates starting in 0. of YY nodes
RESOLUTION             : 0.1    0.03   0.15  !Resolution at each of the defined YY nodes
<end_grid_yy>
```

With these options a grid file is created, with origin in -9ºW and 39ºN (the approach is independent of the type of coordinates), with a constant spacing in the XX direction (resolution of 0.1º) and a non-constant spacing in the YY direction. In the YY direction, 2 segments are defined. The first starts at 0º from the origin and ends at 1.5º from the origin. The resolution of this segment starts at 0.1º and progressively decreases to 0.03º at the end of the segment. The second segment starts at the end of the first and ends at 3º from the origin. This means that the entire grid is 3º long in the YY direction. The resolution of the second segment starts at 0.03º and progressively decreases to 0.15º at the end of the segment.

## Fast start

How to create a grid. Go to MOHID Gis. Create NewData items choosing simple grid, then adjust coordinates.

## See Also

- [How to generate regular grid under MOHID GIS](../io/formats/grid.md)
