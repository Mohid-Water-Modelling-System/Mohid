# Digital Terrain Creator

The **Digital Terrain Creator** is a small executable used to generate bathymetries from baseline data. Its inputs are a [grid file(.grd)](../io/formats/grid.md) and an ascii [xyz points file(.xyz)](../io/formats/xyz-points.md)(or several [bathymetry](../io/formats/bathymetry.md) files for overlapping).

Optionally a land [polygons file(.xy)](../io/formats/polygons.md) can be added in order to force land points. Its output is the [bathymetry](../io/formats/bathymetry.md) [griddata file(.dat)](../io/formats/grid-data.md). The **Digital Terrain Creator** tool reads its option from the CreateBathymetry.dat file.

A new option was added to dredge sections (river section in land or channels in water). This is usefull in MOHID Land projects for river sections with fine resolution DTM's that have bottom erroneusly defined in water bodies or in MOHID Water to dredge specific channels. This options is for sections that are larger than one cell because the definition of section does not make sense inside one cell.

## CreateBathymetry.dat file

Here's what a CreateBathymetry.dat option file looks like:

```
!This is the output bathymetry file:
BATIM_FILE                  : GridData_2.dat

!This is the grid where the bathymetry will be created
GRID_FILE                   : ..\..\CostaPortuguesa12W36N_6p15W45NETV4.grd

!This is an initial bathymetry griddata file (optional)
!BATIM_INI                   : ..\..\CostaPortuguesa12W36N_6p15W45NETV3.dat
!MATCH_IN                    : 0

!Here lie the land polygons files
<BeginLandAreaFiles>
..\..\PortugueseLineCoast_new.xy
<EndLandAreaFiles>

!Here lie the baseline XYZ Points data files
<BeginXYZPointsFiles>
..\..\xyz\ET13W35N_5W46N.xyz
<EndXYZPointsFiles>

!Perform bathymetry smoothing within a specific radius
SMOOTH                        : 0 (1/0)
RADIUS                        : 0.05(dimensions of the grid file, for eg., degrees)

!Here is defined the interpolation method
INTERPOLATION               : Triangulation   !No Interpolation; Triangulation; Spline

!Here are the interpolation method options
POINTS_FOR_INTERPOLATION  : 3         !Method to select base points for interpolation 1-AllPoints; 2-GridCellCenterPoints; 3-PointsInsideGridLimits
EXPAND_GRID_LIMITS        : 1
GRID_LIMITS_PERCENTAGE    : 0.25
```

A new option was added to avoid processing points with Z with a given value (e.g. points that have no data). To connect this use the keyword with the Z value that should not be processed:

```
INVALID_Z_VALUE             : -9999
```

## Overlapping

To overlap several bathymetries, each one with a specific weight on the final bathymetry, add the following section to the options file. Overlapping requires an existent initial bathymetry and at least one bathymetry to overlap with.

```
BATIM_INI                   : Initial.dat    !the initial bathymetry
GRID_FILE                   : Initial.dat    !the grid (must be the same for the initial and others
BATIM_FILE                  : FinalFile.dat
OVERLAPPING                 : 1              !do overlapping
OVERLAPPING_NUMBER          : 1              !number of overlapping bathymetries

(Block for each overlapping bathymetries)
<BeginGridDataInfo>
LEVEL                     : 1
GRIDATA_FILE              : Overlap.dat      !the overlapping bathymetry
PERCENTAGE                : 1                !Fraction (0-1) of importance given to this Grid Data
	<<BeginAreaFiles>>
         overlap_area.xy                     !Path to the polygons files defining area of the
                                             !grid data file to be used. Each line defines one file
	<<EndAreaFiles>>
<EndGridDataInfo>
```

## River/Channel Dredging

Digital Terrain creator can also dredge a part of the DTM creating sections and using triangulation to fill between sections. For this the user provides the bathymetry, the area to dredge, the points collection with section begin and end and the program does the rest (creating points to define the sections and interpolating in between sections).

Also a new feature was added to be able to create more margin points and the user only needs to give the main river changes.

### Requirements

To dredge a river the algorithm needs:

- an original bathymetry (to be dredged).
- an area of original bathymetry to apply the dredging - this can be obtained in GUI by creating a polygon.
- points collection where each pair defines section begin and end (z is not read if SECTION\_MIDDLE\_POINT\_METHOD : 1 or is read and defines middle point Z if SECTION\_MIDDLE\_POINT\_METHOD : 2) - this collection can be obtained in GUI by creating points with any Z (it will be computed).

```
Remember: the first two consecutive points need to belong to one section (one for section start and the other for end) and the next two
points for the next section and so on.
Remember: the points defining section start and end must be outside the polygon area.
Remember: if using the option to create more margin points from original points then all the first points of each pair needs to belong
to one margin and the second points of each pair need to belong to the other.
```

- middle point height (or depth) in relation to section start point - read if SECTION\_MIDDLE\_POINT\_METHOD : 1.
- section shape (if linear the section will be "V" shape, if root will be "U" shape and if power will be "hearts" with flat plain and fast decrease to section center).

### Assumptions

The method assumes that the lower point of the section is exactly at middle between start and end and uses a unique value for all sections (the height of the middle point in relation to first point or section start) if SECTION\_MIDDLE\_POINT\_METHOD : 1. If the user wants to give middle point Z (if SECTION\_MIDDLE\_POINT\_METHOD : 2) define it in GUI when creating the file (define Z of first point of each pair) or after creation, editing xyz file.

### Process

The algorithm starts by finding for each XYZ point (sections start and end) the bathymetry level by searching in the original bathymetry file.

In the next step, the points needed to define each section are added and elevation computed (using section start and end elevation, middle point depth (or level) and section shape). In this step a xyz file with the original points (section begin and end) and the added points is outputed.

In the next step, if interpolation is connected, the defined area for dredging is prepared.

In the next step, if interpolation is connected, the interpolation process (triangulation) will use all the points to create the river/channel in between sections in the defined area. The bathymetry cells outside the defined area will remain as in the original bathymetry.

### Data File Example

```
BATIM_INI                   : Initial.dat    !the initial bathymetry
GRID_FILE                   : Initial.dat    !the grid
BATIM_FILE                  : FinalFile.dat

!Here lie the XYZ Points data files (points need to be even because 1st is section 1 start point, 2nd is section 1 end point
!3rd is section 2 start point, 4th is section 2 end point and so on. z value is undefined (can have any value) and wll be computed by the algorithm.
<BeginXYZPointsFiles>
SectionBoundaries.xyz
<EndXYZPointsFiles>

SECTION_COMPUTE             : 1               !1 - compute sections from given xyz points (pairs of section begin and end)
SECTION_MIDDLE_POINT_METHOD : 1               !1 - use middle height; 2 - use Z in original xyz collection
SECTION_HEIGHT              : 3.0             !section middle height or depth (in reference to first section point z) - read if SECTION_MIDDLE_POINT_METHOD : 1
SECTION_INTERPOL_TYPE       : 3               !1 - linear; 2 - power; 3 - root
SECTION_INTERPOL_ORDER      : 2               !power or root order (2 - power of 2 or squared root)
SECTION_POINTS_TO_ADD       : 7               !points to add to each section (beside section end and start)
!output file with original points (section boundaries) and created points
SECTION_XYZ_FILE            : FinalSections.xyz

CREATE_MARGIN_POINTS   : 1                    !1 - compute more margin points from original xyz points (interpolating)
MARGIN_POINTS_TO_ADD   : 7                    !points to add between each original points
```

!areas to select where to triangulate (if triangulation active)

```
<BeginAreaFiles>
SelectedArea.xy
<EndAreaFiles>

!Here is defined the interpolation method: No Interpolation; Triangulation; Spline
!this is optional, if disconnected only the river section points are created, if
!connected than the section points will be used to create an interpolated river between sections
INTERPOLATION               : Triangulation
FILL_OUTSIDE_POINTS         : 1
WRITE_TRIANGLES             : 1
TRIANGLES_FILE              : Triangles.xy

!Here are the interpolation method options
POINTS_FOR_INTERPOLATION  : 4
EXPAND_GRID_LIMITS        : 1
GRID_LIMITS_PERCENTAGE    : 0.25
```

## Interpolation methods

### No Interpolation

In the case of no interpolation only the cells that have some points inside will have a value (average of the points), the other will be left blank.

### Triangulation

The cells that are left blank can be interpolated by triangulation that is defined in ModuleTriangulation. The triangulation process program used was the developed by Robert Renka, Department of Computer Science, University of North Texas.
