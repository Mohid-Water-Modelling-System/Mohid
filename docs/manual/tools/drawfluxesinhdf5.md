# DrawFluxesInHDF5

**DrawFluxesInHDF5** is a little tool that reads the boxes fluxes time series (.BXF and .BXM) and returns an HDF5 file with vector fields representing the fluxes between the boxes.

## Source-safe whereabouts

The latest source code is in the [codeplex](../archive/web/codeplex.md) inside the "SmallTools\\DrawFluxesInHDF5" folder. An older source code is in the **Mohid\_V4** repository, inside the **SmallTools\\DrawFluxesInHDF5** folder.

## Configuration

### Nomfich.dat

The nomfich file contains a single line pointing to the configuration file.

```
IN_MODEL : CreateHDF5Boxes.dat
```

### CreateHDF5Boxes.dat

Here's a sample configuration file (CreateHDF5Boxes.dat):

```
HDF5BOXES_INFILE : Boxes.hdf5
HDF5BOXES_OUFILE : OutBoxes.hdf5
BOXESSALINITYFILE : salinity.bxm
BOXESWATERFILE : water.bxm
BOXESWATERFLUXFILE : water.bxf
INDEXESFILE : BoxesFluxesIndexes.txt
INDEXESLENGTH : 436
```

1. The first line contains the input boxes.hdf5 file. They represent the index map of the boxes created with the [Box2HDF5](box2hdf5.md) tool.
2. The second line tells the output filename.
3. The third line contains the salinity inside the boxes.
4. The fourth line contains the boxes water volumes timeseries.
5. The fifth line contains the fluxes between boxes timeseries.
6. The sixth line contains three columns displaying the column number, the source box index and the destination box index. Each row represents a flux of the water fluxes timeseries.
7. The seventh and last line contains the total number of fluxes available in the water fluxes timeseries.

## BoxesFluxesIndexes.txt

Here's a sample of the indexes file:

```
9 1 2
10 2 3
11 5 6
12 6 7
13 6 38
14 7 1
15 7 8
16 7 39
17 8 2
18 8 9
19 8 40
20 9 3
```

1. First column: column index number.
2. Second column: source box index.
3. Third column: target box index.

## Example output

## Troubleshoot

**Q: Too many columns?**

A: Add the following preprocessor directives in Module [Mohid Base 1](../water/README.md#mohid-base-1#mohid-base-1):

```
_BIG_MAX_COLUMNS,_EXTRA_LONG_LINE_LENGTH
```

## See also

- [Boxes](../io/formats/boxes.md)
- [Box2HDF5](box2hdf5.md)
