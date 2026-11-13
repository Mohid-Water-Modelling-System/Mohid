# Filter.dat

****Filter.dat**** is the name of the file read by the [FilterBathymetry](../tools/filterbathymetry.md) tool to specify its options.

Input [Bathymetry](formats/bathymetry.md) in Mohid [Grid Data](formats/grid-data.md) format

```
    IN_BATIM        : InputBathymetry.dat
```

Output Bathymetry in Mohid GridData format

```
    OUT_BATIM       : OutputBathymetry.dat
```

Check if the user wants to do an automaitc smoothing in a way that the slope of each cell is below a slope parameter

```
    SIGMA_SMOOTH    : 1
```

Slope parameter in the bibliography this value for sigma models should be lower than 0.2

```
    SLOPE_PARAMETER : 0.2
```

Example for filtering the intertidal area

Minimum value to be filter

```
    H_MIN           : -20
```

Maximum value to be filter

```
    H_MAX           : 0
```
