# Profile

Profile file format is used to define a vertical profile of a certain property (e.g. temperature, salinity, phytoplankton). This format is useful to initialize 3D fields or to impose vertical gradients as a boundary condition.

Profile is an ASCII text file in which data is defined inside a [block](../keywords/blocks.md) structure **<BeginProfile>** and **<EndProfile>**. The block has 2 columns: the left column corresponds to the depth and the right column defines the property values at the correspondent depth.

One profile file corresponds to one set of values, i.e., one file contains only a set of data correspondent to a single property. This means, that one cannot define temperature and salinity in the same file.

## Sample file

```
<BeginProfile>
5 18.246
15 17.93
25 17.361
40 16.077
62.5 14.405
87.5 13.34
112.5 12.885
137.5 12.668
175 12.443
225 12.144
275 11.872
350 11.546
450 11.175
550 10.964
650 10.976
750 11.085
850 11.008
950 10.752
1050 10.449
1150 10.054
1250 9.281
1350 8.211
1450 7.166
1625 5.899
1875 4.44
2250 3.973
2750 3.0
3250 2.7
<EndProfile>
```
