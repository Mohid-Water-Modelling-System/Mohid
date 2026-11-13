# Geometry diagnostic

Once the vertical discretization or [geometry](../base/module-geometry.md) is imposed and the bathymetry is chosen, the bottom layer can yield stability problems when using shaved cells. You can have very thin bottom cell next to a very wide bottom cell. To diagnose the existence of such problematic cells, a **geometry diagnostic** tool was developed.

## Source-safe whereabouts

The latest source code is in the **Mohid\_V4** repository, inside the **small tools\\GeometryDiagnostic** folder.

## Configuration

The inputs required are a bathymetry [grid data](../io/formats/grid-data.md) file, a [geometry](../base/module-geometry.md) file, an hdf5 output filename and a [grid data](../io/formats/grid-data.md) output filename.

### Nomfich.dat

```
IN_MODEL    : GeometryDiagnostics.dat
```

### GeometryDiagnostics.dat

```
BATIM    : Biscay_Level2__.new2
GEOMETRY : D:\Aplica\BiscayAplica\GeneralData\Ocean\Geometry\Level2\Geometry.dat
OUT_FILE       : geomdiag_new2.hdf5
OUT_GRIDDATA3D : geomdiag_new2.dat
```

## Outputs

The outputs are an hdf5 file and a three-dimensional [grid data](../io/formats/grid-data.md) file.

### LowerLayerThickness

For all horizontal cells *i*, *j*, there is one bottom layer *kbottom* that defines the bottom layer thickness *DWZ(i, j, kbottom)*; where *DWZ* is the 3D matrix that defines the layers thickness for every water column.

```
LowerLayerThickness(i,j) = DWZ(i,j, kbottom)
```

### LowerLayerThicknessGradU

Is the ratio of neighboring cells bottom layer thickness along the X axis.

```
aux = DWZ(i,j, kbottom)/DWZ(i,j-1, kbottom)
LowerLayerThicknessGradU(i,j) = min( aux, 1/aux)
```

### LowerLayerThicknessGradV

Is the ratio of neighboring cells bottom layer thickness along the Y axis.

```
aux = DWZ(i,j, kbottom)/DWZ(i-1,j, kbottom)
LowerLayerThicknessGradV(i,j) = min( aux, 1/aux)
```

## See also

- [Module Geometry](../base/module-geometry.md)
