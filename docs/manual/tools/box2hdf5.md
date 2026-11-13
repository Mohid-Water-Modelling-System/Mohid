# Box2HDF5

Box2HDF5 is a little tool that creates a HDF5 file from a [Boxes](../io/formats/boxes.md) definition file, a geometry file and a bathymetry file. This hdf5 file serve for checking that the desired boxes definition is correct in a visual way and also serves as input file for the tool [DrawFluxesInHDF5](drawfluxesinhdf5.md) that returns an HDF5 file with vector fields representing the fluxes between the boxes.

## Source-safe whereabouts

At the moment is not found neither in the sourcesafe nor in the Codeplex

## Configuration

### Nomfich.dat

The nomfich file contains three lines pointing to the bathymetry, geometry and boxes files.

```
IN_BATIM         : TagusRibeiras.dat
IN_GEOM          : Geometry_1.dat
IN_BOXES         : WQBoxesTagusRibeira3D.dat
```

## See also

- [Boxes](../io/formats/boxes.md)
- [DrawFluxesInHDF5](drawfluxesinhdf5.md)
