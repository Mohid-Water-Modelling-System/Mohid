# Module Sediment

## Overview

Module Sediment allows simulating the dynamics of sediments in natural systems under the action of waves and currents. It allows to assess the transport of multiple sediment fractions (non-cohesive and cohesive) and associated morphological evolution. The effects of mixed size fractions on sediment motion are taken into account.The bottom column can be divided into several layers to represent bottom stratigraphy. The consolidation of cohesive sediments can increase the bed density, affecting the resistance to erosion, the thickness of the sediment column and the bathymetry. As cohesive sediment form flocs, the diameter of individual particles becomes irrelevant for their transport. Therefore, cohesive sediment was assigned as a unique class, considering flocculation as a function of suspended sediment concentration and turbulence. In contrast, the number of sand classes is unlimited.

## Main options

In order to consider Module Sediment in [MOHID Water](../mohid-water.md) simulations, the following keyword must be activated in the [Module Model](README.md#module-model#module-model) input data file:

```
SEDIMENT : 1
```

The options of Module Sediment are defined in the input file ([Sediment.dat](../../io/sediment-dat.md)).

The geometry of the sediment column is given in a different file ([SedimentGeometry.dat](../../io/README.md#sedimentgeometry-dat#sedimentgeometry-dat)).

Other modules must be modified: ([InterfaceSedimentWater.dat](../../io/interfacesedimentwater-dat.md)) ([WaterProperties.dat](../../io/waterproperties-dat.md)) ([Free Vertical Movement.dat](../../io/free-vertical-movement-dat.md))

## References

- Franz, G., Leitão, P., Pinto, L., Jauch, E., Fernandes, L. and Neves, R., 2017. Development and validation of a morphological model for multiple sediment classes. International Journal of Sediment Research, 32(4), pp.585-596.
- Franz, G., Delpey, M.T., Brito, D., Pinto, L., Leitão, P. and Neves, R., 2017. Modelling of sediment transport and morphological evolution under the combined action of waves and currents. Ocean Science, 13(5), p.673.
