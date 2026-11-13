# Module HorizontalMap

## Overview

## Mapping variables

### WaterPoints2D

All points where the Bathymetry > -55.

### OpenPoints2D

All points where ComputeFaces2D%U(i, j, k) + ComputeFaces2D%U(i, j+1, k) + ComputeFaces2D%V(i, j, k) + ComputeFaces2D%V(i+1, j, k) > 0

### LandPoints2D

All points where the Bathymetry < -90.

### BoundaryPoints2D

All points which are WaterPoints2D at the limit of the domain and WaterPoints2D close to Bathymetry points with a value of -80.

### ExteriorPoints2D

All points where:

- \-90 < Bathymetry < -55
- I = Size%ILB or I = Size%IUB
- J = Size%JLB or J = Size%JUB

### ComputeFaces2D

All faces which, which simultaneously:

- on both sides have Waterpoints
- at least one of the Waterpoints have a waterlevel above H\_Min
- are NO face between boundary points

### ExteriorBoundaryFaces

All faces which Exterior Faces of open boundary points and faces between boundary points and points of the bathymetry with the value of -80.

### BoundaryFaces

All faces which have in one side a interior point and in another a boundary point

### IMin(:), IMax(:)

Using these arrays it is possible to make:

- DO I = IMin(J), IMax(J)
- IMin is the first grid point with water
- IMax is the last grid point with water
