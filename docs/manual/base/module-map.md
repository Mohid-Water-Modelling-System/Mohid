# Module Map

## Overview

## Mapping variables

### WaterPoints3D

All points where WaterPoints2D = 1 AND k >= kfloorZ

### LandPoints3D

All points where LandPoints2D = 1 and WaterPoints2D = 1 and k < kfloorZ

### OpenPoints3D

All points where:

```
ComputeFaces3D%U(i, j, k) + ComputeFaces3D%U(i, j+1, k) +
ComputeFaces3D%V(i, j, k) + ComputeFaces3D%V(i+1, j, k) +
ComputeFaces3D%W(i, j, k) + ComputeFaces3D%W(i, j, k+1) > 0
```

### ComputeFaces3D

#### U

All faces where ComputeFace2D%U == 1 and k >= kfloorU

#### V

All faces where ComputeFace2D%V == 1 and k >= kfloorV

#### W

All faces in (k+1) where:

```
ComputeFaces3D%U(i, j, k) + ComputeFaces3D%U(i, j+1, k) +
ComputeFaces3D%V(i, j, k) + ComputeFaces3D%V(i+1, j, k) > 0
```

and expect where BoudaryPoints2D(i, j) == 1

### LandBoundaryFaces3D

All faces which have one side of the face a WaterPoints3D and on the other side a LandPoints3D

### ImposedNormalFaces

All faces which are used to define the boundary conditions normal to the boundary

### ImposedTangentialFaces

All faces which are used to define the boundary conditions tangential to the boundary
