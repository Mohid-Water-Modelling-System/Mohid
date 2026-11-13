# Open Boundary conditions

## WaterProperties\_x.dat

For each property it's possible to define a specific boundary condition. The available options are:

1 - 1D horizontal mass balance normal to the boundary

2 - An imposed value

3 - 1D vertical diffusion is added to option 1

4 - Null gradient

5 - SubModel, a value interpolated directly from the father is imposed

6 - Orlanski type boundary condition

8 - Cyclic boundary condition

By default, option 1 is assumed.

```
    <beginproperty>
     BOUNDARY_CONDITION      : 1
    <endproperty>
```

By default the value assumed for the boundary is the one assumed in the intialization process of the property field (INTERIOR) the user can also specify a value for the boundary (EXTERIOR) to be set in the DEFAULTBOUNDARY keyword.

```
    <beginproperty>
     BOUNDARY_INITIALIZATION     : EXTERIOR
     DEFAULTBOUNDARY             : 0.
    <endproperty>
```
