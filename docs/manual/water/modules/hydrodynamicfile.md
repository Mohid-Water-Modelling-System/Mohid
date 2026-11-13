# Module HydrodynamicFile

## Overview

Module HydrodynamicFile can be seen as an auxiliary module, which allows the user to integrate the hydrodynamic solution in space and in time and store this solution in a file. This file can later be used to simulate longer periods. The usage of the hydrodynamic file module has shown that the errors introduced in the integrated solution are usually small, and not significant for long-term water quality simulations.

## Methodology

### Spatial integration

The spacial integration consists in the integration of several grid cells into one single cell. This grouping can be done for any quadratic group of grid cells, like 2x2 or 3x3. This grouping results in a drastic reduction of computing points. In the case of space integration 2x2, the resulting domain will just contain one forth of grid points, each with the double size. For a space integration of 2x2, the number of the overall grid points reduces from 16 to 4.

The space integration is divided into to steps:

- Compute an integrated bathymetry, base on the bathymetry with a higher resolution
- In each integration step, the sum of all water fluxes along a given face is kept in the output file. The average surface elevation is also kept in the output file.

### Bathymetry integration

The integration of the bathymetry can be done in two different ways. In both ways, before the bathymetry is integrated, the land points with are filled with the minimum depth of the cells which are to be integrated. The first way of the integration is designated as "Mean Integration". This methodology calculates the average depth of the cells to be integrated, using this average as the depth for the new bathymetry. The second way of integration is designated as "Maximum Integration". This methodology uses the maximum depth of all cells, multiplied by the area of each cell and then divided by the total area of the water points.

### Integration of the water fluxes

The water flux between two cells is calculated from the average water flux between these cells. These fluxes can be integrated over several time steps.

### Time integration

The time integration consists in the integration along several time steps of the hydrodynamic solution. The time integration can be directly connected to the space integration, once the larger grid spacing obtained by the space integration, allows the model to run with a larger time step. In the case of mass transport, the celerity which controls the stability of the model isn’t the propagation of the pressure wave, but the maximum flow velocity. In the resulting file, the hydrodynamic solution of the *"in time"* and *"in space"* integrated hydrodynamic solution is stored, which in posterior simulation can be used to obtain results more quickly. Regarding time, there are two different ways of storing this information: the information can be stored as a integrated solution with a given start date and a given end date, or as a solution which repeats itself (one tidal cycle).

## References
