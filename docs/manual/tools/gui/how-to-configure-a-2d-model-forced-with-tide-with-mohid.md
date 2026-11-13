# How to configure a 2D model forced with tide with MOHID

Before you start creating your own simulation, you may want to create a [bathymetry](how-to-create-a-model-gridded-bathymetry.md) and [tidal gauges](how-to-generate-tide-for-mohid.md) first.

## Step 1 - Create a new mohid GUI project

- Create a new [Mohid GUI](mohid-gui.md) project.

- Name the project and select a location

## Step 2 - Insert a new simulation and a new run

- Right-click and select Insert Simulation

- Name the simulation, browse and select a bathymetry file, browse and select a tidal gauge data

- Right-click on the simulation and select *New Run*

- Edit at will the start date, the end date and the simulation time-step, then click *Ok*

In order to choose an adequate time-step, $\Delta t$, we recommend to aim to a 1D [Courant number](http://en.wikipedia.org/wiki/Courant_number) of 15; $Cr = c \frac{\Delta t}{\Delta x}$. In estuarine and shelf waters, you may estimate the system maximum velocity, 

$$
c
$$

, with $c = \sqrt{g\,H} </mathtex>
where <math>H$ is the maximum depth of your domain and $g$ is local gravitational field (approximately worth 10 m/s2). For example, with a 200 m resolution domain where maximum depth is 400 m, the approximate time step would be around 50 s.

## Step 3 - Configure the modules

To configure the modules simply double-click on the module icon in the GUI.

Warning: the configuration files are made of a list pairs of keywords and values. **Only whitespaces are admitted.** Any tabulation will yield in error during the reading of the keywords.

### [Module Model](../../water/modules/README.md#module-model#module-model)

The module model determines the start and end times, as well as the model's time-step. It is probably already configured and should look like this:

```
START                         : 2008 10 28 12 0 0
END                           : 2008 10 31 12 0 0
DT                            : 50.
VARIABLEDT                    : 0
GMTREFERENCE                  : 0
SPLITTING                     : Double_Splitting
```

### [Module Atmosphere](../../water/modules/atmosphere.md)

The module atmosphere configures the atmospheric forcing such as wind, solar radiation, latent and sensible heat sea-air exchange, infrared radiation. In this test-case, we have no atmospheric forcing, thus the file should be left in blank.

### [Module Geometry](../../base/module-geometry.md)

The module geometry describes the vertical discretization of the domain. In this test-case we have a single vertical sigma layer, since the model is 2D.

```
MINIMUMDEPTH            : 0.1000

<begindomain>
ID                      : 1
TYPE                    : SIGMA
LAYERS                  : 1
LAYERTHICKNESS          : 1.0000
DOMAINDEPTH             : -99.00
TOLERANCEDEPTH          : 0.0500
<enddomain>
```

### [Module InterfaceWaterAir](../../water/modules/interfacewaterair.md)

This module parameterizes the exchanges of momentum, radiation, heat and salt between the atmosphere and the water. In this case, we simply need to define the surface rugosity.

```
<begin_rugosity>
INITIALIZATION_METHOD        : CONSTANT
DEFAULTVALUE                 : 0.0025
<end_rugosity>
```

### [Module InterfaceSedimentWater](../../water/modules/interfacesedimentwater.md)

This module parameterizes the exchanges of momentum, mass and heat between the sediments and the water. In this case, we simply need to define the bottom rugosity.

```
<begin_rugosity>
INITIALIZATION_METHOD        : CONSTANT
DEFAULTVALUE                 : 0.0025
<end_rugosity>
```

### [Module Hydrodynamic](../../water/modules/hydrodynamic.md)

This module parameterizes the numerical model solving the advection-diffusion equations of momentum combined with the continuity equation. It includes also boundary conditions definitions, as well as initialization conditions. Finally, assimilation schemes of U, V or water level. The momentum equations are also modularized, in the sense that each of its components may, or may not, be deactivated.

In this case we are simply choosing a TVD spatial numerical scheme, while adding the Coriolis force to the momentum equations. Since the model is 2D and has homogeneous density, the baroclinic term in the pressure gradient is deactivated, for performance sake. The model starts at rest and has the water level forced with tide at the open-boundaries. The model outputs results every 3 hours. Furthermore it will also output residual (averaged) hydrodynamic results and a domain-integrated energy time-serie.

```
ADV_METHOD_H            : 4
ADV_METHOD_V            : 4
TVD_LIMIT_H             : 4
TVD_LIMIT_V             : 4

BAROCLINIC              : 0

CORIOLIS                : 1

TIDE                    : 1

OUTPUT_TIME             : 0. 10800.

RESIDUAL                : 1
ENERGY                  : 1
```

### [Module Turbulence](../../water/modules/turbulence.md)

The module Turbulence configures how the horizontal and vertical sub-grid scale turbulent viscosities are calculated. In this case we opt for constant sub-grid scale viscosities.

```
!Unit: m2/s
VISCOSITY_V             : 0.0010
VISCOSITY_H             : 20
```

### [Module WaterProperties](../../water/modules/waterproperties.md)

The module Waterproperties applies the advection-diffusion equation to any eulerian tracer. It could be the salinity, the temperature, or any other property that is advected and diffused in the water. In this case we have a homogeneous density, thus we don't need to calculate the density field, thus there's no need to transport salinity or temperature or any other waterproperty for that matter. Hence we leave a blank file in the module Waterproperties in this test-case.

## Step 4 - Running the model

If all is correctly configured, and all the input files are double-checked, we can hopefully start running our model!

- Select *Tools* from menu, then *Launch Mohid...*

- Select the *Run\_1*, then select *Run Mohid*, then select *Launch after Ok*, then click on *Ok*

- Supervise the screen-outputs to check for any error, or to simply read the estimated end time of simulation

## Step 5 - Inspect the results

Once the model has finished running, (or while the model is running, after it has created a few HDF5 outputs), you will want to switch in analysis mode and graphically inspect the results to see if all went well (or is going well). We will proceed to visualize the instantaneous velocity as a vector field, and the instantaneous water level as a colormap field. The proceedings are exactly the same for any other variables.

- Click on HDF file, then click on Grid HDF group

- Map the GridY, GridX, VerticalZ and Mapping with the Group items Latitude, Longitude, VerticalZ and OpenPoints, respectively.

- Map the time HDF Group to the Time selection button

- Map the Velocity U, Velocity V and Water level HDF results group's items to the Vector X, Vector Y and Color selection buttons, respectively

- From the *Action menu*, select *Run OpenGL*

- Select *settings*

- Select *Color*

- Define the minimum water level (m) and the maximum water level (m). Then click on *Ok*.

- Select *Vectors*

- Define the vector *Scale* to 0.2 and the vector stepping along X and along Y to 6. Then click on *Ok*.

- Render the results, frame by frame

If all went well, you should see results with patterns and values as expected for the domain. Now you can try to add time series, or to improve the open boundary conditions ...

## Next steps

- [How to create MOHID timeseries outputs](how-to-create-mohid-timeseries-outputs.md)
- [How to create discharges in MOHID](how-to-create-discharges-in-mohid.md)
- How to use the Blumberg-Kantha radiation method
- [How to couple atmosphere forcing](../../water/coupling-water-atmosphere-user-manual.md)

## Previous steps

- [How to create a model gridded bathymetry?](how-to-create-a-model-gridded-bathymetry.md)
- [How to generate tide for Mohid?](how-to-generate-tide-for-mohid.md)

## About the screen capturing tool

The screen capturing tool used for this tutorial was the [Jing Project](http://www.jingproject.com/).
