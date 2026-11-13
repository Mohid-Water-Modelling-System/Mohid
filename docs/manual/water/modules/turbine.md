# Module Turbine

## Overview

Module Turbine is in charge of calculating the force exerted by turbines in the water flow as well as the power and energy extraction. This implementation is available for 2D and 3D simulations being, in the 3D model, only discretised in the vertical domain. There is no discretisation of the turbine geometry (and force) in the horizontal domain, the model computes the force as a punctual force in the cell where the turbine is placed.

For the moment, the model is thought for horizontal bidirectional turbines, with pitch control and free rotation in the vertical axis (the perpendicularity between flow and turbine is implicit in the model as the force is introduced into the hydrodynamic model as a punctual force).

## Implementation

### Forces

**\\text{F}\_{T} = \\frac{1}{2}\\rho A\_{T} C\_{T} U^{2}**

**\\text{P}\_{T} = \\frac{1}{2}\\rho A\_{T} C\_{P} U^{3}**

**Thrust force produced by the turbine rotor due to the energy extraction**

*F\_{T} = (1/2)ρA\_{T}C\_{T}U^{2}*

**The power extracting by the turbine**

*P\_{T}=(1/2)ρA\_{T}C\_{P}U^{3}*

Where *ρ* is the water density, *A\_{T}* is the area swept by the blades, *C\_{T}* is the thrust coefficient that quantifies the force exerted by the turbine to the flow and *C\_{P}* is the power coefficient that quantifies the amount of power extracted from the flow.

### Coefficient parameterisation

*C\_{T}=*

- *0* if U≤U\_{C}
- *C\_{T\_0}* if U\_{C}<U≤U\_{D}
- *C\_{T\_0}U\_{D}^{3}/U^{3}* if U>U\_{D}

*C\_{P}=*

- *0* if U≤U\_{C}
- *C\_{P\_0}* if U\_{C}<U≤U\_{D}
- *C\_{P\_0}U\_{D}^{3}/U^{3}* if U>U\_{D}

Where C\_{T\_0} and C\_{P\_0} are the design values for each coefficient and U\_{C} and U\_{D} are the cut-in and design speed, respectively

### Vertical discretisation

The vertical discretization is made through the equation

- *A\_{T\_k}=r^{2}θ/2-rsin(θ/2)d-∑\_{k=1}A\_{T\_k-1}*

Where *A\_{T\_k}* is the area swept by the blades in layer *k*, *r* is the blade length, *θ* is the angle between two consecutive blades and *d* is the vertical distance between the centre of the blade and the end of the blade.

The force excreted by the turbine is calculated with the velocity of the flow in each layer while for the power, the velocity is an average value between the layers that contain the turbine. So, the fores can be rewritten as follow:

- ***F\_{T\_k}**\=(1/2)ρA\_{T\_k}C\_{T}U\_{k}**U\_{k}***

- *P\_{T}=(1/2)ρA\_{T}C\_{P}U\_{AV}^3*

Where ***U\_{k}*** is the velocity vector of the turbine in an Arakawa C manner, *U\_{k}* is the velocity modulus and *U\_{AV}* is the average modulus velocity of the k layers of the cells in the coordinates i,j that contain the turbine calculats as:

*U\_{AV}=∑\_{k}A\_{T\_k}U\_{k}/∑\_{k}A\_{T\_k}*

## Limitations of the model

The origin of the main limitations of this implementation is the non-discretisation in the horizontal plane.

The implementation is not programed to work with really small grid size, the results will be always inappropriate when the cell size is under the turbine diameter thought there is no discretisation. Also it is not prepared to work with really large grid size, because visual data doesn’t show significant changes due to the resolution and also the calculation are not really accurate.

Another limitation derived of the horizontal discretisation is the layout array studies. Actually there are some limitations with the distances and positions between the turbines because the minimal distance between them will be, at least, as big as the grid cell size diameter.

## Input data format

## Configuration files

In order to activate the implementation, some modifications need to be done to the main configuration files. Two files need to be changed, the hydrodynamic.dat file and the nomfich.dat. In the hydrodynamic file the keyword TURBINE needs to be written with a value of 1 for activating the implementation. If it is not written or with a value of 0 it will not work.

image

The file Nomfich.dat is where the routes of the files needed for the simulations are specified. So, in order to make the implementation work and that the model can read the input data of the turbines, this route should be written in this file. The keyword for this file is TURBINE.

image

## Thesis link
