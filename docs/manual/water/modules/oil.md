# Module Oil

## Overview

The prediction and simulation of the trajectory and weathering of oil spills are essential to the development of pollution response and contingency plans, as well as to the evaluation of environmental impact assessments. In order to predict the behaviour of the oil products spilled in coastal zones, an oil weathering model was developed, which predicts the evolution and behaviour of the processes (transport, spreading and behaviour) and properties of the oil product spilled in water. Some pollution response methods are also integrated in the model.

## Implementation

Oil density and viscosity, and many different processes are included in oil module, such as oil spreading, evaporation, dispersion, sedimentation, dissolution, emulsification, oil beaching and removal techniques. Different alternative methods were coded for the prediction of some processes like oil spreading, evaporation, dispersion, sedimentation and emulsification. Therefore, when using the model, there is more than one way of simulating the same process, depending, for example, on the characteristics of the computational mesh or on the magnitude of the spill. The oil weathering module (OWM) uses mainly the 3D hydrodynamics and 3D lagrangian transport modules. The hydrodynamic module simulates the velocity field necessary for the lagrangian module to calculate oil trajectories. These oil trajectories are computed assuming that oil can be idealized as a large number of particles that independently move in water. Water properties and atmospheric conditions are introduced in lagrangian module and used by oil module for determination of oil processes and properties. Excepting spreading and oil-beaching, all weathering processes and properties are assumed uniform for all tracers, like water properties and atmospheric conditions, which are considered equal to these environmental conditions determined in accident origin. As it was already mentioned, the movement of the oil tracers can be influenced by the velocity field from the hydrodynamic module, by the wind from the surface module, by the spreading velocity from oil module and by random velocity. Oil temperature is assumed equal to water temperature, neglecting solar radiation or any other energy transfer process that may influence oil temperature.

## Equations

### Spreading

### Density

### Viscosity

### Evaporation

### Emulsification

### Dispersion

### Sedimentation

### Dissolution

### Oil-Beaching

### Removal techniques

#### Chemical Dispersion

#### Mechanical Cleanup

## User manual

### Main options

### Input data file

## References
