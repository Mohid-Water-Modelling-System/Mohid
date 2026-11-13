# Relaxation

Relaxation schemes are used in [MOHID Water](mohid-water.md) to relax computed solutions to a reference solution. The standard method follows the equation:

$\frac{\partial \phi}{\partial t} =  \frac{ \phi_{ref} - \phi}{\tau}$

where $\phi$ is the relaxed variable, $\phi_{ref}$ is the reference solution and $\tau$ relaxation time scale.

In MOHID Water, the reference fields and the relaxation time scales are read and handled by [Module Assimilation](modules/assimilation.md). The relaxation algorithm itself is implemented in [Module Hydrodynamic](modules/hydrodynamic.md) and [Module WaterProperties](modules/waterproperties.md) to apply the relaxation scheme to properties such as velocity U, velocity V, water level, temperature, salinity, etc. It is normally applied in the boundary condition computation, usually when [nesting models](../reference/nesting-models.md) or using a climatological reference field, via the definition of a "sponge layer" around the open boundary. However it can also be applied to the entire simulation domain.
