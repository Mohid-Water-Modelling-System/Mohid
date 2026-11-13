# Module SedimentProperties

## Overview

Module SedimentProperties is responsible by handling the sediment compartment. The sediment column consists of a saturated porous media, formed by sediments and by water that fills the interstices between the sediments. Properties in this compartment can also be dissolved - in the porewater-, or particulate – adsorbed on to sediments.

The sediment column model is basically a set of 1D vertical models defined below the 3D water column model. Both models share the same horizontal discretization, but compute independent vertical coordinates.

## Dissolved properties

Transport of dissolved properties in porewater is computed only in the vertical axis, as horizontal gradients are not considered in this study. Therefore, the transport equation can be written:

```
Equation
```

Where Cd is the concentration (kg.mwater-3), w (m/s) is the porewater velocity due to compaction, kz (m2s-1) is the diffusivity coefficient. The molecular diffusion coefficients must be corrected with tortuosity parameterization, to account for the increase in the solute pathways due to difficulty presented by the sediment particles for diffusion to occur. Two different formulations were included in model, the first following formulation by Berner (1980), tortuosity dependent, which on the other hand is porosity dependent:

```
Equation
```

The second by Soetaert (1996), is dependent of porosity square.

```
Equation
```

Where, DINF (m2s-1) is the molecular diffusion coefficient in solution, and Dm (m2s-1) is the corrected molecular diffusion coefficient.

```
Equation
```

Bioturbation is computed as a diffusion coefficient, which is present until a certain depth, and decreases exponentially with it. This pretends to simulate benthic fauna activity, which is most of the times present in sediments upper 10-15cm. Below this, bio-activity can be considered negligible. Thus, the bioturbation diffusion coefficient, Db (m2s-1), can be computed by:

```
Equation
```

In which zb is the depth limit for maximum biological activity and α is a decay coefficient (m) to account for the decrease of bioturbation with depth.

```
Equation
```

Thus, the diffusivity coefficient, kz, becomes the sum of the molecular and bioturbation diffusion coefficients:

## Particulate properties

Particulate properties (kgproperty/kgsediment) vary in time due to sinks and sources, namely adsorption/desorption, and due to bioturbation mixing effect.

## Adsorption/Desorption

[Adsorption/Desorption](../adsorption-desorption.md) processes are simulated with a similar approach as in the water column.
