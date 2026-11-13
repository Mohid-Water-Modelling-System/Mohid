# Module Consolidation

## Overview

## Sediment physical processes and properties

The sediment compartment is constituted by a module which computes the sediment geometry (variations due to erosion and consolidation), namely dry sediment volumes and interstitial water volumes. In terms of vertical referential, it is located below the water column until a certain defined depth.

The construction of the domain is made by means of a depths file, similar to the bathymetry for the water column model. This way, sediments’ upper layer is located at the same coordinate of the water column model bathymetric value, with a certain depth, usually 10 to 30 cm.

This compartment is considered to be a saturated porous media, so a key variable is porosity (Ф), which represents the fraction of volume occupied by interstitial water. Porosity decreases with depth and relates to tortuosity, a parameter which reflects the influence of porous media geometry in the transport phenomena, namely diffusion. Tortuosity can be seen as an extension of the path a solute has to take in the porewater, due to the fact that, it has to follow a complex structure of micro-channels in the available spaces between sediment particles. Boudreau (1996), finds a good agreement between tortuosity and porosity:

```
Equation
```

A decay of porosity can be computed, accounting for the consolidation process.

```
Equation
```

Where Ф∞ is the porosity of a fully consolidation sediment and λ is decay factor (s). This consolidation process has a time scale several times higher than the erosion/deposition processes and in this study is neglected. However, it is included in the model, and can be useful in long term simulations, has when consolidating, interstitial water is pushed upwards, therefore advecting solutes through the sediment column, and even through the water-sediment interface onto the overlying water column. These fluxes can also be accounted as a source of contaminants to the water column.

The sediment compartment boundary conditions consist on the erosion and consolidation fluxes, and are controlled by the [sediment-water interface module](interfacesedimentwater.md). Erosion is made, by removing material from the sediments’ upper layer. As sediment layers are being scoured, critical shear stress increases, due to the fact that sediments compaction level increases with depth. Therefore, critical shear stress can be computed such as:

```
Equation
```

Where z is the depth (m) and Ψ is a decay coefficient (m). A specific new algorithm was developed to solve discretization problems of a complex vertical domain, like the sediment compartment. The vertical resolution must be high enough to solve properly the sharp concentration gradients (contaminants, organic matter, oxygen, etc) existing in estuarine sediments. Two main problems can be found: the sediment top layer is constantly eroded until it disappears; or the deposition flux is so high that the top layer thickness increases to a level that it cannot be assumed that properties inside the layer are constant. Thus, in order to handle these problems two thickness limitations were imposed: a minimum and a maximum layer thickness.

When erosion fluxes remove material from the sediments’ compartment upper layer, this flux is limited so that, in one iteration, the layer does not exceed the minimum layer thickness. When this happens, the upper layer collapses and becomes part of the lower layer, which then, becomes the top layer. When consolidation fluxes raise the top layer thickness so that it exceeds the maximum layer thickness, a new layer is created, splitting the upper layer into two. The new upper layer is initialized has having the minimum layer thickness allowed.

To overcome these problems, a new vertical coordinate system was created to account for collapsing and splitting of layers. A two-dimensional mapping variable monitors which is the index of the top layer, above which, all water and sediment volumes are null, as well as all processes. The model must always be started with a certain number of empty top layers to account for possible creation of new layers, if consolidation occurs. If the initial number of layers is exceed, the model stops. The same happens when all sediment layers are eroded and collapsed. The layers collapsing and splitting is followed by mass conserving algorithms applied to each of the sediment properties, both dissolved and particulate.
