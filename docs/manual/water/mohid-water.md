# Mohid Water

## Overview

MOHID Water consists of a three-dimensional numerical model to simulate surface water bodies such as rivers, reservoirs, estuaries, coastal areas or the ocean. It is one of the main programs in [MOHID Water Modelling System](../overview/README.md#mohid-water-modelling-system#mohid-water-modelling-system), written in FORTRAN 95 using an [object oriented programming](../development/programming/object-oriented-programming-in-mohid.md) phylosophy. It is composed of a series of [modules](modules/modules.md) built on top of [Mohid Base 1](README.md#mohid-base-1#mohid-base-1) and [Mohid Base 2](README.md#mohid-base-2#mohid-base-2). These modules are mainly responsible for computing physical or biogeochemical process, e.g. [Module Hydrodynamic](modules/hydrodynamic.md), [Module WaterProperties](modules/waterproperties.md), [Module Turbulence](modules/turbulence.md), [Module Lagrangian](modules/lagrangian.md), etc.

MOHID Water was designed in order to be able to simulate aquatic systems dividing them into three compartments or media: air, water and land. Thus it was constructed assuming: one model consisting of two main interfaces: the water-sediment interface and the water-air interface, dividing three well defined compartments, the atmosphere, the water column and the sediment. The two interfaces should be able to communicate by handling the fluxes between the three compartments. To do this, two modules were created: [Module InterfaceSedimentWater](modules/interfacesedimentwater.md) and [Module InterfaceWaterAir](modules/interfacewaterair.md).

Thus, the model was set up in a way that each of the three compartments did not have access to another compartment, being all communications made by the interface modules. That is, with this organization, [Module SedimentProperties](modules/sedimentproperties.md) and [Module Consolidation](modules/consolidation.md) (representative of the sediment column), communicate only with [Module InterfaceSedimentWater](modules/interfacesedimentwater.md), as well as modules Hydrodynamic, WaterProperties, Turbulence and Lagrangian (representative of the water column). These modules also communicate with [Module InterfaceWaterAir](modules/interfacewaterair.md) which provides the connection with [Module Atmosphere](modules/atmosphere.md). [Module Model](modules/README.md#module-model#module-model) controls and manages all the information fluxes.

## Physics

[Mohid Water Physics](mohid-water-physics.md)

## Biogeochemistry

[Mohid Water Biogeochemistry](../overview/README.md#mohid-water-biogeochemistry#mohid-water-biogeochemistry)

[Heavy metals](../reference/heavy-metals.md)

## Sediments

[Module Sediment](modules/sediment.md)

## Other features

[Parallel processing](../development/programming/parallel-processing.md)

## Support Tools

[Mohid Water support tools](../tools/mohid-support-tools.md)

## Others

[Mohid Water software](../reference/mohidwater.md)

[Modeling the water properties transport - a step by step methodology](modeling-waterproperties-transport.md)

[Relaxation](relaxation.md) - using [Module Assimilation](modules/assimilation.md)

[Nesting models](../reference/nesting-models.md)

[New: MOHID Water Vectorial and Angle Properties](../reference/new-mohid-water-vectorial-and-angle-properties.md)
