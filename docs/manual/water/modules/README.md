# Modules

Optional physical, biogeochemical and numerical modules that extend the base model.

## Pages

- [Module Assimilation](assimilation.md) — Module Assimilation is responsible for handling information used in relaxation schemes. It mainly reads the reference fields and the relaxation coeffi…
- [Module Atmosphere](atmosphere.md) — The atmosphere module is responsible for meteorological data needed to compute processes occurring at the water-air interface, such as computing wind…
- [Module Benthos](benthos.md) — Module Benthos is a zero-dimensional biogeochemical model and was designed to compute, using simple model formulations, some of the main benthic bioge…
- [Module CEQUALW2](cequalw2.md) — Module CEQUALW2 is a module using the CE-QUAL-W2 ecological model formulations. CE-QUAL-W2 is 2D laterally average hydrodynamic and ecological model d…
- [Module Consolidation](consolidation.md) — The sediment compartment is constituted by a module which computes the sediment geometry (variations due to erosion and consolidation), namely dry sed…
- [Module Discharges](discharges.md) — Module Discharges handles all the point discharges in Mohid Water Modelling System. These discharges can be simple water discharges but can also have…
- [Module DrainageNetwork](drainagenetwork.md) — This module can be used by two models: MOHID Land, and MOHID River Network. Like other modules, it has a specific input file, called DrainageNetwork\X…
- [Module FreeVerticalMovement](freeverticalmovement.md) — The free vertical movement module computes particulate water properties vertical fluxes. It is normally used to computed settling velocity in cohesive…
- [Module GOTM](gotm.md) — Module GOTM is a "wrap-up" module containing GOTM routines which were coupled into MOHID in 2001. This routines consist of a set of turbulence models…
- [Module Hydrodynamic](hydrodynamic.md) — The hydrodynamic model solves the primitive continuity and momentum equations for the surface elevation and 3D velocity field for incompressible flows…
- [Module HydrodynamicFile](hydrodynamicfile.md) — Module HydrodynamicFile can be seen as an auxiliary module, which allows the user to integrate the hydrodynamic solution in space and in time and stor…
- [Module InterfaceSedimentWater](interfacesedimentwater.md) — The sediment-water interface module computes boundary conditions at the bottom of the water column. It computes shear stress as a boundary condition t…
- [Module InterfaceWaterAir](interfacewaterair.md) — The water-air interface module is responsible by processes occurring at the water-air interface, such as computing wind shear stress, radiation balanc…
- [Module Jet](jet.md) — Module JET, also known as MOHIDJET integral model, aims to simulate the initial dilution associated to outfalls jets. The model is used as an initial…
- [Module Lagrangian](lagrangian.md) — Lagrangian transport models are very useful to simulate localized processes with sharp gradients (submarine outfalls, sediment erosion due to dredging…
- [Module LagrangianGlobal](lagrangianglobal.md) — The LangrangianGlobal Module is a deep upgrade of the Lagrangian Module. These upgrade has three main goals:
- [Module Life](life.md) — Module Life is a multiparameter biogeochemical module based on ERSEM. This is a more complex model compared to Module WaterQuality, which is able to s…
- [Module MacroAlgae](macroalgae.md) — Module MacroAlgae is a zero-dimensional biogeochemical model and was designed to compute the main biogeochemical processes related with macroalgae in…
- [Module Oil](oil.md) — The prediction and simulation of the trajectory and weathering of oil spills are essential to the development of pollution response and contingency pl…
- [Module Reservoirs](reservoirs.md) — Module Reservoirs allows to compute reservois dynamic that change river flows and concentrations. Module Reservoirs recieves flows and concentrations…
- [Module Sand](sand.md) — Module Sand was the result of the accumulated experience gain by Hidromod in the framework of many engineering projects [
- [Module Sediment](sediment.md) — Module Sediment allows simulating the dynamics of sediments in natural systems under the action of waves and currents. It allows to assess the transpo…
- [Module SedimentProperties](sedimentproperties.md) — Module SedimentProperties is responsible by handling the sediment compartment. The sediment column consists of a saturated porous media, formed by sed…
- [Module SedimentQuality](sedimentquality.md) — This is the 0D module that handles sediment biological processes made by bacteria as organic matter mineralization, nitrification, denitrification and…
- [Module Statistics](statistics.md) — It is better described in the support tool HDF5Statistics, but can be used in the MohidWater program.
- [Module TimeSerie](timeserie.md) — Module TimeSerie manages all operations involving Time Series files, namely in terms of input and output.
- [Module Turbine](turbine.md) — Module Turbine is in charge of calculating the force exerted by turbines in the water flow as well as the power and energy extraction. This implementa…
- [Module Turbulence](turbulence.md) — Burchard et al, 1999 – (http://www.gotm.net), a turbulence models library coupled with MOHID, including a k-ε model and Mellor-Yamada second order tur…
- [Module WaterProperties](waterproperties.md) — Module WaterProperties is the 3D eulerian transport module included in MOHID. Module WaterProperties is responsible for computing the properties evolu…
- [Module WaterQuality](waterquality.md) — Water Quality diagram that illustrates the document intitled MOHID WaterQuality Module - Technical Guide - September 2005.
- [Module Waves](waves.md) — Module Waves is responsible to simulate surface waves in MOHID. It can act as a wave database, reading wave data from files, or it can predict wave pa…
- [ModuleNETCDF](modulenetcdf.md) — ModuleNETCDF is an interface module between the netcdf libraries and Mohid Framework, and it was developed in order to read and write netcdf files in…
- [Modules](modules.md)

## Short entries

### Module AdvectionDiffusion

#### Overview

#### Advection

##### First order upwind

##### Second order upwind

##### Third order upwind

##### Total Variation Diminishing (TVD)

##### Central differences

#### Diffusion

#### Other remarks

### Module Bivalve

For more information please download the Module Bivalve Manual 

$$
[Download](http://www.mohid.com/PublicData/Products/Manuals/MOHID_Module_Bivalve_Manual.pdf)
$$

### Module GridData

Module GridData handles all information regarding reading and loading [Grid Data](../../io/formats/grid-data.md) files in MOHID, including [Bathymetry](../../io/formats/bathymetry.md) or Topography files.

### Module HorizontalGrid

#### Overview

Module HorizontalGrid manages all the information in MOHID regarding the horizontal discretization of a simulation domain. This information is read from a [grid](../../io/formats/grid.md) file.

#### Types of coordinates supported

##### Geographic

##### UTM

##### Portuguese Military

##### Simple geographic

##### Non-referenced

##### Circular

##### Netherlands RD

#### Coordinates conversion

#### References

### Module Model

#### Overview

Module Model handles all the information of a single model (in [Mohid Water](../mohid-water.md)). It is responsible for constructing, modifying and destructing each model and for controlling information fluxes between the different modules, and time and mapping evolution.

#### See also

[Choose the model time step](../../reference/README.md#choose-the-model-time-step#choose-the-model-time-step)

#### References

### Module Shell

Module Shell is the shell of a standard Mohid module, consisting of an almost empty module, which contains only module management code, used to fulfill the object oriented philosophy in which [Mohid Framework](../../development/programming/README.md#mohid-framework#mohid-framework) is based.

### Module Time

#### Overview

Module Time handles the time functions in MOHID.

#### See also

[MOHID Time format](../../io/formats/README.md#mohid-time-format#mohid-time-format)

### Module TurbGOTM

#### Overview

Module TurbGOTM is an interface module used to wrap [Module GOTM](gotm.md) and couple it to [Module Turbulence](turbulence.md). As [GOTM](http://www.gotm.net) original code is programmed using a different philosophy than MOHID's, this module was designed in order to communicate with [Module Turbulence](turbulence.md) in a "MOHID compliant way".
