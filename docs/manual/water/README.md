# MOHID Water

The three-dimensional hydrodynamic and water-quality model for estuaries, coastal seas and reservoirs.

## Sections

- [Modules](modules/README.md) (41 pages)

## Pages

- [Adsorption/Desorption](adsorption-desorption.md) — Adsorption relates to the process where a solute in a liquid phase becomes bonded to the surface of a solid (Linde, 2002). It can occur in three major…
- [Advection numerical scheme](advection-numerical-scheme.md) — 1 - First order upwind is the more robust but the less precise;
- [Atmospheric boundary condition](atmospheric-boundary-condition.md) — For each water property the user needs to define if it want or not to impose atmospheric fluxes. Here's an example of atmospheric properties that need…
- [Automatic Data Interpolation](automatic-data-interpolation.md) — A model is setup to run every week. However he will OBCs and initial conditions for u, v, S, T and waterlevel. Thus, the external data needs to be int…
- [Blumberg-Kantha radiation scheme](blumberg-kantha-radiation-scheme.md) — In order to activate the Blumberg-Kantha radiation scheme keyword
- [Bottom boundary condition](bottom-boundary-condition.md) — For each water property the user needs to define if it wants or not to impose sediment-water fluxes of particulate or dissolved matter. An example of…
- [Coupling Water-Atmosphere User Manual](coupling-water-atmosphere-user-manual.md) — This user manual intends to help MOHID users to construct properties in the atmospheric and interface water-air compartments, in a way that his modeli…
- [Diffusion values](diffusion-values.md) — By default the horizontal turbulent viscosity is assumed equal to the horizontal turbulent diffusion (or diffusivity). In the vertical, the same is as…
- [Energy](energy.md) — Kinetic energy
- [Evolution](evolution.md) — There are several ways of running MOHID water hydrodynamics. The keyword that controls these several ways is EVOLUTION from the module Hydrodynamic fi…
- [Hydrolight](hydrolight.md) — Hydrolight is a radiative transfer numerical model that solves the time - independent radiative transfer equation. In order to compute the radiance, H…
- [ImposedSolution](imposedsolution.md) — The idea behind the imposed solution, is the ability to load pre-calculated hydrodynamic and waterproperties values from hdf5 files. Typically is used…
- [Mercator Automatic Data Acquisition](mercator-automatic-data-acquisition.md) — This article describes the automatic data acquisition system built to get large-scale circulation model gridded data from Mercator-Océan ftp server. A…
- [Model setup methodology](model-setup-methodology.md)
- [Modeling waterproperties transport](modeling-waterproperties-transport.md) — Check the bathymetry:
- [Mohid Water](mohid-water.md) — MOHID Water consists of a three-dimensional numerical model to simulate surface water bodies such as rivers, reservoirs, estuaries, coastal areas or t…
- [Mohid Water Physics](mohid-water-physics.md) — Module Lagrangian
- [MPIConfig.dat](mpiconfig-dat.md) — The MPIConfig.dat file is no longer used in MPICH 1.2.x. You must write a command line with every option you intend to use.
- [Open Boundary conditions](open-boundary-conditions.md) — For each property it's possible to define a specific boundary condition. The available options are:
- [OpModel TimeLine](opmodel-timeline.md) — This article describes the Mercator-IST operational system timeline.
- [Point Origin](point-origin.md) — To define a point origin in the lagrangian module you need to define the coordinates where the tracers are to be released. This can be done in several…
- [Relaxation](relaxation.md) — Relaxation schemes are used in MOHID Water to relax computed solutions to a reference solution. The standard method follows the equation:

## Short entries

### Main

#### Overview

#### Sequence

##### Simple simulations

##### Parallel nested simulations

#### References

### Modelacao ambiental

#### Links exteriores

- [Página modelação ambiental IST](https://fenix.ist.utl.pt/homepage/ist146730/modelacao-ambiental-2007-2008)

### Mohid Base 1

MOHID Base 1 is a set of modules in [Mohid Framework](../development/programming/README.md#mohid-framework#mohid-framework) that handle specific operations used by several numerical models and tools in MOHID Water Modelling System.

### Mohid Base 2

MOHID Base 2 is a set of modules in [Mohid Framework](../development/programming/README.md#mohid-framework#mohid-framework) that manage gridded information in MOHID Water Modelling System.
