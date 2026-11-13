# Authors and contributions

[Mohid Water Modelling System](../overview/README.md#mohid-water-modelling-system#mohid-water-modelling-system) is currently developed by [MARETEC](http://www.maretec.mohid.com/), a research group at [Instituto Superior Técnico](http://www.ist.utl.pt/), Technical University of Lisbon, Portugal, in cooperation with ([Hidromod](http://www.hidromod.pt/)).

The architecture of the model is presently coordinated by Ramiro Neves, Paulo Leitão ([Hidromod](http://www.hidromod.pt/)), Frank Braunschweig and Luis Fernandes. The engineering of the model is mainly a responsibility of Ricardo Miranda (software concepts and architecture), Paulo Leitão (numerical algorithms), Frank Braunschweig (software architecture / implementation and GUI development) and Luis Fernandes (software architecture / implementation and GUI development).

Mohid's modular structure enables a large team of developers working simultaneously on the same source code version. Currently the team is constituted by Pedro Chambel Leitão (soil hydrodynamics), David Brito (vegetation, property transport and transformation), Sofia Saraiva (biogeochemical models), Rodrigo Fernandes (operational model, oil spills), Guillaume Riflet (hydrodynamics), Angela Canas (data assimilation and [Mohid Water support tools](../tools/mohid-support-tools.md)), Pedro Pina (biogeochemical models), Pedro Galvão (biogeochemical models, soil hydrodynamics, GUI development), Henrique Coelho [Hidromod](http://www.hidromod.pt/) (turbulence), Marcos Mateus (biogeochemical models), João Nogueira (data assimilation), Ana Rosa Trancoso (river hydrodynamics, biogeochemical models) and Ricardo Miranda.

Contributions from Hernani Theias (non-hydrostatic processes), Miguel Carmo (sand transport formulations), Manuel Villarreal, Pedro Montero and Pablo Carracedo García are also included in Mohid.

The main contributors for the concepts included in the model and to its development are Ramiro Neves, Adélio Silva, José Leitão, Flávio Martins and Aires dos Santos.

A lot of other people gave contributions that all together have strongly influenced the development of Mohid.

## Source code from other projects included in Mohid

> -   [GOTM](http://www.gotm.net/) - turbulence model library included in [Module Turbulence](../water/modules/turbulence.md)
> -   WGS84 to NLRD coordinate conversion C code written by E.J.O Schrama 

$$
e.j.o.schrama@lr.tudelft.nl
$$

 was ported to FORTRAN95 into [Module HorizontalGrid](../water/modules/README.md#module-horizontalgrid#module-horizontalgrid).

## Software used in Mohid

> -   [HDF5](http://hdf.ncsa.uiuc.edu/HDF5/)
> -   [MPI](http://www.mcs.anl.gov/mpi/)
