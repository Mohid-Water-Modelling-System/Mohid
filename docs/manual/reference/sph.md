# SPH

**Smoothed-Particle Hydrodynamics** is a computational technique to simulate fluid dynamics using a mesh-free lagrangian-tracer approach. The mesh-free lagrangian tracers can represent any property, from granulate to continuum in nature, from sand particles, to density, to fluid momentum and vorticity. **SPH** estimates any property value on any point of the domain based on the lagrangian tracers in the vicinity, weighted by a kernel function. The lagrangian tracers, in turn, are transported by the very same velocity field they generate. They are also subjected to any physical process whose continuum field they previously generated. This allows SPH to generate highly non-linear continuum fields, most unlike the traditional structured-mesh algorithms, and it can realistically simulate, for example, wave breaking, flooding, two-phasic flow interaction.

Because of the lagrangian-tracer nature of SPH, it has been successfully parallelized in both MPI and GPU architectures, with performance increases in the order of 100 for the latter; thus making this technique more and more interesting for realistic fluid simulation scenarios.

The mesh-free nature of SPH allows it simulate extremely turbulent and non-linear flow, filling and emptying extremely complex volumes.

Both of these interesting features, many million particle simulation at competitive performance thanks to GPU computind and non-linear flow filling and emptying very complex volumes, makes SPH very attractive as a service to maritime industry. Some companies, such as HydrOcean, estimate an increase of its use, at the request of clients in the industry, in the forthcoming years.

## History

## State-of-the-art

## External references

- [Smoothed-Particle Hydrodynamics on Wikipedia](http://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)

## Major publications on SPH
