# Energy

**Energy** calculates a domain averaged timeseries of the following diagnostic properties:

Kinetic energy 

is the mass with squared velocity modulus product integration.

Potential energy 

is the mass with gravitic constant and hydrographic zero distance product integration.

Enstrophy 

Is the squared velocity curl integral.

Volume 

is the integral of all wet [T-cells](../base/module-geometry.md) volumes.

Mass 

is the density with [T-cells](../base/module-geometry.md) volume product integration.

The diagnostic properties computation can be found in the routine **ComputeSystemEnergy** inside the [Module Hydrodynamic](modules/hydrodynamic.md). All the properties were interpolated to type [T-Cells](../base/module-geometry.md).

## Usage

```
Energy: 1
```

## See also

- [Module Hydrodynamic](modules/hydrodynamic.md)
