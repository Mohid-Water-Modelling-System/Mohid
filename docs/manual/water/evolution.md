# Evolution

There are several ways of running MOHID water hydrodynamics. The keyword that controls these several ways is **EVOLUTION** from the [module Hydrodynamic](modules/hydrodynamic.md) file.

Here's the list of possible options for this keyword:

1. Solve\_Equations : Default value, it solves the Primitive Ocean Equations.
2. Read\_File : Reads the solution from a file
3. No\_hydrodynamic : Doesn't calculates the hydrodynamic at all
4. Residual\_hydrodynamic : ??
5. Run\_Off : ??
6. [ImposedSolution](imposedsolution.md) : Uses a reference solution for U, V and level.
7. Vertical1D : Used when modelling a simple 1D vertical model. It takes into account the Ekman spiral, provided the model is forced with wind stress at the surface, or with friction at the bottom.

## Example

```
EVOLUTION : Vertical1D
```

## External links

- [Mohid Keywords database](http://www.mohid.com/IHS)
