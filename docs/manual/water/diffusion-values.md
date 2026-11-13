# Diffusion values

By default the horizontal turbulent viscosity is assumed equal to the horizontal turbulent diffusion (or diffusivity). In the vertical, the same is assumed except if the turbulent model (GOTM) option is on.

## Turbulence\_x.dat

By default it is assume a constant turbulent horizontal viscosity (m2/s)

```
    VISCOSITY_H             : 10
```

By default it is assume a constant turbulent vertical viscosity (m2/s)

```
    VISCOSITY_V             : 1e-3
```

The user can activate the turbulence model (GOTM package). In this case the user must look to the turbulence module output results to see the fields of vertical turbulent viscosity and diffusivity.

```
    MODTURB                 : turbulence_equation
```
