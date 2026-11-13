# Upwelling

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Upwelling is quite well described [here](http://www.pfeg.noaa.gov/products/PFEL/modeled/indices/upwelling/NA/what_is_upwell.html).

## Upwelling index

The upwelling index (*UI*) is defined as the *Ekman Transport per 100 m of coastline* induced by wind parallel to the coast. In the Western Iberian coast. This corresponds to the zonal component of the wind.

The volume Ekman transport is defined in *Kundu* as the momentum of a steady-state horizontally homogeneous viscid geophysical flow integrated over depth. Thus, assuming a steady state constant zonal wind stress $\tau$

$UI=-100\,\frac{\tau}{\rho\,f}$ m^3^/s/100m.

### Wind stress

To calculate the wind stress $\tau$ from the wind velocity modulus $U_{10}$ (air speed at 10 m height from the free surface) *Ekman* used the bulk formula:

$\tau=\rho_{air}C_D\,U_{10}^2$

$\rho_{air}=1.25$ *kg/m^3^* is an usual value for the air density.

$C_D$ is an air drag coefficient. It is mainly dependent from surface rugosity length *z~0~* and from the Richardson number. However a large case of linear parameterizations is available throughout the litterature, with no well defined standard:

$C_D = a + b\,U_{10}$

- *Large & Pond 1981* use *C~D~ = .44 + .63 U~10~*,
- *Smith & Banke 1975* use *C~D~ = .63 + .66U~10~*.

## Viscid, horizontally homogeneous, steady-state, geophysical fluid equations

The flow equations of motion are given in *Kundu* by

$-fv=\nu \frac{d^2u}{d z^2},$

$fu=\nu \frac{d^2v}{d z^2}.$

## The Ekman spiral

*Kundu* shows that when acted by a constant wind stress $\tau$ along the *x*\-direction, the solution to the above equations yield

$u=\frac{\tau/\rho}{\sqrt{f\nu}}e^{z/\delta}\cos\left(-\frac{z}{8}+\frac{\pi}{4}\right),$

$v=-\frac{\tau/\rho}{\sqrt{f\nu}}e^{z/\delta}\sin\left(-\frac{z}{8}+\frac{\pi}{4}\right),$

where $\delta\equiv\sqrt{\frac{2 \nu}{f}}.$

The velocity vector rotates clockwise (looking down) with dept. Its magnitude decays exponentially as well with decay constant $\delta$. This spiral shaped flow is called the *Ekman spiral*. In the works of Ekman, Ekman defined a thickness of decay of the spiraled flow based on the flow geometry. This thickness later defined the *Ekman layer* and is usually around 50 m for the open ocean(*Stewart*). *Kundu* however defines the thickness of the Ekman layer simply by its decay constant $\delta$.

## The Ekman transport

When acted by a constant wind stress $\tau$ along the *x*\-direction, the volume transport in the Ekman layer is given by

$\int_{-\infty}^{0}u\,dz=0$,

$\int_{-\infty}^{0}v\,dz=-\frac{\tau}{\rho\,f}$.

Thus the *net transport is to the right of the wind stress in the northern hemisphere* and is independent of viscosity. This is due because near the surface $\rho \nu \frac{du}{dz}=\tau,\quad z=0$. Over statistical flow averages, the *Ekman transport* calculation fits very well with experimental data and, thus, is considered a reliable calculation.

## References

1. *Kundu*
2. *Stewart*
3. *Pietrzak2002*
4. *Large & Pond 1981*
5. *Smith & Banke 1975*
