# Advection numerical scheme

## WaterProperties\_x.dat

1 - First order upwind is the more robust but the less precise;

2 - Second order upwind has good precision but is highly unstable;

3 - Third order upwind has excellent precision, more or less stable and not;

4 - The Transport Variation Diminishing method(TVD) is in fact a hydrid scheme between a first order and a third order method where a ponderation factor is computed so the method is able to insure positiveness. The method 4 is advised for experienced users. By default the method is 1 to increase robustness of the runs.

5 - Central differences has good precision but is highly unstable

Horizontal keyword option

```
ADV_METHOD_H  : 4
```

Vertical keyword option

```
ADV_METHOD_V  : 4
```

There are several methodologies to compute the TVD ponderation factor. MinMod = 1, VanLeer = 2, Muscl = 3, Superbee = 4, PDM = 5. By default the option admitted is SuperBee.

Horizontal keyword option

```
TVD_LIMIT_H : 4
```

Vertical keyword option

```
TVD_LIMIT_V : 4
```
