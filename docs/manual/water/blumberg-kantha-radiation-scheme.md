# Blumberg-Kantha radiation scheme

## Overview

## How to use in MOHID

### Input file

In order to activate the Blumberg-Kantha radiation scheme keyword

```
RADIATION   : 3
TLAG_FILE   : ..\..\GeneralData\Tlag.dat
```

must be defined. In this file, the relaxation times used in the radiation boundary conditionmust be specified, one time for each open boundary cell where one wants to consider this boundary condition.

The general structure of the file must be a block (**<TlagBegin>**...**<TlagEnd>**) containing the cells i and j location and the relaxation time in seconds. One line per cell, as in the following example.

```
<TlagBegin>
1   1   200
1   2   200
1   3   200
1   4   200
1   5   200
1   6   200
1   7   200
1   8   200
1   9   200
...
<TlagEnd>
```

### Relaxation times

The relaxation time considered for each cell acts like a weight in the ponderation of the internal solution of the domain (obtained by model dynamics) and the reference solution one considers at the open boundary (usually tide).

This makes that the relaxation time must be chosen according with the domain horizontal grid resolution, bathymetry characteristics, modelled dynamics and the reference solution characteristics.

No deterministic rule for the calculation exists: the relaxation times have to be specified considering each modelling case. Nevertheless, some considerations can help the specification of these times:

- A big relaxation time means that the reference solution is being given a smaller weight than in the case of a small relaxation time;
- Relaxation time should be higher in low depth cells than in high depth cells: i.e. near coast one expects higher relaxation times than in the open deep ocean; this consideration happens to take into account the cell dimension available to accomodate the perturbation induced by the reference solution (one wants to avoid instabilities) and specially to allow the noise which propagates along the coast as Kelvin waves to escape the domain;

In the experience obtained with MOHID, usually an approximate relation of 10:1 leads to good results: e.g. 3600s near coast and 360s in the deep ocean.

There must not be an abrupt change of the relaxation time considered in adjacent cells, to avoid model instabilities; this is particularly critical in low depths and small horizontal resolution cells.

The ideal is to considered the same relaxation time in most of the boundary. However if near coast and deep ocean situations exist in the boundary, the account of the previous considerations requires the existence of a transition area, where relaxation times change.

To allow an uniform relaxation times situation in the coast, it is valuable to consider a set of several cells near the coast with the same relaxation time. After this set, the relaxation times change till the deep ocean value is reached. The lenght (in cell number) of the coastal buffer and the transition area must be chosen according with the horizontal grid resolution and depths, to allow a smooth transition.

The specific values of relaxation times used are set usually by experience. It is very valuable to take into account the values used in modelling experiences which provided good results in similar domains. These values can be used as a first approach and according with the results obtained the relaxation times set can then be improved.

As an example, in a domain containing all Portuguese Continental and Galician coast, with a horizontal resolution in the boundary of about 4km, is considered 1800s in the coast and 200s in the deep ocean. The 1800s value is maintained for the 11 cells nearer the coast and then the relaxation time is gradually reduced in the next 15 cells till the 200s value is reached. The depth near coast is about 30m and in the maximum depth in the opean ocean is about 5000m. This relaxation time scheme provides acceptable results in a 2D barotrophic model with tide and wind.

Is highly recommended that the initial values considered for the relaxation times to be tested and improved till the required model results quality standards are reached, since these are very specific of each modelling study. It is not at all advisable to use as definite a set of relaxation times defined for another study, even if very similar, without assessment of the results.

## References

Blumberg, A.F. & Kantha, L.H., Open boundary condition for circulation models, Journal of hydraulic engineering, 111(2):237-55, 1985
