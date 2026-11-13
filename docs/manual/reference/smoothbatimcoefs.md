# SmoothBatimCoefs

This program creates smoothing coefficients for a bathimetry to allow a [smooth transition](../tools/smoothbatimnesting.md) in the boundary between the the bathymetry and a larger scale one. See also [SmoothBatimNesting](../tools/smoothbatimnesting.md).

## Quick start

1. Create the son bathymetry (you may use [Mohid GIS](../tools/gui/mohid-gis.md)).
2. Create and edit the options *SmoothBatimCoef.dat* file and save it in the same path as the executable.
3. Run the executable.

## Options file

Here's what the options file *SmoothBatimCoef.dat* looks like:

```
SON_BATIM                   : ..\GridData_2.dat
NEW_COEF                    : SmoothCoefs.dat
NUMBER_CELLS                : 10
```

The *Number\_cells* are the number of grid points taken from the border in, that will yield non-null smoothing coefficients.

## Algorithms

So far there is only one type of algorithm for the smoothing coefficient $C_{i\,j}$ defined [here](../tools/smoothbatimnesting.md). It's basically a quadratic growth from 0 to 1 from the border to *Ncells* inside, where *Ncells* is a parameter given by the modeller. Here's the formula:

- *ILB*: lower bound in the *i* direction.
- *IUB*: upper bound in the *i* direction.
- *JLB*: lower bound in the *j* direction.
- *JUB*: upper bound in the *j* direction.

$\forall\;_{i\,j}\quad C_{i\,j} = \begin{cases} \frac{(ILB-i)^2}{Ncells^2}, & if\; (ILB-i)^2 < Ncells^2 \\ \frac{(IUB-i)^2}{Ncells^2}, & if\; (IUB-i)^2 < Ncells^2\end{cases}$

$\begin{cases} \frac{(JLB-j)^2}{Ncells^2}, & if\; (JLB-j)^2 < Ncells^2 \\ \frac{(JUB-j)^2}{Ncells^2}, & if\; (JUB-j)^2 < Ncells^2 \\ 1., & otherwise \end{cases}$

### Future algorithms

Please place here the algorithms you'd wish be implemented:

- ...
