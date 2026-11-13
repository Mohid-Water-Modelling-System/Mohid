# Input and output

Every input and output file format used by the models, plus the keyword database.

## Sections

- [Formats](formats/README.md) (13 pages)
- [Keywords](keywords/README.md) (11 pages)

## Pages

- [Filter.dat](filter-dat.md) — Input Bathymetry in Mohid Grid Data format
- [Free Vertical Movement.dat](free-vertical-movement-dat.md)
- [Input Data Errors](input-data-errors.md) — User input data errors in Mohid can be estimated coarsely in near 90% of total execution errors, and can have two origins: user distraction or user un…
- [InterfaceSedimentWater.dat](interfacesedimentwater-dat.md)
- [Land boundary condition](land-boundary-condition.md) — The land boundary presently in the MOHID system is imposed in the form of water and mass fluxes.
- [Land boundary condition;](land-boundary-condition.md) — The land boundary presently in the MOHID system is imposed in the form of water and mass fluxes.
- [Nomfich.dat file](nomfich-dat-file.md) — A nomfich.dat file is a file commonly used in the MOHID Framework to indicate to a MOHID program the name and path of the input (instructions) file(s)…
- [Outputs.dat](outputs-dat.md) — SWAT\WITH\CHANGES  : 1
- [Sediment.dat](sediment-dat.md)
- [Tree.dat](tree-dat.md) — The Tree.dat file specifies the hierarchical structure of a Mohid Water nested models simulation. This file is expected to be placed in the working di…
- [WaterProperties.dat](waterproperties-dat.md)

## Short entries

### Boundary conditions input examples

### Imposing baroclinic radiation

Hydrodynamic\_x.dat

```
!ON = 1, OFF = 0
BAROCLINIC            : 1
!ON = 1, OFF = 0
BAROCLINIC_RADIATION  : 1
!1 = Orlansky; 2 – Constant; 3 – Oey&Chen
CELERITY_TYPE         : 2
!celerity value assumed in the option  CELERITY_TYPE : 2
INTERNAL_CELERITY     : 2
```

### SedimentGeometry.dat

```
<begindomain>
ID                        : 1
TYPE                      : FIXSEDIMENT
LAYERS                    : 5
EMPTY_TOP_LAYERS          : 1
<<beginlayers>>
0.
0.25
0.25
0.25
0.25
<<endlayers>>
<enddomain>
```
