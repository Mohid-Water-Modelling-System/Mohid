# AssimilationZones

The **AssimilationZones** is an application used for generating **[assimilation](../water/modules/assimilation.md) decay coefficients** (also called **relaxation times**) files, usable in MOHID Water to perform relaxation of hydrodynamic and water properties fields to a reference field (e.g. a father model domain solution, a climatological field).

Running options for this application are specified by the user in a input file named [AssimilationZones.dat](assimilationzones.md).

## Introduction

**Typical use:**

3D decay coefficients files are desirable to perform **temperature assimilation** on MOHID Water at the surface layer based on climatology fields, emulating the consideration of heat fluxes in the water-air interface which are in equilibrium relatively to the climatological fields.

This application can be used for less orthodox ends such as in the design of **viscosity sponges** files to attenuate the gradient of velocity (aim is to reduce the reflection of perturbations to the inner domain due to boundary).

**Data input requirements:**

**Ouput:**

## Input file (AssimilationZones.dat)

```
IN_BATIM                 : ... (path/name of grid data bathymetry file of domain to which
                                decay coefficients are to be referred)
OUT_FILE                 : ... (path/name of coefficients/grid data files to be created)

SURROUND                 : 0/1 (1(default)=define assimilation zone in grid outer limits)
(if SURROUND : 1:)
NUMBER_CELLS             : ... (integer number of assimilation cells in XX direction)
ASSIMILATION_COEFS       : ... (decay coefficients (outside-inside) in seconds, one per column)

SURFACE_FORCING          : 0/1 (1=assimilate values in the surface layer; 0 = default)
(if SURFACE_FORCING : 1:)
SURFACE_COEF             : ... (decay coefficient (seconds) for the surface layer)
LAYERS_NUMBER            : ... (number of layers for 3D decay coefficients)

TYPE_Z                   : 0/1 (1(default)=output refers to cells center)
TYPE_U                   : 0/1 (1=output refers to cells faces in XX direction; 0 = default)
TYPE_V                   : 0/1 (1=output refers to cells faces in YY direction; 0 = default)

WEST                     : 0/1 (1(default)=decay values in the west boundary)
EAST                     : 0/1 (1(default)=decay values in the east boundary)
SOUTH                    : 0/1 (1(default)=decay values in the south boundary)
NORTH                    : 0/1 (1(default)=decay values in the north boundary)

DEFAULTVALUE             : ... (default value in seconds assumed for decay coefficients;
                                1e32 s. = default)

MINIMUM_DOMINATES        : 0/1 (1(default)=lower decay coefficients dominate over the higher
                           values)

NO_ASSIMILATION_ZONE     : 0/1 (1=use a polygon file to define where assimilation is not to be
                                defined)
(if NO_ASSIMILATION_ZONE : 1:)
NO_ASSIMILATION_FILE     : ... (path/name to polygon file defining no assimilation zone)

WRITE_AS_GRID_DATA       : 0/1 (1(default)=write grid data regular format; only if TYPE_Z : 1;
                                this grid data file can be viewed in MOHID GIS to check
                                coefficient values)
WRITE_IN_INDICES_FORMAT  : 0/1 (1(default)=write as [i],[j],[value] format)
```

```
Remarks:

- for any coefficients files to be generated either SURROUND : 1 or SURFACE_FORCING : 1
  options must be chosen;

- in case of SURROUND : 1 option the keywords NUMBER_CELLS and ASSIMILATION_COEFS must be
  specified;

- if SURFACE_FORCING : 1 then 3D decay coefficients files are produced with the number of
  vertical layers specified in LAYERS_NUMBER; in practice all layers have the same values for
  decay coefficients except for the surface layer which takes into account the info provided
  on SURFACE_COEF and MINIMUM_DOMINATES;

- for 3D decay coefficients files to be produced should always be specified SURFACE_FORCING :
  1; if no special coefficients are intended for the surface layer then convenient
  instructions must be provided through SURFACE_COEF and MINIMUM_DOMINATES keywords (see
  previous remark);

- the value specified in DEFAULTVALUE is used for decay coefficients in cells not specified
  by the SURROUND or SURFACE_FORCING keywords.
```

## Sample

```
IN_BATIM                 : Outfall.dat
OUT_FILE                 : OutfallCoefs.dat

SURROUND                 : 1
NUMBER_CELLS             : 11
ASSIMILATION_COEFS       : 3e4 3e4 1e5 3e5 1e6 3e6 1e7 3e7 1e8 3e8 1e9

TYPE_Z                   : 1
TYPE_U                   : 0
TYPE_V                   : 0

NO_ASSIMILATION_ZONE     : 0

WRITE_AS_GRID_DATA       : 1
WRITE_IN_INDICES_FORMAT  : 1

DEFAULTVALUE             : 1e9

SURFACE_FORCING          : 1
SURFACE_COEF             : 1e6

LAYERS_NUMBER            : 38
```
