# Geostrophic Initialization

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

'**Warning: this page is under construction - do not use it directely'**

## Description

In this initialization methodology all forces are switched (except tidal force) at the first iteration. The velocity and sea surface height initial fields have non zero values consistent with the salinity and temperature initial fields, trough the geostrophic equilibrium. If the initialization is made in a level 2 embedded in a 2D tidal model (level 1), tidal force must be switched progressively along a selected period (in level 1).

## Example File

```
Hydrodynamic_1.dat

OUTPUT_TIME               : 0  3600
SURFACE_OUTPUT_TIME       : 0  3600

GEOST_INITIALIZATION      : 1   ! Select geostrophic initialization
SUBMODEL_FATHER_HOT_START : 0   ! must be zero for embedded models to prevent the imposition of father values

! imposed water level initial solution
<begin_waterlevel>
NAME                      : water level
DEFAULTVALUE              : 0
INITIALIZATION_METHOD     : HDF
FILENAME                  : ..\..\GeneralData\Mercator\Mercator_ANA_Level2_Jun.hdf5
TYPE_ZUV                  : z
ADDING_FACTOR             : 0.0  ! This factor adds to the mean level of the Initial Solution
<end_waterlevel>

TIDE                      : 0
TIDEPOTENTIAL             : 1    ! If model run embedded in a 2D tide model
CONTINUOUS                : 0    ! The model cannot be a continous run

DATA_ASSIMILATION         : 1
BRFORCE                   : 1
BRGEOMETRY                : 1
BRWATERLEVEL              : 0

RADIATION                 : 2
LOCAL_SOLUTION            : 3    ! 5 or 7 if model runs embedded
SUBMODEL                  : 0    ! 1 if model runs embedded
MISSING_NULL              : 1

WIND                      : 1
ATM_PRESSURE              : 1
ATM_PRESSURE_TYPE         : 2    ! can be 1 or 2
ATM_RAMP                  : 0
ATM_PERIOD                : 86400
IMPOSE_INVERTED_BAROMETER : 1

BIHARMONIC                : 1
BIHARMONIC_COEF           : 1e10
```
