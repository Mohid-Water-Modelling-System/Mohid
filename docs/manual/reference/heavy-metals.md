# Heavy metals

## Introduction

In MOHID the distribution of contaminants is parameterized by using the ratio between the adsorbed particulate concentration and the dissolved concentration (see [Adsorption/Desorption](../water/adsorption-desorption.md) ). Physically, the partition coefficient (Johansson et al., 2001) illustrates particle affinity and represents the chemical equilibrium of numerous processes such adsorption onto particulate matter, precipitation and dissolution. Pollutants transport in the adsorbed phase is strongly connected with cohesive sediment transport due to the affinity of many contaminants for the solid phase. Once adsorbed to the sediments, these substances are transported by the sediments (Fernandes L., 2005).

In MOHID the partition coefficient is defined through a keyword (keyword PARTITION\_FRACTION). The concentration of the contaminant in the water is described by the following equations (1) and (2):

```
            (1)
```

```
     (2)
```

where:

Cd is the concentration of the contaminant in the dissolved phase (mol/l)

Cp is the concentration of the contaminant in the particulate phase (mol/l)

r is the transfer rate from particulate to dissolved phase (1/sec) (keyword PARTITION\_RATE)

pd is the dissolved fraction (dimensionless)

ps is the solid fraction (dimensionless)

For each contaminant phase it is necessary to specify opposite phase through the keyword PARTITION\_COUPLE. For example, for the contaminant dissolved arsenic, the opposite phase is particulate arsenic. It is assumed that particulate and dissolved phases must have equal partition rates. Also, the sum of particulate and dissolved phases’ fractions must be equal to 1.

In MOHID, the transfer rate depends on the cohesive sediment concentration. A coefficient Rs is defined as (eq.3):

```
          (3)
```

SedimentRefConcentration is a cohesive sediment reference concentration (keyword SED\_REF\_CONC). Particulate and dissolved phases of the same contaminant must have equal cohesive sediment reference concentration. The transfer rate is equal to r for Rs > 1 and equal to r x Rs for Rs<1

## Eulerian Applications

### Example input files

#### File WaterProperty\_X.dat

```
OUTPUT_TIME               : 0.  3600.

<beginproperty>
 NAME                      : particulate copper
 UNITS                     : mg/l
 DESCRIPTION               : No description was given.
 PARTICULATE               : 1
 BOUNDARY_INITIALIZATION   : EXTERIOR
 DEFAULTBOUNDARY           : 0.0
 BOUNDARY_CONDITION        : 2
 ADVECTION_DIFFUSION       : 1
 ADVECTION_H_IMP_EXP       : 0
 INITIALIZATION_METHOD     : CONSTANT
 DEFAULTVALUE              : 0.
 IS_COEF                   : 0.001
 PARTITION                 : 1
 DISCHARGES                : 1
 VERTICAL_MOVEMENT         : 1
 BOTTOM_FLUXES             : 1
 OUTPUT_HDF                : 1
 TIME_SERIE                : 0
 BOX_TIME_SERIE            : 0
 PARTITION                 : 1
 PARTITION_FRACTION        : 0.7
 PARTITION_RATE            : 5e-5
 PARTITION_COUPLE          : dissolved copper
 USE_SED_REF_CONC          : 1
 SED_REF_CONC              : 2.
<endproperty>

<beginproperty>
 NAME                      : dissolved copper
 UNITS                     : mg/l
 DESCRIPTION               : No description was given.
 PARTICULATE               : 0
 BOUNDARY_INITIALIZATION   : EXTERIOR
 DEFAULTBOUNDARY           : 0.02
 BOUNDARY_CONDITION        : 2
 ADVECTION_DIFFUSION       : 1
 ADVECTION_H_IMP_EXP       : 0
 INITIALIZATION_METHOD     : CONSTANT
 DEFAULTVALUE              : 0.02
 IS_COEF                   : 0.001
 PARTITION                 : 1
 DISCHARGES                : 1
 VERTICAL_MOVEMENT         : 0
 BOTTOM_FLUXES             : 0
 OUTPUT_HDF                : 1
 TIME_SERIE                : 0
 BOX_TIME_SERIE            : 0
 PARTITION                 : 1
 PARTITION_FRACTION        : 0.3
 PARTITION_RATE            : 5e-5
 PARTITION_COUPLE          : particulate copper
 USE_SED_REF_CONC          : 1
 SED_REF_CONC              : 2.
<endproperty>

<beginproperty>
 NAME                      : particulate zinc
 UNITS                     : mg/l
 DESCRIPTION               : No description was given.
 PARTICULATE               : 1
 BOUNDARY_INITIALIZATION   : EXTERIOR
 DEFAULTBOUNDARY           : 0.0
 BOUNDARY_CONDITION        : 2
 ADVECTION_DIFFUSION       : 1
 ADVECTION_H_IMP_EXP       : 0
 INITIALIZATION_METHOD     : CONSTANT
 DEFAULTVALUE              : 0.
 IS_COEF                   : 0.001
 PARTITION                 : 1
 DISCHARGES                : 1
 VERTICAL_MOVEMENT         : 1
 BOTTOM_FLUXES             : 1
 OUTPUT_HDF                : 1
 TIME_SERIE                : 0
 BOX_TIME_SERIE            : 0
 PARTITION                 : 1
 PARTITION_FRACTION        : 0.7
 PARTITION_RATE            : 5e-5
 PARTITION_COUPLE          : dissolved zinc
 USE_SED_REF_CONC          : 1
 SED_REF_CONC              : 2.
<endproperty>

<beginproperty>
 NAME                      : dissolved zinc
 UNITS                     : mg/l
 DESCRIPTION               : No description was given.
 PARTICULATE               : 0
 BOUNDARY_INITIALIZATION   : EXTERIOR
 DEFAULTBOUNDARY           : 0.02
 BOUNDARY_CONDITION        : 2
 ADVECTION_DIFFUSION       : 1
 ADVECTION_H_IMP_EXP       : 0
 INITIALIZATION_METHOD     : CONSTANT
 DEFAULTVALUE              : 0.02
 IS_COEF                   : 0.001
 PARTITION                 : 1
 DISCHARGES                : 1
 VERTICAL_MOVEMENT         : 0
 BOTTOM_FLUXES             : 0
 OUTPUT_HDF                : 1
 TIME_SERIE                : 0
 BOX_TIME_SERIE            : 0
 PARTITION                 : 1
 PARTITION_FRACTION        : 0.3
 PARTITION_RATE            : 5e-5
 PARTITION_COUPLE          : particulate zinc
 USE_SED_REF_CONC          : 1
 SED_REF_CONC              : 2.
<endproperty>

<beginproperty>
 NAME                      : particulate metal
 UNITS                     : mg/l
 DESCRIPTION               : No description was given.
 PARTICULATE               : 1
 BOUNDARY_INITIALIZATION   : EXTERIOR
 DEFAULTBOUNDARY           : 0.0
 BOUNDARY_CONDITION        : 2
 ADVECTION_DIFFUSION       : 1
 ADVECTION_H_IMP_EXP       : 0
 INITIALIZATION_METHOD     : CONSTANT
 DEFAULTVALUE              : 0.
 IS_COEF                   : 0.001
 PARTITION                 : 1
 DISCHARGES                : 1
 VERTICAL_MOVEMENT         : 1
 BOTTOM_FLUXES             : 1
 OUTPUT_HDF                : 1
 TIME_SERIE                : 0
 BOX_TIME_SERIE            : 0
 PARTITION                 : 1
 PARTITION_FRACTION        : 0.7
 PARTITION_RATE            : 5e-5
 PARTITION_COUPLE          : dissolved metal
 USE_SED_REF_CONC          : 1
 SED_REF_CONC              : 2.
<endproperty>

<beginproperty>
 NAME                      : dissolved metal
 UNITS                     : mg/l
 DESCRIPTION               : No description was given.
 PARTICULATE               : 0
 BOUNDARY_INITIALIZATION   : EXTERIOR
 DEFAULTBOUNDARY           : 0.02
 BOUNDARY_CONDITION        : 2
 ADVECTION_DIFFUSION       : 1
 ADVECTION_H_IMP_EXP       : 0
 INITIALIZATION_METHOD     : CONSTANT
 DEFAULTVALUE              : 0.02
 IS_COEF                   : 0.001
 PARTITION                 : 1
 DISCHARGES                : 1
 VERTICAL_MOVEMENT         : 0
 BOTTOM_FLUXES             : 0
 OUTPUT_HDF                : 1
 TIME_SERIE                : 0
 BOX_TIME_SERIE            : 0
 PARTITION                 : 1
 PARTITION_FRACTION        : 0.3
 PARTITION_RATE            : 5e-5
 PARTITION_COUPLE          : particulate metal
 USE_SED_REF_CONC          : 1
 SED_REF_CONC              : 2.
<endproperty>
```

#### File InterfaceSedimentWater\_X.dat

```
OUTPUT_TIME                  : 0  3600

<begin_rugosity>
 INITIALIZATION_METHOD        : CONSTANT
 DEFAULTVALUE                 : 0.0025
<end_rugosity>

<begin_critical_shear_erosion>
 INITIALIZATION_METHOD        : CONSTANT
 DEFAULTVALUE                 : 0.2
<end_critical_shear_erosion>

<begin_critical_shear_deposition>
 INITIALIZATION_METHOD        : CONSTANT
 DEFAULTVALUE                 : 0.1
<end_critical_shear_deposition>

<begin_erosion_rate>
 INITIALIZATION_METHOD        : CONSTANT
 DEFAULTVALUE                 : 0.00005
<end_erosion_rate>

<beginproperty>
 NAME                         : cohesive sediment
 UNITS                        : kg/m2
 DESCRIPTION                  : cohesive sediment mass
 PARTICULATE                  : 1
 WATER_FLUXES                 : 1
 EROSION                      : 1
 DEPOSITION                   : 1
 MASS_LIMITATION              : 1
 MASS_MIN                     : 0.
 INITIALIZATION_METHOD        : CONSTANT
 DEFAULTVALUE                 : 0.00001
 TIME_SERIE                   : 0
 BOX_TIME_SERIE               : 0
 OUTPUT_HDF                   : 1
<endproperty>

<beginproperty>
 NAME                         : particulate copper
 UNITS                        : kg/m2
 DESCRIPTION                  : PON mass
 PARTICULATE                  : 1
 WATER_FLUXES                 : 1
 EROSION                      : 1
 DEPOSITION                   : 1
 MASS_LIMITATION              : 1
 MASS_MIN                     : 0.
 INITIALIZATION_METHOD        : CONSTANT
 DEFAULTVALUE                 : 0.
 TIME_SERIE                   : 0
 BOX_TIME_SERIE               : 0
 OUTPUT_HDF                   : 1
<endproperty>

<beginproperty>
 NAME                         : particulate zinc
 UNITS                        : kg/m2
 DESCRIPTION                  : PON mass
 PARTICULATE                  : 1
 WATER_FLUXES                 : 1
 EROSION                      : 1
 DEPOSITION                   : 1
 MASS_LIMITATION              : 1
 MASS_MIN                     : 0.
 INITIALIZATION_METHOD        : CONSTANT
 DEFAULTVALUE                 : 0.
 TIME_SERIE                   : 0
 BOX_TIME_SERIE               : 0
 OUTPUT_HDF                   : 1
<endproperty>

<beginproperty>
 NAME                         : particulate metal
 UNITS                        : kg/m2
 DESCRIPTION                  : PON mass
 PARTICULATE                  : 1
 WATER_FLUXES                 : 1
 EROSION                      : 1
 DEPOSITION                   : 1
 MASS_LIMITATION              : 1
 MASS_MIN                     : 0.
 INITIALIZATION_METHOD        : CONSTANT
 DEFAULTVALUE                 : 0.
 TIME_SERIE                   : 0
 BOX_TIME_SERIE               : 0
 OUTPUT_HDF                   : 1
<endproperty>
```

#### File FreeVerticalMovement\_X.dat

```
<beginproperty>
 NAME                     : cohesive sediment
 WS_VALUE                 : 0.00001
 WS_TYPE                  : 1
 DEPOSITION               : 1
<endproperty>

<beginproperty>
 NAME                     : particulate copper
 WS_VALUE                 : 0.00001
 WS_TYPE                  : 1
 DEPOSITION               : 1
<endproperty>

<beginproperty>
 NAME                     : particulate zinc
 WS_VALUE                 : 0.00001
 WS_TYPE                  : 1
 DEPOSITION               : 1
<endproperty>

<beginproperty>
 NAME                     : particulate metal
 WS_VALUE                 : 0.00001
 WS_TYPE                  : 1
 DEPOSITION               : 1
<endproperty>
```

#### File Discharges\_X.dat

```
<begindischarge>
 NAME                          : Discharge
 I_CELL                        : 28
 J_CELL                        : 8
 K_CELL                        : 1
 ALTERNATIVE_LOCATIONS         : 1
 DEFAULT_FLOW_VALUE            : 50

<<beginproperty>>
 NAME                          : cohesive sediment
 UNITS                         : psu
 DEFAULTVALUE                  : 50
<<endproperty>>

<<beginproperty>>
 NAME                          : particulate copper
 UNITS                         : mg/l
 DEFAULTVALUE                  : 0.02
<<endproperty>>

<<beginproperty>>
 NAME                          : dissolved copper
 UNITS                         : mg/l
 DEFAULTVALUE                  : 0.02
<<endproperty>>

<<beginproperty>>
 NAME                          : particulate zinc
 UNITS                         : mg/l
 DEFAULTVALUE                  : 0.02
<<endproperty>>

<<beginproperty>>
 NAME                          : dissolved zinc
 UNITS                         : mg/l
 DEFAULTVALUE                  : 0.02
<<endproperty>>

<<beginproperty>>
 NAME                          : particulate metal
 UNITS                         : mg/l
 DEFAULTVALUE                  : 0.02
<<endproperty>>

<<beginproperty>>
 NAME                          : dissolved metal
 UNITS                         : mg/l
 DEFAULTVALUE                  : 0.02
<<endproperty>>

<enddischarge>
```

## Lagrangian Applications

### Example input files

#### File Model\_X.dat

```
START                         : 2008 6 1 12 0 0
END                           : 2008 6 5 12 0 0
DT                            : 30.
VARIABLEDT                    : 0
GMTREFERENCE                  : 0
SPLITTING                     : Double_Splitting
LAGRANGEANE                   : 1
LAGRANGIAN                    : 1
```

#### File Lagrangian\_X.dat

```
OUTPUT_TIME             : 0 3600.
DT_PARTIC               : 30.
OUTPUT_CONC             : 2

<BeginOrigin>
 ORIGIN_NAME             : Emission (60um)
 OLD                     : 0
 GROUP_ID                : 2
 EMISSION_SPATIAL        : Point
 EMISSION_TEMPORAL       : Continuous
 DT_EMIT                 : 150
 FLOW                    : 1
 FLOAT                   : 0
 MOVEMENT                : SullivanAllen
 VARVELHX                : 0.0
 VARVELH                 : 0.03
 TURB_V                  : Profile
 POSITION_CELLS          : 7.5 27.5
 DEPTH_METERS            : 0.
 NBR_PARTIC              : 1
 SEDIMENTATION           : Stokes
 D50                     : 0.06
 MIN_SED_VELOCITY        : 0.0
 DEPOSITION              : 1
 TAU_ERO                 : 0.2
 TAU_DEP                 : 0.1
 BOTTOM_DISTANCE         : 0.1
 BOTTOM_EMISSION         : 0
 TIME_DECAY              : 172800.
 EROSION_RATE            : 0.05

<<BeginProperty>>
 NAME                    : sediment
 UNITS                   : mg/l
 CONCENTRATION           : 10.
 AMBIENT_CONC            : 1
<<EndProperty>>

<<BeginProperty>>
 NAME                    : particulate zinc
 UNITS                   : mg/l
 CONCENTRATION           : 0.02
 AMBIENT_CONC            : 0
 PARTITION_WATER         : 1
 PARTITION_COEF_WATER    : 0.7
 PARTITION_RATE_WATER    : 5e-5
 PARTITION_COUPLE_WATER  : 0.02
 PARTITION_SED           : 1
 PARTITION_COEF_SED      : 0.9
 PARTITION_RATE_SED      : 5e-5
 PARTITION_COUPLE_SED    : 0.02
<<EndProperty>>

<EndOrigin>

<BeginOrigin>
 ORIGIN_NAME             : Emission (40um)
 OLD                     : 0
 GROUP_ID                : 1
 EMISSION_SPATIAL        : Point
 EMISSION_TEMPORAL       : Continuous
 DT_EMIT                 : 150
 FLOW                    : 1
 FLOAT                   : 0
 MOVEMENT                : SullivanAllen
 VARVELHX                : 0.0
 VARVELH                 : 0.03
 TURB_V                  : Profile
 POSITION_CELLS          : 7.5 27.5
 DEPTH_METERS            : 0.
 NBR_PARTIC              : 1
 SEDIMENTATION           : Stokes
 D50                     : 0.04
 MIN_SED_VELOCITY        : 0.0
 DEPOSITION              : 1
 TAU_ERO                 : 0.2
 TAU_DEP                 : 0.1
 BOTTOM_DISTANCE         : 0.1
 BOTTOM_EMISSION         : 0
 TIME_DECAY              : 172800.
 EROSION_RATE            : 0.05

<<BeginProperty>>
 NAME                    : sediment
 UNITS                   : mg/l
 CONCENTRATION           : 10.
 AMBIENT_CONC            : 1
<<EndProperty>>

<<BeginProperty>>
 NAME                    : particulate copper
 UNITS                   : mg/l
 CONCENTRATION           : 0.02
 AMBIENT_CONC            : 0
 PARTITION_WATER         : 1
 PARTITION_COEF_WATER    : 0.7
 PARTITION_RATE_WATER    : 5e-5
 PARTITION_COUPLE_WATER  : 0.02
 PARTITION_SED           : 1
 PARTITION_COEF_SED      : 0.9
 PARTITION_RATE_SED      : 5e-5
 PARTITION_COUPLE_SED    : 0.02
<<EndProperty>>

<EndOrigin>

<BeginOrigin>
 ORIGIN_NAME             : Emission (20um)
 OLD                     : 0
 GROUP_ID                : 1
 EMISSION_SPATIAL        : Point
 EMISSION_TEMPORAL       : Continuous
 DT_EMIT                 : 150
 FLOW                    : 1
 FLOAT                   : 0
 MOVEMENT                : SullivanAllen
 VARVELHX                : 0.0
 VARVELH                 : 0.03
 TURB_V                  : Profile
 POSITION_CELLS          : 7.5 27.5
 DEPTH_METERS            : 0.
 NBR_PARTIC              : 1
 SEDIMENTATION           : Stokes
 D50                     : 0.02
 MIN_SED_VELOCITY        : 0.0
 DEPOSITION              : 1
 TAU_ERO                 : 0.2
 TAU_DEP                 : 0.1
 BOTTOM_DISTANCE         : 0.1
 BOTTOM_EMISSION         : 0
 TIME_DECAY              : 172800.
 EROSION_RATE            : 0.05

<<BeginProperty>>
 NAME                    : sediment
 UNITS                   : mg/l
 CONCENTRATION           : 10.
 AMBIENT_CONC            : 1
<<EndProperty>>

<<BeginProperty>>
 NAME                    : particulate copper
 UNITS                   : mg/l
 CONCENTRATION           : 0.02
 AMBIENT_CONC            : 0
 PARTITION_WATER         : 1
 PARTITION_COEF_WATER    : 0.7
 PARTITION_RATE_WATER    : 5e-5
 PARTITION_COUPLE_WATER  : 0.02
 PARTITION_SED           : 1
 PARTITION_COEF_SED      : 0.9
 PARTITION_RATE_SED      : 5e-5
 PARTITION_COUPLE_SED    : 0.02
<<EndProperty>>

<EndOrigin>

<BeginOrigin>
 ORIGIN_NAME             : Emission (10um)
 OLD                     : 0
 GROUP_ID                : 1
 EMISSION_SPATIAL        : Point
 EMISSION_TEMPORAL       : Continuous
 DT_EMIT                 : 150
 FLOW                    : 1
 FLOAT                   : 0
 MOVEMENT                : SullivanAllen
 VARVELHX                : 0.0
 VARVELH                 : 0.03
 TURB_V                  : Profile
 POSITION_CELLS          : 7.5 27.5
 DEPTH_METERS            : 0.
 NBR_PARTIC              : 1
 SEDIMENTATION           : Stokes
 D50                     : 0.01
 MIN_SED_VELOCITY        : 0.0
 DEPOSITION              : 1
 TAU_ERO                 : 0.2
 TAU_DEP                 : 0.1
 BOTTOM_DISTANCE         : 0.1
 BOTTOM_EMISSION         : 0
 TIME_DECAY              : 172800.
 EROSION_RATE            : 0.05

<<BeginProperty>>
 NAME                    : sediment
 UNITS                   : mg/l
 CONCENTRATION           : 10.
 AMBIENT_CONC            : 1
<<EndProperty>>

<<BeginProperty>>
 NAME                    : particulate copper
 UNITS                   : mg/l
 CONCENTRATION           : 0.02
 AMBIENT_CONC            : 0
 PARTITION_WATER         : 1
 PARTITION_COEF_WATER    : 0.7
 PARTITION_RATE_WATER    : 5e-5
 PARTITION_COUPLE_WATER  : 0.02
 PARTITION_SED           : 1
 PARTITION_COEF_SED      : 0.9
 PARTITION_RATE_SED      : 5e-5
 PARTITION_COUPLE_SED    : 0.02
<<EndProperty>>

<EndOrigin>

<BeginOrigin>
 ORIGIN_NAME             : Emission (2um)
 OLD                     : 0
 GROUP_ID                : 1
 EMISSION_SPATIAL        : Point
 EMISSION_TEMPORAL       : Continuous
 DT_EMIT                 : 150
 FLOW                    : 1
 FLOAT                   : 0
 MOVEMENT                : SullivanAllen
 VARVELHX                : 0.0
 VARVELH                 : 0.03
 TURB_V                  : Profile
 POSITION_CELLS          : 7.5 27.5
 DEPTH_METERS            : 0.
 NBR_PARTIC              : 1
 SEDIMENTATION           : Stokes
 D50                     : 0.02
 MIN_SED_VELOCITY        : 0.0
 DEPOSITION              : 1
 TAU_ERO                 : 0.2
 TAU_DEP                 : 0.1
 BOTTOM_DISTANCE         : 0.1
 BOTTOM_EMISSION         : 0
 TIME_DECAY              : 172800.
 EROSION_RATE            : 0.05

<<BeginProperty>>
 NAME                    : sediment
 UNITS                   : mg/l
 CONCENTRATION           : 10.
 AMBIENT_CONC            : 1
<<EndProperty>>

<<BeginProperty>>
 NAME                    : particulate copper
 UNITS                   : mg/l
 CONCENTRATION           : 0.02
 AMBIENT_CONC            : 0
 PARTITION_WATER         : 1
 PARTITION_COEF_WATER    : 0.7
 PARTITION_RATE_WATER    : 5e-5
 PARTITION_COUPLE_WATER  : 0.02
 PARTITION_SED           : 1
 PARTITION_COEF_SED      : 0.9
 PARTITION_RATE_SED      : 5e-5
 PARTITION_COUPLE_SED    : 0.02
<<EndProperty>>

<EndOrigin>
```

## References

Johasson, H., Lindstrom, M., Hakanson, L., (2001). On the Modelling of particulate and dissolved fractions of substances in aquatic systems – sedimentological and ecological interactions, Ecological Modelling, 137, 225-240

Fernandes L., (2005). Modelling of arsenic Dynamics in the Tagus Estuary. MSc dissertation thesis. Technical University of Lisbon.
