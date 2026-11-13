# Module Atmosphere

## Overview

The atmosphere module is responsible for meteorological data needed to compute processes occurring at the water-air interface, such as computing wind shear stress, radiation balances, latent and sensible heat fluxes.

A user manual is available in the link at the bottom of the page. This user manual intends to help the user to couple atmosphere to water, activating in MOHID wind forcing, heat fluxes and mass fluxes between this interface.

## Concepts

## Main processes

### Air temperature

Air temperature is used for latent heat, sensible heat and downward longwave radiation computation (in the atmoshphere-water interface).

### Solar radiation

Heat flux with solar origin hitting the water surface. In absence of field data it can be computed based on two approaches that depend both on Top Of Atmosphere Radiation, sun height and cloud cover:

i) MOHID (default method)

ii) CE-QUAL based method.

Solar radiation is used to compute radiation that enters the water surface (after reflection). One method for cloud cover computation also uses radiation as input.

### Wind

Wind is described by its velocity decomposed on x and y components. If not available it can be computed from wind modulus and wind direction.

Wind is used to compute latent heat, sensible heat and wind stress (in the water-atmoshphere interface).

### Specific humidity

Specific humidity is the ratio of water mass present in the atmosphere against the air mass. Moisture mixture is the ratio of water mass against that of dry air. Relative humidity can be computed from specific humidity, pressure and air temperature when not given.

### Relative humidity

Relative humidity is a term used to describe the amount of water vapor that exists in a gaseous mixture of air and water. The relative humidity of an air-water mixture is defined as the ratio of the partial pressure of water vapor in the mixture to the saturated vapor pressure of water at a given temperature.

Relative humidity is used for latent heat computation (in the water-atmoshphere interface).

### Precipitation

Mass Flux from rain.

### Sun hours

### Cloud cover

Fraction of sky covered with clouds. Represents the radiation absorption in the atmosphere. It can be computed from three approaches:

i) comparison between measured sun hours (if available) with potential sun hours;

ii) random solution (this option was removed from June 2015, it is shown here for user help).

iii) comparison between solar radiation (measured or computed) and top of atmosphere radiation (default method);

Cloud cover is used for downward long wave radiation computation (in the atmoshphere-water interface). One method in solar radiation computation (if no data available) uses cloud cover as input.

**Remarks**:

- Cloud cover can not be computed from computed radiation because cloud cover will use computed radiation and radiation uses computed cloud cover ("pescadinha de rabo na boca" situation). This returns always zero cloud cover and zero radiation in the results (\*).

- Cloud cover computation from measured radiation has the advantage that is not a random calculation so it was used as default method after June 2015.

(\*)This happens because the model is initialized with zero radiation; then, if cloud cover is computed from radiation it returns zero cloud cover in next time step which in turn results in zero transmissivity and zero radiation, going on cycle.

### Irrigation

Mass flux from water added by anthropogenic sources.

### Other Properties

When using **[Module PorousMediaProperties](../../land/modules/porousmediaproperties.md)** or **Module RunOffProperties**, it's possible to include property concentrations in the Atmosphere data file, that will be used as "top boundary condition" for these properties.

Ex.:

```
<beginproperty>
  NAME                               : solution ammonia
  UNITS                              : mg/L
  DESCRIPTION                        : Ammonia in precipitation
  PRECIPITATION                      : 1
  FILE_IN_TIME                       : TIMESERIE
  FILENAME                           : ..\GeneralData\Boundary Conditions\ammonia_in_precipitation_timeseries.dat
  DATA_COLUMN                        : 2
  NO_INTERPOLATION_OR_ACCUMULATION   : 1
  DEFAULTVALUE                       : 0
  REMAIN_CONSTANT                    : 0
  TIME_SERIE                         : 1
<endproperty>
```

## Other features

### Random components

## User Manual

[Coupling Water-Atmosphere User Manual](../coupling-water-atmosphere-user-manual.md)

## Links

- **[Module\_InterfaceWaterAir](interfacewaterair.md)**
- **[ConvertToHDF5](../../tools/converttohdf5.md)**
- **[FillMatrix](../../tools/fillmatrix.md)**
