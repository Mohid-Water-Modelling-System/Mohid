# Module InterfaceWaterAir

## Overview

The water-air interface module is responsible by processes occurring at the water-air interface, such as computing [wind shear stress](../../reference/README.md#wind-surface-stress#wind-surface-stress), radiation balances, latent and sensible heat fluxes. This modules uses [Module Atmosphere](atmosphere.md) as a database for meteorological data and combines it with information from, for example, [Module Hydrodynamic](hydrodynamic.md) and [Module WaterProperties](waterproperties.md) to compute mass, heat and momentum fluxes across the water-air interface.

A user manual is available in the link at the bottom of the page. This user manual intends to help the user to couple atmosphere to water, activating in MOHID wind forcing, heat fluxes and mass fluxes between this interface.

## Main Processes

### Momentum fluxes

#### Surface rugosity

#### Wind shear stress

#### Wind shear velocity

#### Turbulent kinetic energy

### Heat fluxes

In MOHID heat fluxes in water column are computed with:

i) a heat source: solar radiation that enters through the water surface (surface radiation - see description) suffering a decay with depth (light extinction).

ii) a boundary condition for surface: non solar flux (see description)

iii) a boundary condition for bottom: no flux, all radiation reaching the bottom is transformed in heat

#### Short wave and long wave radiation

Solar Spectrum at the surface is composed from ultra-violet (UV), visible and infrared (IR) bands from 250nm to 2500nm (Ohlman and Siegel 2000b). From the total solar spectrum, around 250-400nm is UV, around 400-700nm is visible and from 700-2500nm is infrared (Monteith and Unsworth, 1990).

Earth emits radiation (IR radiation) in the 3000-100000nm (or 3-100um) (Monteith and Unsworth, 1990) which means that emitted radiation has different spectra than that coming from solar origin.

Atmosphere absorbs part of emitted radiation from earth (green house gases) and emitts in the same spectrum (from 3 to 100 um)- (Monteith and Unsworth, 1990).

It is common to find the term "short wave" associated with solar radiation and "long wave" with terrestrial and atmosphere radiation because of the different well defined spectral bands. In oceanography best care is taken to UV and visible bands because are the bands most sensible for fitoplankton and which penetrate in depth; infrared is rapidily attenuated in the first centimeters of water (Ohlman et al. 1996). As so, a "short wave solar radiation" term can be used integrating UV and visible bands and a "long wave solar radiation" for the correspondent IR band.

In MOHID it will be used the term "short wave" for UV and visible bands and longwave for IR band. This means that solar radiation has short wave (around 250-700nm) and long wave bands (around 700-2500nm) and that terrestrial (upward radiation) and atmosphere (downward radiation) are longwave bands as well (3-100um).

##### Short wave radiation

- Surface radiation

Surface radiation represents the solar radiation entering the water surface (after reflection - albedo). A fraction of solar radiation is short wave radiation (60% of surface radiation by default) .

The shortwave solar spectrum represents the UV and visible bands from solar radiation (around 250 - 700nm) which are sensible for biology and penetrate in depth in the water column. The fraction of this spectra to total solar radiation is variable with gases concentration (O3, O2, water vapour, CO2) in atmosphere, solar zenith, turbidity, and mainly, by cloudiness (Ohlman and Siegel 2000a). The 60% of surface radiation used in MOHID as default for solar short wave radiation corresponds to clear sky conditions. Accordingly to Ohlman and Siegel, 2000a UV and visible spectra can undergo to 70% and 80% of total solar radiation when 40% and 90% of radiation is reduced due to clouds, respectively.

Solar radiation is an atmosphere property (see module Atmosphere for details).

##### Long wave radiation

- Surface radiation

Surface radiation represents the solar radiation entering the water surface (after reflection - albedo). A fraction of solar radiation is long wave radiation (40% of surface radiation by default) .

The longwave solar spectrum represents the IR bands from solar radiation (around 700 - 2500nm) which are rapidilly attenuated in the first centimeters fo the water column. The fraction of this spectra to total solar radiation is variable with gases concentration (O3, O2, water vapour, CO2) in atmosphere, solar zenith, turbidity, and mainly, by cloudiness(Ohlman and Siegel 2000a). The 40% of surface radiation used in MOHID as default for solar long wave radiation corresponds to clear sky conditions. Adapted from Ohlman and Siegel, 2000a IR spectra from the solar radiation can undergo to 30% and 20% of total solar radiation when 40% and 90% of radiation is reduced due to clouds, respectively.

Solar radiation is an atmosphere property (see module Atmosphere for details).

- Upward long wave radiation

Heat emission from water to air related to its temperature. Represents the heat loss from water (to air) by radiation.

If computed, depends on water temperature.

- Downward long wave radiation

Heat emission from air to water related to its temperature. Represents the heat loss from air (to water) by radiation.

If computed, depends on air temperature and cloud cover.

- Net long wave radiation

Represents the balance between upward long wave radiation and downward long wave radiation. According to Monteith and Unsworth, 1990 net long wave radiation has exit direction from earth becuase not all upward radiation is absorbed and re-emitted by green house gases.

#### Latent heat

Absorbed or removed heat from water when changing phase at constant temperature. In other words, is the heat that affects the physical sate of water. Represents heat exchange between water and air in terms of evaporation and condensation.

If computed, depends on water temperature, air temperature, humidity and wind velocity.

#### Sensible heat

Absorbed or removed heat from water due to a change in temperature. No changes in physical state occur. In other words, is the heat that affects the temperature of water. Represents the heat transfer between water and air in terms of conduction and convection.

If computed, depends on water temperature, air temperature and wind velocity.

#### Non solar flux

The heat flux between water and air interface that has not solar origin. Represents the balance between latent heat, sensible heat and net long wave radiation.

Non solar flux is the boundary condition at the surface for water temperature calculation.

#### Using COARE algorithm

COARE3.0 algorithm for fluxes calculation can also be used. More information in [Coare](../../reference/coare.md).

### Mass fluxes

#### Oxygen

#### Carbon dioxide

#### Surface water fluxes

Mass flux between the water-air interface. Represents the balance between precipitation, evaporation and irrigation.

- Evaporation

Mass flux from water to air occurring at constant temperature. If computed, depends on latent heat, latent heat of vaporization (constant) and water density.

- Precipitation

Mass flux from rain.

- Irrigation

Mass flux from water removed by anthropogenic sources.

## User manual

[Coupling Water-Atmosphere User Manual](../coupling-water-atmosphere-user-manual.md)

## References

- Ohlman J.C., Siegel D.A. (2000a) - Ocean Radiant Heating. Part I: Optical Influences. Journal of Physical Oceanography, Volume 30 August 2000.
- Ohlman J.C., Siegel D.A. (2000b) - Ocean Radiant Heating. Part II: Parameterizing Solar Radiation Transmission through the Upper Ocean. Journal of Physical Oceanography, Volume 30 August 2000.
- Ohlman J.C., Siegel D.A., Gautier, C. (1996) - Ocean Mixed Layer Radiant heating and Solar Penetration: A Global Analysis. Journal of Climate, Volume 9, October 1996.
- Monteith, J.L., Unsworth, M.H. (1990) - Principles of Environmental Physics. Second Edition. Arnold Press, London. ISBN: 0 7131 2931 X

## Links

- [Module Atmosphere](atmosphere.md)
- [Coupling Water-Atmosphere Manual](../coupling-water-atmosphere-user-manual.md)
- [ConvertToHDF5](../../tools/converttohdf5.md)
- [Meteo](../../reference/meteo.md) - Available meteorological data
