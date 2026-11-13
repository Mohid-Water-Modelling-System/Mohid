# Evapotranspiration

Evapotranspiration is the sum of water evaporated from soil, canopy storage, surface water column and from plant transpiration.

## Reference Evapotranspiration

Reference evapotranspiration is an input for MOHID Land and can be provided by user or calculated. In any case, a property named *reference evapotranspiration* must exist in the Basin data file. This property can have a constant value in time, can be set to be calculated or data can be provided with a *Timeseries* file, a *HDF* file or any other valid input data file (see **[Module FillMatrix](../base/module-fillmatrix.md)**).

```
<beginproperty>
  NAME                               : reference evapotranspiration
  UNITS                              : mm/h
  DESCRIPTION                        : fao evapotranspiration
  DEFAULTVALUE                       : 0.2
  REMAIN_CONSTANT                    : 1
  NO_INTERPOLATION_OR_ACCUMULATION   : 1
<endproperty>
```

The above example sets a constant reference evapotranspiration of 0.2 mm/h during the entire simulation period.

To make MOHID Land calculate the reference evapotranspiration, the REMAIN\_CONSTANT keyword must be set to 0 (false) and no input data files should be provided. MOHID Land will then calculate the reference evapotranspiration using a set of properties that must be set in the Atmosphere data file (see **[Module Atmosphere](../water/modules/atmosphere.md)**). Below is an example of the reference evapotranspiration property set to be calculated:

```
<beginproperty>
  NAME                               : reference evapotranspiration
  UNITS                              : mm/h
  DESCRIPTION                        : fao evapotranspiration
  DEFAULTVALUE                       : 0.0
  REMAIN_CONSTANT                    : 0
  NO_INTERPOLATION_OR_ACCUMULATION   : 1
<endproperty>
```

Evapotranspiration is calculated using the 'FAO crop reference evapotranspiration' method, a standardized Penman-Monteith equation.

All necessary information can be obtained from the FAO website: [

$$
1
$$

](http://www.fao.org/docrep/X0490E/x0490e00.htm#Contents)

**ATTENTION**: The "UNITS" chosen for the reference evapotranspiration property are important. Internally MOHID Land will convert the values provided with the user units to m/s (meters by second). The units chosen for the property must match that of the input file if one is used. The units are always in LENGTH/TIME format.

```
The valid units for length are: m, cm, mm
The valid units for time are  : d (for days), h (for hours), s (for seconds)
```

## Crop Evapotranspiration

If **[Module Vegetation](../land/modules/vegetation.md)** is used, then the reference evapotranspiration will be adjusted using a crop coefficient (Kc) that must be provided in the Vegetation data file (*crop coefficient* property).

The Kc can be provided as a constant value or can came from one of the many input options at the disposal of the user through the **[Module FillMatrix](../base/module-fillmatrix.md)**, like as a timeserie file or as a HDF map, for example.

The crop evapotranspiration is calculated using the equation below:

$CropEvapotrans=RefEvapotrans \times K_{c}$

### Crop Coefficient Adjustment

Optionally, if the user has a timeseries or map of measured LAI and Potential LAI (LAI for maximum plant growth), this can be provided using the *leaf area index* and *potential leaf area index* properties in the Vegetation data file (see **[Module Vegetation](../land/modules/vegetation.md)**). When both properties are provided, the Kc can be "adjusted" using the LAI and Potential LAI. The equation used to adjust the Kc based on the Potential and Actual value of LAI can be seen below:

$K_{c}=K_{c}-[1.0-(\frac{LAI}{LAI_{dense}})^{0.5}]$

Where $LAI_{dense}$ is the potential LAI (*potential leaf area index* property see **[Module Vegetation](../land/modules/vegetation.md)**)

FOR THIS OPTION TO BE ACTIVE *potential leaf area index* MUST BE IN MODULE VEGETATION FILE

From 'Necessidades de Água e Métodos de Rega', Luis Santos Pereira, 2004 Publicações Europa-América, pg 85-86

### Crop Coefficient Adjustment related keywords

When using the LAI/Potential LAI to adjust the crop coefficient value, it's possible to set a minimum Kc and a value for Kc when LAI is ZERO. This is done using these keywords in the Basin data file:

```
DEFAULT_KC_WHEN_LAI_ZERO
KC_MIN
```

The **DEFAULT\_KC\_WHEN\_LAI\_ZERO** keyword tells that Kc must be set to this value when LAI is ZERO. The **KC\_MIN** keyword tells that if the Kc value found is lower than it, so the min value set by user will be used instead. It's a good idea to provide at least a **KC\_MIN** value (usually 0.3 for soil without vegetation).

## Potencial Evaporation & Potencial Transpiration

If the keyword **EVAPOTRANSPIRATION\_METHOD** in Basin data file was set to 2 (separete evaporation and transpiration), then the crop evapotranspiration will be separated into a Potential Evaporation and a Potencial Transpiration.

This is done using LAI (leaf area index), according with the equations below:

$PotentialTranspiration=CropEvapotrans \times (1.0-e^{-0.463 \times LAI})$

$PotentialEvaporation=CropEvapotrans-PotentialTranspiration$

## Actual Evapotranspiration

The actual (or real) evapotranspiration is calculated using the crop evapotranspiration (if **EVAPOTRANSPIRATION\_METHOD** = 1) or using the PotentialEvaporation and PotentialTranspiration (if **EVAPOTRANSPIRATION\_METHOD** = 2) and is subjected to the soil water content, stress factors and other limitations (like Head limits).

## Basin Output Timeserie

Depending on the value of the **EVAPOTRANSPIRATION\_METHOD** keyword, the 2D output timeserie file for basin (extension srb) will have a set or all of the columns listed next:

```
EvapoTranspiration_Rate_[mm/hour]   => Actual Evapotranspiration rate (actual transpiration + actual evaporation)
Potential_Crop_EVTP_[mm/hour]       => Crop Evapotranspiration rate   (reference evapotranspiration rate corrected using Kc)
Potential_Evaporation_[mm/h]        => Potential Evaporation rate     (only if EVAPOTRANSPIRATION_METHOD = 2)
Potential_Transpiration_[mm/h]      => Potential Transpiration rate   (only if EVAPOTRANSPIRATION_METHOD = 2)
Actual_Evaporation_[mm/h]           => Actual Evaporation rate        (only if EVAPOTRANSPIRATION_METHOD = 2)
Actual_Transpiration_[mm/h]         => Actual Transpiration rate      (only if EVAPOTRANSPIRATION_METHOD = 2)
Reference_Evapotranspiration_[mm/h] => Reference Evapotranspiration rate
```

## Links

- **[Module Basin](../land/modules/basin.md)**

- **[Module FillMatrix](../base/module-fillmatrix.md)**

- **[Module Atmosphere](../water/modules/atmosphere.md)**

- **[Module Vegetation](../land/modules/vegetation.md)**
