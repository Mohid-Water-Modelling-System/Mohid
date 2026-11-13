# Module Basin

## Overview

Module Basin works as an interface among the different modules of Mohid-Land. Indeed it manages fluxes between modules as precipitation, evapotranspiration, infiltration, etc and updates water column and concentration after each module call. This module is able to compute a water and mass balance for each property transported in all mediums.

## Main Processes

The processes made in the Module Basin can be summarized as following:

- Reading entering data and grid construction
- Atmospheric processes (precipitation, leaf interception, leaf drainage, evaporation) in order to obtain the potential water column
- Call of [Module PorousMedia](porousmedia.md) giving potential water column and obtain the infiltration rate
- Update of the water column and send it to ModuleRunoff (the holder of water column)
- Call of [Module PorousMediaProperties](porousmediaproperties.md) and update of water column concentrations send it to the ModuleRunoffProperties
- Call of [Module Runoff](runoff.md) giving the remaining water columns to be transported
- Call of [Module RunoffProperties](runoffproperties.md)

(When Module Runoff and RunoffProperties run as they are the holders of water column and water column concentration, no update is needed).

- Call of [Module DrainageNetwork](../../water/modules/drainagenetwork.md) to route the water in the river and the new transfered from groundwater and from runoff.
- Output of the different components of the water and property flux

## Evapotranspiration

Some water may be extracted from the soil because of the evaporation and transpiration processes, which become a sink in soil water profile. These two processes are currently named Evapotranspiration and Potential Evapotranspiration may be modeled using the Penmann Monteith equation:

$\overset{\text{Energy flux rate}}{\lambda_v E=\frac{\Delta R_n   +   \rho_a c_p  \left(  \delta q  \right) g_a }{\Delta  + \gamma \left (    1 + g_a / g_s    \right)}}\Leftrightarrow$ $\overset{\text{Volume flux rate}}{ET_o=\frac{\Delta R_n   +   \rho_a c_p  \left(  \delta q  \right) g_a } { \left(   \Delta  + \gamma \left (    1 + g_a / g_s    \right)    \right) \lambda_v }}$

*λ*~v~ = Latent heat of vaporization. Energy required per unit mass of water vaporized. (J/g)

*L*~v~ = Volumetric latent heat of vaporization. Energy required per water volume vaporized. (*L*~v~ = 2453 MJ m^\-3^)

*E* = Mass water evapotranspiration rate (g s^\-1^ m^\-2^)

*ET*~o~ = Water volume evapotranspired (m^3^ s^\-1^ m^\-2^)

Δ = Rate of change of saturation specific humidity with air temperature. (Pa K^\-1^)

*R*~n~ = Net irradiance (W m^\-2^), the external source of energy flux

*c*~p~ = Specific heat capacity of air (J kg^\-1^ K^\-1^)

*ρ*~a~ = dry air density (kg m^\-3^)

δ*e* = vapor pressure deficit, or specific humidity (Pa)

*g*~a~ = Hydraulic conductivity of air, atmospheric conductance (m s^\-1^)

*g*~s~ = Conductivity of stoma, surface conductance (m s^\-1^)

*γ* = Psychrometric constant (*γ* ≈ 66 Pa K^\-1^)

Mohid Land has a sub-daily time step. According with FAO-56 [Hourly Time Step ET](http://www.fao.org/docrep/X0490E/x0490e08.htm) is calulated with the equation below. This is the equation used in Mohid-Land.

$\ ET_o=\frac{ \left  0.408\Delta(R_n - G)   +   \gamma \frac{37}{T_h+273}u_2 \left ( e^{o}*T_h-e_a \right )} {    \Delta  + \gamma \left (1 + 0.34u_2) }$

ETo = reference evapotranspiration 

$$
mm hour-1
$$



Rn = net radiation at the grass surface 

$$
MJ m-2 hour-1
$$



G = soil heat flux density 

$$
MJ m-2 hour-1
$$



T~h~ = mean hourly air temperature at 2m hight 

$$
°C
$$



Δ = saturation slope vapour pressure curve at Th 

$$
kPa °C-1
$$



*γ* = psychrometric constant 

$$
kPa °C-1
$$



e°(Th) = saturation vapour pressure at air temperature Th 

$$
kPa
$$



e~a~ = average hourly actual vapour pressure 

$$
kPa
$$



u~2~ = average hourly wind speed at 2m hight 

$$
m s-1
$$

.

Calculate Psicrometric constant:

$\gamma =\frac{ \left( c_p \right)_{air} *  P }{ \lambda_v * MW_{ratio} }$

$\gamma =$ psychrometric constant 

$$
kPa °C^\-1^
$$

,

P = atmospheric pressure 

$$
kPa
$$

,

$\lambda_v =$ latent heat of water vaporization, 2.45 

$$
MJ kg^\-1^
$$

,

$c_p =$ specific heat of air at constant pressure, 1.013 10-3 

$$
MJ kg^\-1^ °C^\-1^
$$

,

$MW_{ratio} =$ ratio molecular weight of water vapor/dry air = 0.622.

Equation of $\gamma$ becomes:

$\gamma =  0.665*10^{-3} * P$

Calculation of the atmospheric pressure based on the heigth simplification of the ideal gas law

$P = 101.3  \left ( \frac{293.-0.0065 * Elevation}{293.}\right  ) ^{5.26}$

G beneath a dense cover of grass does not correlate well with air temperature. Hourly G can be approximated during daylight periods as:

$G = 0.1 Rn$

and during nighttime periods as:

$G = 0.5 Rn$

The net radiation (Rn) is the difference between the incoming net shortwave radiation (Rns) and the outgoing net longwave radiation (Rnl):

$Rn = Rns - Rnl = [(1-0.23) * SolarRadiation] - [5.669e-08 * (T_h + 273.15)^{4} * LwradCorrection]$

$LwradCorrection =   (0.34 - 0.14 * VP ^{0.5}) * (1.35 * ATMTransmitivity  - 0.35)$

As saturation vapour pressure is related to air temperature, it can be calculated from the air temperature. The relationship is expressed by:

$e_{0}(T_h) = 0.6108 * e^{17.27 * T_h / (T_h + 237.3)}$

The slope of the relationship between saturation vapour pressure and temperature:

$\Delta = 4098.* \frac {e_{0}(T_h)} {(T_h + 237.3)^{2}}$

The actual vapour pressure can also be calculated from the relative humidity.

$e_a = RH_h * e_{0}(T_h)$

$RH_h$ = Relative Humidity 

$$
-
$$



The Penman Montheith Potential Evapotranspiration computation will be active if in basin file

```
EVAPOTRANSPIRATION : 1
```

and the property evapotranspiration is not readed from file.

#### If the user is running with vegetation

than Crop Evapotranspiration is obtained from Potential Evapotranspiration using crop coefficient from Module Vegetation (dependent on crop).

$CropEvapoTransp = PotentialEvapoTransp * CropCoefficient$

Also if the user is running with vegetation a differentiation in Crop Evapotranspiration between Potential Transpiration and Potential Evaporation may be done using LAI:

```
EVAPOTRANSPIRATION_METHOD : 2
```

$PotentialTrans = CropEvapoTransp * (1 - e^{-0.463 * LAI})$

$PotentialEvap = CropEvapoTransp - PotentialTrans$

if the user is running with vegetation and does not want the Crop Evapotranspiration to be separated in Potential Transpiration and Potential Evaporation then

```
EVAPOTRANSPIRATION_METHOD : 1
```

there will not be soil evaporation (in Module PorousMedia) and Crop Evapotranspiration is all transpiration (in Module Vegetation).

#### If the user is not running with vegetation

than Potential Evapotranspiration is all evaporation in the soil surface and will be done in Module PorousMedia.

#### Fluxes to other modules

If they are used, Potential transpiration computed in Module Basin is given to Module Vegetation to compute the effective transpiration. Potential Evaporation computed in Module Basin is given to Module Porous Media to compute the effective evaporation.

## Other Features

## Outputs

### .srb TimeSeries

water\_column\_

$$
m
$$

 - runoff water column height in top of soil (height from soil surface)

water\_level\_

$$
m
$$

 - runoff water altitude (in reference to mean sea level)

Infil.\_Rate\_

$$
mm/hour
$$

 - infiltration rate in the selected cell

Precipitation\_Rate\_

$$
mm/hour
$$

 - precipitation rate in the selected cell

Throughfall\_Rate\_

$$
mm/hour
$$

 - precipitation rate that arrives to soil surface (after leaf interception) in the selected cell

Canopy\_Capacity\_

$$
m
$$

 - leaf capacity to accumulate water (in height) in the selected cell

Canopy\_Storage\_

$$
m
$$

 - leaf accumulated water (in height) in the selected cell

Canopy\_Drainage\_

$$
m
$$

 - leaf dip water (in height) in the selected cell

Reference\_Evapotranspiration\_

$$
mm/h
$$

 - reference evapotranspiration (computed by Penman Motheith or given by the user)

Potential\_Crop\_EVTP\_

$$
mm/h
$$

 - reference evapotranspiration \* Kc

Potential\_Evaporation\_

$$
mm/h
$$

 - potential crop EVTP separated in evaporation (potential as may not occur depending on soil water availability)

Potential\_Transpiration\_

$$
mm/h
$$

 - potential crop EVTP separated in transpiration (potential as may not occur depending on soil water availability)

EvapoTranspiration\_Rate\_

$$
mm/hour
$$

 - actual evapotranspiration in the selected cell (the effevtive water removed)

Actual\_Evaporation\_

$$
mm/h
$$

 - effective evaporation removed (from potential evaporation) - in surface soil

Actual\_Transpiration\_

$$
mm/h
$$

 - effective transpiration removed (from potential transpiration) - along vegetation roots

[Check Mohid Land Heights and Levels to understand the difference between level and height results](http://screencast.com/t/I8We0B4j2uJ)

### WaterBalance TimeSeries

Basin as the module that controls all others is also the responsible for computing mass balances (for water and properties)

There is a file with properties in m3 (Basin Water Balance) and for simplifying computation of acumulated mass balances there is other in m3/s (Basin Water Balance Fluxes).

```
Basin Water Balance Fluxes
Evap_From_Vegetation_m3/s               - evaporation from leafs
Evap_From_Ground_m3/s                   - evaporation from runoff water column
Evap_From_Soil_m3/s                     - evaporation from soil
Rain_Above_Leafs_m3/s                   - total rain
Rain_Uncovered_m3/s                     - rain in uncovered areas
Rain_Covered_m3/s                       - rain in leafs
Leaf_Drainage_m3/s                      - water that drains from leafs to runoff water column
Rain_Below_Leafs_m3/s                   - amount of rain uncovered + leaf drainage
Out_Volume_Channel_m3/s                 - volume that exits watershed in river – outlet and evaporation
Out_Volume_Overland_m3/s                - volume that exits watershed in runoff (if boundary condition on water can exit to outside, else no exit)
Infiltration_m3/s                       - water that infiltrates in soil
OL_Volume_To_Channels_m3/s              - flow from runoff to river (positive from runoff to river, negative from river to runoff)
GW_Volume_To_Channels_m3/s              - flow from porous media to river (positive from porousmedia to river, negative from river to porousmedia)
```

The accumulation in the watershed would be a difference between the inputs to the watershed and outputs from the watershed Accumulation = Sum ((In – Out)\*LocalTimeStep) In -> Rain (Rain\_Above\_Leafs) + Discharges IN (there is no global output of the discharges, if you have input discharges you have to compute it with the values inserted) Out-> Water loss from soil or surface - Evapotranspiration (Evap\_From\_Vegetation + Evap\_From\_Ground + Evap\_From\_Soil) plus water lost from river (Out\_Volume\_Channel) pus water lost from runoff (Out\_Volume\_Overland)

Note: Negative discharges in river or discharges in runoff are not yet accounted in the mass balance (for mass error results). All that discharges have to be accounted if doing the accumulation computation yourself.

## References

## Data File

### Keywords

```
ATMOSPHERE                    : 0/1              [1]          !Use Module Atmosphere
POROUS_MEDIA                  : 0/1              [1]          !Use Module Porous Media
POROUS_MEDIA_PROPERTIES       : 0/1              [1]          !Use Module Porous Media Properties
RUN_OFF                       : 0/1              [1]          !Use Module RunOff
RUN_OFF_PROPERTIES            : 0/1              [1]          !Use Module RunOff Properties
DRAINAGE_NET                  : 0/1              [1]          !Use Module Drainage Netork
DT_DURING_RAIN                : sec.             [60.]        !initial dt that is tried when rains in the time step
EVAPOTRANSPIRATION_METHOD     : integer          [1]          !1- evapotranspiration - everything in crops; 2- separate
                                                              !between evaporation (surface soil) and transpiration (in crops)
EVAP_FROM_CANOPY              : 0/1              [1]          !Evaporate from plant leafs
EVAP_FROM_WATER_COLUMN        : 0/1              [0]          !Evaporate from water column.
EVAP_METHOD                   : integer          [2]          !Evaporation from canopy and water column method
                                                              !1 - latent heat method, 2 - uses potential evaporation
                                                              !3 - no evaporation from canopy and water column
OUTPUT_TIME                   : sec. sec. sec.    -           !Output Time
TIME_SERIE_LOCATION           : char              -           !Path to time serie location file

Computed or user defined potential evapotranspiration
<beginproperty>
NAME                      : reference evapotranspiration
UNITS                     : mm/h
DESCRIPTION               : fao evapotranspiration
DEFAULTVALUE              : 0.0
REMAIN_CONSTANT           : 0
<endproperty>
```

### Sample

```
ATMOSPHERE                : 1
EVAPOTRANSPIRATION        : 1
EVAPOTRANSPIRATION_METHOD : 1
VEGETATION                : 0
POROUS_MEDIA              : 1
POROUS_MEDIA_PROPERTIES   : 1
RUN_OFF                   : 1
RUN_OFF_PROPERTIES        : 1
DRAINAGE_NET              : 0

!OUTPUT_TIME              : 0 7200
TIME_SERIE_LOCATION       : ..\General Data\TimeSeries\TimeSeriesLocation2D_2.dat
VERIFY_MASS               : 1
CONTINUOUS                : 0
DT_DURING_RAIN            : 3600.

EVAP_FROM_CANOPY          : 1
EVAP_FROM_WATER_COLUMN    : 0

<beginproperty>
NAME                      : reference evapotranspiration
UNITS                     : mm/h
DESCRIPTION               : fao evapotranspiration
DEFAULTVALUE              : 0.0
REMAIN_CONSTANT           : 0
<endproperty>
```
