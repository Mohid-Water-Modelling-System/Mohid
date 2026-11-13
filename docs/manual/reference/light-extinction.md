# Light extinction

This is a transcript of the [MOHID Light Extinction Technical Manual](http://www.mohid.com/PublicData/Products/Manuals/manual_light_param_mohid_v1.pdf) dated of December, the 6th, 2011.

In this wiki:

- Original Figures are missing
- Some tables must be formatted

## Document outline

The purpose of this document is to describe how solar light is modelled and the light-related processes in the MOHID Water Modelling System, in terms of rationale, options, and parameterizations. Modules, routines, keywords and input files are described and their role explained, so the user can have supplementary insight on the model functioning when modelling scenarios where processes such as surface fluxes, primary production or the effect of light in the decay of microorganisms is used. The processes involving light radiation in the water are divided in three main sections in this document: solar radiation reaching the top of the water mass, light ambient in the water column and light limitation for primary producers. Although linked, these processes are controlled in distinct modules inside MOHID. The document has the following arrangement:

- Section 2 – link between the modules that deal with light radiation in MOHID;
- Section 3 - modelling options for the solar radiation that reaches the top of the water mass
- Section 4 - model options to estimate the extinction coefficient for longwave and shortwave radiation, and model output for light radiation and extinction factor in the water;
- Section 5 - light limitation and control on primary production and algorithms used by the model;
- Section 6 – list of references;
- Section 7 – additional information on the model parameterization.

Users are advised to consult the source code for insight. Additional or detailed information on the KEYWORDS mentioned in this document can be found at www.mohid.com.

## Light ambient in the water

The solar energy available to the oceans is only a part of the total energy irradiated by the sun. From the total incoming solar energy (100%), the 31% is reflected by clouds and atmosphere (Figure 2 1), the 16% is absorbed by atmosphere, water vapour and aerosols, and the 4% is reflected by the surface of the water (albedo). The remaining 49% is available to the oceans (also the 16% of energy stored in the atmosphere might become eventually available to the ocean). The heat of the sun warms the upper ocean waters with different intensity depending on the latitude, the hour of the day, cloud cover and the season. The intensity of the light decreases exponentially with the increasing depth because of absorption and scattering. The sunlight that enters the ocean surface contains different wavelengths from ultraviolet to infrared. Short wavelength penetrate deeper in the water. The rate at which sunlight is attenuated determines the depth which is lighted and heated by the sun.

*Figure 2 1 Simplified schema of heat budget.*

MOHID water modelling system includes the calculation of [processes occurring at the water-air interface](../water/coupling-water-atmosphere-user-manual.md), such as computing wind shear stress, radiation balances, latent and sensible heat fluxes.

### Definitions

Understanding the MOHID water modelling system includes the learning of some specific terms and understanding of information flow between the different units of the modelling system. This section provides some definitions that are useful to understand the light-related processes. Water property: a property of water such as temperature, salinity, biogeochemical quantities (i.e. phytoplankton, zooplankton and nutrients) used as state variables in MOHID; Module: A module is a separately compiled program unit of MOHID that contains definitions and data to share between programs and algorithms to solve specific problems. Several modules handle the light-related processes:

- ModuleWaterProperties: computes water properties processes. In this module the light extinction along the water column is computed;
- ModuleLightExtinction: calculates the light extinction coefficients along the water column;
- ModuleWaterQuality: a zero-dimensional model for primary production, nitrogen and phosphorus cycle;
- ModuleMacroalgae: a zero-dimensional model for macroalgae;
- ModuleLife: a multi-element biogeochemical water-column model with variable stoichiometry;
- ModuleAtmosphere: read from input files and calculate values for atmosphere;
- ModuleInterfaceWaterAir: compute mass and energy fluxes at the water-air interface;
- ModuleFunctions: a set of functions shared among modules to calculate factors, coefficients and other quantities used in MOHID, such as the phytoplankton light limitation factor, among others.

### MOHID information flow

Modelling light-related processes in MOHID implies the export and import of information between modules. In order to respond to the light regime, the water quality module has to receive the instant radiation value for a specific cell. The scheme shown in Figure 2 2 illustrates the relations between the modules in MOHID in terms of flux of information. ModuleWaterProperties gets the radiation values at the top of the water body from ModuleAtmosphere, and then estimates its decay in depth by getting the extinction factor coefficient values from ModuleLightExtinction. If the selected scheme to calculate the extinction coefficient in ModuleLightExtinction needs information such as chlorophyll or cohesive sediment concentrations, this will be provided by ModuleWaterProperties. Ecological modules (ModuleWaterQuality, ModuleLife, ModuleMacroAlgae) are called inside ModuleWaterProperties, and this, in turn, receive the information of available radiation from ModuleWaterProperties. In addition, ModuleWaterQuality calls the function Phytolightextinctionfactor in ModuleFunctions to calculate the light limiting factor. To estimate the downward and upward fluxes of energy at the interface water-air, ModuleWaterQuality calls ModuleInterfaceWaterAir. Feedback mechanisms are involved because radiation levels can, in turn, be affected by the concentration of properties in the water quality modules, such as chlorophyll, DOM, POM or cohesive sediments. In order to address this feedback mechanism, ModuleWaterQuality and ecological models have to export the concentration values to ModuleLightExtinction, which then uses them to calculate the light extinction in water. Then available radiation is sent back to the ecological modules to be used either to estimate the light limitation factor (ModuleWaterQuality) or chlorophyll synthesis (ModuleLife).

*Figure 2 2 Relation between the modules involved in light parameterization in the MOHID. For details see text.*

## Solar radiation at the top of the water body

The parameterization of light in MOHID starts with the amount of radiation that reaches the top layer of the water body, and the associated downward and upward fluxes of energy. There are two modules that control these processes, ModuleAtmosphere and ModuleInterfaceWaterAir, hence there are two input files to setup this processes, Atmosphere.dat and InterfaceWaterAir.dat. The first step to include solar radiation in the simulation consists in declaring its use in Atmosphere.dat, as exemplified in Table 3 1. The effect of cloud cover can also be considered, either as constant or a random value. By default the model will generate its own time series for solar radiation, based on the geographical coordinates of the grid in the header of the bathymetric file (Grid data file). Users must be aware of the time zone of the modelled site and change the GMT reference in Model.dat input file accordingly, if needed (Table 3 2).

*Table 3 1 Example of keyword definition for solar radiation and cloud cover in the Atmosphere.dat input file.*

```
<beginproperty>
NAME                    : solar radiation
UNITS                   : W/m^2
DESCRIPTION             : climatologic solar radiation

FILE_IN_TIME            : TIMESERIE
FILENAME                : ..\..\GeneralData\radiation_timeseries.dat
DATA_COLUMN             : 2

ALBEDO                  : 0.05
FILE_IN_TIME            : NONE
REMAIN_CONSTANT         : 0
<endproperty>

<beginproperty>
NAME                    : cloud cover
UNITS                   : fraction
DESCRIPTION             : cloud cover
MAIN_SOURCE             : MODEL
REMAIN_CONSTANT         : 1
<endproperty>

```

*Table 3 2 GMTREFERENCE definition in the Model.dat input file.*

```
START                         : 2009 5 12 0 0 0
END                           : 2009 8 12 0 0 0
DT                            : 120.
VARIABLEDT                    : 0
GMTREFERENCE                  : 0
SPLITTING                     : Double_Splitting
GMTREFERENCE                  : -3
```

Alternatively, users can define solar radiation with an input time series (as in Table 3 1). Temporal resolution for this forcing is defined inside the time series file. It is also need to setup the InterfaceWaterAir.dat input file properly for the model to deal with the radiation fluxes at the water-air interface. An example of the input data is presented in Table 3 3.

*Table 3 3 Example of property definition for solar radiation in the InterfaceWaterAir.dat input file.*

```
<beginproperty>
NAME                    : surface radiation
UNITS                   : W/m^2
DESCRIPTION             : climatologic solar radiation
INITIALIZATION_METHOD   : CONSTANT
FILE_IN_TIME            : NONE
ALBEDO                  : 0.05
<endproperty>

<beginproperty>
NAME                    : latent heat
UNITS                   : W/m^2
DESCRIPTION             : Calculated latent heat
FILE_IN_TIME            : NONE
<endproperty>

<beginproperty>
NAME                    : sensible heat
UNITS                   : W/m^2
DESCRIPTION             : Calculated sensible heat
FILE_IN_TIME            : NONE
<endproperty>

<beginproperty>
NAME                    : net long wave radiation
UNITS                   : W/m^2
DESCRIPTION             : Calculated net long wave radiation
FILE_IN_TIME            : NONE
<endproperty>

<beginproperty>
NAME                    : upward long wave radiation
UNITS                   : W/m^2
DESCRIPTION             : Calculated upward long wave radiation
FILE_IN_TIME            : NONE
<endproperty>

<beginproperty>
NAME                    : downward long wave radiation
UNITS                   : W/m^2
DESCRIPTION             : Calculated downward long wave radiation
FILE_IN_TIME            : NONE
<endproperty>

<beginproperty>
NAME                    : non solar flux
UNITS                   : W/m^2
DESCRIPTION             : Calculated infrared radiation
FILE_IN_TIME            : NONE
<endproperty>

```

## Available options to calculate light extinction

Light energy reaching the water surface, I0 (Wm-2), can be estimated by MOHID based on a set of coordinates (Lat, Long) or imposed with a time-series. Light is divided into a shortwave component S and a longwave component L, so that . The short wave component of the light is used by the photosynthetic organisms and represents the 40-60% of the incoming radiation I0. In MOHID, long wave and shortwave radiation are calculated separately and can have different extinction coefficients. The model assumes that the light extinction in depth follows the decay given by Lambert-Beer’s law (Figure 4 1), which states that light intensity at depth z is:

```


$$
{{I}_{z}}={{I}_{0}}\,{{e}^{\left( -k\,z \right)}}
$$


          Equation 1
```

where Io is the incident radiation at the surface, Iz the available radiation at depth z and k the extinction coefficient. This law is applied to both shortwave and longwave radiation, so that Equation 3 can be rewritten as:

```
 

$$
{{S}_{z}}={{S}_{0}}\,{{e}^{\left( -{{k}_{s}}z \right)}}
$$


 for shortwave radiation     Equation 2
```

```
 

$$
{{L}_{z}}={{L}_{0}}\,{{e}^{\left( -{{k}_{l}}z \right)}}
$$


 for longwave radiation      Equation 3
```

where ks and kl are the shortwave radiation extinction coefficient and the long wave radiation extinction coefficient, respectively, S0 the shortwave radiation component and L0 is the long wave radiation component at the sea surface. These, in turn, are expressed as:

```


$$
{{S}_{0}}\,=SW\times {{I}_{0}}
$$


       Equation 4
```

```


$$
{{L}_{0}}\,=LW\times {{I}_{0}}
$$


       Equation 5
```

LW and SW are percentages and their sum must be 1. In MOHID, the percentages LW and SW are defined in the input file WaterProperties\_x.dat, by means of the two keywords LW\_PERCENTAGE and SW\_PERCENTAGE (Table 4 1). Because of their distinct characteristics, longwave and the shortwave components of the radiation can have different extinction coefficients. Both are calculated by ModuleLightExtinction. In order to activate the exchanges of heat at the interface water-air, it is necessary to activate the keyword SURFACE\_FLUXES in the file WaterProperties\_x.dat for the property Temperature, as shown in Table 4 2.

*Figure 4 1 Light attenuation in the water column estimated according to Lambert-Beer’s law for optimal conditions (clear water). Water absorption coefficient = 0.04; Io = 300 W m-2.*

*Table 4 1 Example of keyword definition for LW and SW (in blue) in the header of WaterProperties.dat input file.*

```
REFERENCE_DENSITY       : 1027
OUTPUT_TIME             : 0.   7200.
DT_OUTPUT_TIME          : 3600

LW_PERCENTAGE           : .40
SW_PERCENTAGE           : .60
LW_EXTINCTION_COEF      : .67
SW_EXTINCTION_COEF      : 0.04
SW_EXTINCTION_TYPE      : 6
SW_KW                   : 0.04
```

*Table 4 2 Example of keyword definition (in blue) to use in WaterProperties.dat input file.*

```
<beginproperty>
NAME                     : temperature
UNITS                    : ºC
DESCRIPTION              : Temperature of the water
IS_COEF                  : 1
PARTICULATE              : 0
INITIALIZATION_METHOD    : CONSTANT
DEFAULTVALUE             : 20.
DEFAULTBOUNDARY          : 20.
BOUNDARY_INITIALIZATION  : EXTERIOR
BOUNDARY_CONDITION       : 4
OLD                      : 0
ADVECTION_DIFFUSION      : 0
ADVECTION_H_IMP_EXP      : 0
DISCHARGES               : 0
SURFACE_FLUXES           : 1
TIME_SERIE               : 1
OUTPUT_HDF               : 0
<endproperty>
```

### Longwave radiation

Longwave radiation fluxes (upward longwave radiation ﬂux emitted by the sea surface to the atmosphere, and downward longwave radiation ﬂux reaching the sea from the atmosphere) play a major role in the thermal structure of the water body. Unlike shortwave radiation, long wave radiation is not used as an energy source for primary producers, and does not have the same penetrative effect on the water. However, the information given by the longwave radiation is used for the calculation of the temperature in the water. As such, there are only a few options available for longwave radiation extinction in MOHID, specified by the value of keyword LW\_EXTINCTION\_TYPE (Table 4 3).

*Table 4 3 Options for the longwave radiation extinction coefficient in MOHID.*

```
Option	LW_EXTINCTION_TYPE	Description
Constant	1	Constant extinction coefficient
Ascii file	5	Value of the extinction coefficient defined in a time series files
```

### Shortwave radiation

Longwave radiation, alongside with evaporation and conduction, are all surface heat transfer effects that heat occur only at the surface of the water. Shortwave radiation, on the other hand, is a penetrative effect that distributes its heat through a significant range of the water column. Depending on the vertical geometry of the modelled domain (2D or 3D setting), the model estimates the extinction of shortwave radiation within each layer. Since the model uses the control volume approach, the available light for each control volume (or cell) is calculated separately. To calculate the specific amount of solar radiation available for a control volume the model estimates the absorption of light in the water column above (layers) and within the control volume itself, as in Figure 4 2, where the curve represents the extinction of light along the water column according to the Lambert-Beer’s law. This is obtained knowing the total absorption coefficient and the height of the water column above. If the control volume is on the surface, then the incident radiation is considered. There are six methods to calculate shortwave extinction coefficient in MOHID (Table 4 4). The method can be selected by specifying the value of the keyword SW\_EXTINCTION\_TYPE in the file WaterProperties\_x.dat (see Table 4 1). The comparison between the analytical solution of each method and the estimations made by the model in setups with different vertical resolution (number of layers) is plotted in Figure 4 3. The results show that the model is able to reproduce radiation decay in depth satisfactorily. However, users must be aware that low vertical resolution will lead to greater deviation from the analytical solution, meaning that in application where light ambient in the top layers is crucial (e.g., focus on ecological processes), the number of layer must be increased if feasible in terms of computational time.

*Table 4 4 Methods to calculate shortwave light extinction coefficient in MOHID.*

```
Option	SW_EXTINCTION_TYPE	Description
Constant	1	Constant extinction coefficient
Parsons Ocean	2

Portela - Tagus Estuary	3

Combined Parsons-Portela	4

Ascii file	5	ks defined in a time series files
Multiparameter	6

Chla = chlorophyll-a concentration (converted from phytoplankton carbon biomass in ModuleWaterQuality); SPM = suspended sediments, assumed to be cohesive sediments (mg/L) in MOHID; DOC = dissolved organic carbon; POC = particulate organic carbon; Kw = water extinction coefficient; kxn = extinction coefficient for property Xn (n stands for the range of properties that can be included);   = chlorophyll content in producers (ModuleLife).
```

*Figure 4 2 Scheme for calculating the available shortwave radiation in each layer. For simplicity it is assumed that the extinction coefficient ks is constant along the water column.*

#### Constant

A constant value for the extinction factor is used during the entire simulation period. In this case, the WaterProperties input file will have the information described in Table 4 5. The keyword SW\_EXTINCTION\_COEF is the shortwave radiation extinction coefficient in the water (m-1). If the value of the coefficient is not defined, MOHID uses a default value equal to 0.05 m-1. The keyword LW\_EXTINCTION\_COEF is the longwave radiation extinction coefficient in the water (m-1). If the value of the coefficient is not defined, MOHID uses a default value equal to 0.33 m-1. The keywords SW\_PERCENTAGE (dimensionless) and LW\_PERCENTAGE (dimensionless) described in Section 4, if not specified, are set by default to 0.6 and 0.4, respectively.

*Table 4 5 Example of keyword definition to use in WaterProperties.dat input file to implement a constant value for the light extinction coefficient in the water column.*

```
LIGHT_EXTINCTION:1

SW_EXTINCTION_TYPE	:1
SW_EXTINCTION_COEF     :0.04
SW_PERCENTAGE          :0.6
LW_EXTINCTION_TYPE     :1
LW_EXTINCTION_COEF     :0.33
LW_PERCENTAGE          :0.4

<beginproperty>
NAME                     : temperature
UNITS                    : ºC
DESCRIPTION              : Water temperature
IS_COEF                  : 1
INITIALIZATION_METHOD    : CONSTANT
DEFAULTVALUE             : 20.
DEFAULTBOUNDARY          : 20.
BOUNDARY_INITIALIZATION  : EXTERIOR
BOUNDARY_CONDITION       : 4
ADVECTION_DIFFUSION      : 0
ADVECTION_H_IMP_EXP      : 0
DISCHARGES               : 0
SURFACE_FLUXES           : 1
TIME_SERIE               : 0
OUTPUT_HDF               : 1
<endproperty>
```

#### Parsons Ocean

The formulation for this method is taken from Parson et al. 

$$
1
$$

. Only absorption by water molecules and chlorophyll pigments is considered. In this case, in the WaterProperties input file it is necessary to define also the water property phytoplankton as shown in Table 4 6. The constant (0.04) is close to the average extinction coefficient for visible light in clean water. To use this equation in other ambient where dissolved and particulate materials are important to the light extinction, which are not directly related to chlorophyll, the default value of this constant can be changed with the SW\_EXTINCTION\_COEF keyword. A typical value of 0.3 m-1 for reservoirs is described by Chapra 

$$
2
$$

.

Since the Parsons formula accepts values of Chl expressed in mgChla/m3, but the model uses the concentration of phytoplankton expressed as mgC/l, a unit conversion from mgC/l to mg Chla/m3 is done in the model by using a default ratio C:Chla equal to 60. This ratio varies between about 10 and 100 for total phytoplankton 

$$
1
$$

, with 10 corresponding to low light intensities and high temperatures, and 100 corresponding to high light intensities and low temperatures. The default value can be changed using the RATIO\_C\_CHLA keyword in the ModuleWaterProperties.

*Table 4 6 Example of keyword definition to use in WaterProperties\_x.dat input file to implement the Parsons Ocean shortwave radiation extinction coefficient in the water column (in this case the longwave radiation extinction coefficient is constant).*

```
LIGHT_EXTINCTION:1

SW_EXTINCTION_TYPE     :2
SW_EXTINCTION_COEF     :0.04
SW_PERCENTAGE          :0.6

LW_EXTINCTION_TYPE     :1
LW_EXTINCTION_COEF     :0.33
LW_PERCENTAGE          :0.4

<beginproperty>
NAME                     : temperature
UNITS                    : ºC
DESCRIPTION              : Temperature in the water
IS_COEF                  : 1
PARTICULATE              : 0
INITIALIZATION_METHOD    : CONSTANT
DEFAULTVALUE             : 20.
DEFAULTBOUNDARY          : 20.
BOUNDARY_INITIALIZATION  : EXTERIOR
BOUNDARY_CONDITION       : 4
OLD                      : 0
ADVECTION_DIFFUSION      : 0
ADVECTION_H_IMP_EXP      : 0
DISCHARGES               : 0
SURFACE_FLUXES           : 1
TIME_SERIE               : 0
OUTPUT_HDF               : 1
<endproperty>

<beginproperty>
NAME                         : phytoplankton
UNITS                        : mgC/l
DESCRIPTION                  : phytoplankton concentration
IS_COEF                      : 0.001
PARTICULATE                  : 0
………
<endproperty>
```

#### Portela - Tagus Estuary

This method uses a site-specific formula derived from a model calibration with field data from the Tagus estuary, Portugal 

$$
3
$$

. User must be aware that this method can be adequate for estuaries with similar characteristics, and may be used critically to other sites. Only cohesive sediments and water absorption of radiation is taken into account. Please notice that cohesive sediments are assumed to be suspended particulate matter (SPM) as appear in the formulation. In this case, in the WaterProperties input file it is necessary to define also the water properties phytoplankton and cohesive sediment, as shown in Table 4 7, otherwise the model will stop and give an error message.

*Table 4 7 Example of keyword definition to use in WaterProperties\_x.dat input file to implement the Portela-Tagus estuary shortwave radiation extinction coefficient in the water column (in this case the longwave radiation extinction coefficient is constant).*

```
LIGHT_EXTINCTION:1

SW_EXTINCTION_TYPE	 :3
SW_EXTINCTION_COEF     :0.04
SW_PERCENTAGE          :0.6

LW_EXTINCTION_TYPE     :1
LW_EXTINCTION_COEF     :0.33
LW_PERCENTAGE          :0.4

<beginproperty>
NAME                     : temperature
UNITS                    : ºC
DESCRIPTION              : Temperature in the water
IS_COEF                  : 1
INITIALIZATION_METHOD    : CONSTANT
DEFAULTVALUE             : 20.
DEFAULTBOUNDARY          : 20.
BOUNDARY_INITIALIZATION  : EXTERIOR
BOUNDARY_CONDITION       : 4
ADVECTION_DIFFUSION      : 0
ADVECTION_H_IMP_EXP      : 0
DISCHARGES               : 0
SURFACE_FLUXES           : 1
TIME_SERIE               : 0
OUTPUT_HDF               : 1
<endproperty>

<beginproperty>
NAME                         : phytoplankton
UNITS                        : mgC/l
DESCRIPTION                  : phytoplankton concentration
IS_COEF                      : 0.001
PARTICULATE                  : 0
………
<endproperty>

<beginproperty>
NAME                         : cohesive sediment
UNITS                        : mg/l
DESCRIPTION                  : No description was given.
IS_COEF                      : .001
PARTICULATE                  : 1
……………
<endproperty>
```

#### Combined Parsons-Portela

This method combines the formulation of the Parsons Ocean and Portela – Tagus estuary methods described above. This method has been developed to combine the effects of cohesive sediments and chlorophyll in the same expression. The consideration presented in Portela – Tagus estuary method must be considered also when using this method.

#### Ascii file

This option enables the user to define the extinction coefficient by using a single value (constant for the entire simulation) or a time-series with variable values (temporal resolution also defined by the user). In this case it is necessary to define the path of the input file by using a keyword SW\_LW\_EXTINCTION\_FILE. In the file the time series will be organized in column and the column where the shortwave radiation extinction coefficient is present, should be defined with the keyword SW\_EXTINCTION\_COLUMN. The column where the longwave radiation extinction coefficient is present should be defined with the keyword LW\_EXTINCTION\_COLUMN see Table 4 8

*Table 4 8 Example of keyword definition to use in WaterProperties\_x.dat input file to use known values for the shortwave radiation extinction coefficient in the water column (in the example the longwave radiation extinction coefficient is constant).*

```
LIGHT_EXTINCTION:1

SW_EXTINCTION_TYPE     : 5
SW_LW_EXTINCTION_FILE  : C:\ TestSimulation\GeneralData\RadiationTimeSerie.dat
SW_EXTINCTION_COLUMN   : 2

LW_EXTINCTION_TYPE     :1
LW_EXTINCTION_COEF     :0.33
LW_PERCENTAGE          :0.4

<beginproperty>
NAME                     : temperature
UNITS                    : ºC
DESCRIPTION              : Temperature in the water
IS_COEF                  : 1
PARTICULATE              : 0
INITIALIZATION_METHOD    : CONSTANT
DEFAULTVALUE             : 20.
DEFAULTBOUNDARY          : 20.
BOUNDARY_INITIALIZATION  : EXTERIOR
BOUNDARY_CONDITION       : 4
ADVECTION_DIFFUSION      : 0
ADVECTION_H_IMP_EXP      : 0
DISCHARGES               : 0
SURFACE_FLUXES           : 1
TIME_SERIE               : 0
OUTPUT_HDF               : 1
<endproperty>
```

#### Multiparameter

The multiparameter approach to estimate the extinction coefficient has been developed specifically for ModuleLife 

$$
4
$$

, but can be used with ModuleWaterProperties and ModuleWaterQuality. The rationale in this method consists in summing up the individual contribution of each component that may affect light (cohesive sediments, DOC, POC, etc.). These components (state-variables) are defined by the user, by means of keywords for activating the property (LIGHT\_EXTINCTION) and for defining its specific extinction coefficient (EXTINCTION\_PARAMETER) (see example in Table 4 9). There must be a matching in units of the state-variables and respective extinction coefficients must be the same, for example, mg/L and m-1 / mg/L, respectively.

*Table 4 9 Example of keyword definition (in blue) to use multiparameter option in WaterProperties.dat input file (to use with ModuleLife).*

```
<beginproperty>
NAME                     : diatoms chlorophyll
UNITS                    : mg Chl / m3
DESCRIPTION              : Chla content in diatoms
ADVECTION_DIFFUSION      : 1
DISCHARGES               : 0
VERTICAL_MOVEMENT        : 0
PARTITION                : 0
LIFE                     : 1
LIGHT_EXTINCTION         : 1
EXTINCTION_PARAMETER     : 0.05	! m-1/(mgChla/m3)
TIME_SERIE               : 1
OUTPUT_HDF               : 1
<endproperty>

```

### Model output

If the model is setup to use radiation fluxes or to calculate shortwave solar radiation extinction in the water column for ecological or faecal decay model, then the output for these parameters is made by default. The model results for shortwave solar radiation and shortwave solar radiation extinction coefficient for each cell can be seen in the HDF5 result file (Figure 4 4). Also, in the \*.srw time series file for water properties results it is possible to plot shortwave radiation for each cell (average and at the top of the cell), and shortwave radiation extinction coefficient.

Simulation conditions S0 = 300 Wm-2; SW percentage = 60; Albedo = 0; Cloud cover = 0; Max. depth = 25 m; Phytoplankton = 0.03 mgC/L; ; C:Chla = 60; Cohesive sediments = 10 mg/L

*Figure 4 3 Shortwave radiation attenuation in depth calculated with the option currently available in the MOHID model – comparison between the analytical solution and estimations made with different number of layers.*

*Figure 4 4 Snapshot of the HDF5 results window. The options for plotting shortwave solar radiation and shortwave solar radiation extinction can be seen in the menu in the left panel.*

*Figure 4 5 Snapshot of the time-series (\*.srw) results window, with the selected options for plotting average shortwave radiation (per layer), shortwave extinction coefficient and shortwave radiation at the top of the cell.*

## The role of light in ecological models

The ecological models that are built-in in MOHID have two rather different approaches to address the role of light in primary production. The differences arise from the fact that the models (WaterQuality and Life) have been developed around different modelling philosophies. Model WaterQuality is based on the NPZD (Nutrient-Phytoplankton-Zooplankton-Detritus) concept where N is used as the currency between state-variables 

$$
6
$$

. These early ecological model usually don’t consider Chla as a state-variable, and phytoplankton is expressed in C-biomass. Light limitation is usually a function of light availability and an optimal value. Model Life, on the other hand, is based on a more complex philosophy, where organism biomass is defined in terms of elemental composition (C, N and P). For primary producers there are some additional components such as Si and Chla. These models have their origin in the ERSEM model 

$$
7-10
$$

, and among the developments they have undergone is the introduction of algorithm for Chla synthesis. The result is that light control the growth of producers via the production of Chla and the instant ratio of C:Chla of each producer group.

### ModuleWaterQuality

Phytoplankton growth (μ) in this model is a function of a maximum growth rate (μmax) and the limitation imposed by light, nutrient and temperature, according to:

```


$$
\mu ={{\mu }_{\max }}\,{{\Omega }_{T}}\,{{\Omega }_{N}}\,{{\Omega }_{L}}\,
$$


       Equation 6
```

The dimensionless limitation factors vary between 0 (total limitation) and 1 (no limitation). The light limitation factor (ΩL) defines the relationship between the ambient light levels and primary producers’ photosynthetic rate. The photosynthetic rate increases with the light intensity until a maximum photosynthetic rate Pmax is reached (at optimal shortwave radiation Sopt). For values higher than Sopt, the photosynthetic rate decreases. The photosynthetic response to the light is given by the Steele’s formula:

```


$$
\frac{P}{{{P}_{\max }}}=\frac{{{S}_{z}}}{{{S}_{opt}}}{{e}^{\left[ 1-\frac{{{S}_{z}}}{{{S}_{opt}}} \right]}}
$$


       Equation 7
```

where P is the photosynthetic rate (d-1), Pmax the maximum photosynthetic rate (d-1), Sz the shortwave radiation (Wm-2) at depth z (m). From this, the light limiting factor (ΩL) is calculated according to Equation 8 in function PhytoLightLimitationFactor inside ModuleFunction, being called by ModuleWater Quality.

```
LightLim =\frac{1}{z}\int\limits_{0}^{z}{\frac{P}{{{P}_{\max }}}dz
 =\frac{{{e}^{1}}}{z{{k}_{s}}}\left[ {{e}^{\left( -\frac{{{S}_{0}}}{{{S}_{opt}}}{{e}^{-{{k}_{s}}z}} \right)}}-{{e}^{\left( -\frac{{{S}_{0}}}{{{S}_{opt}}} \right)}} \right]}
       Equation 8
```

This formula to calculate the limitation of the light on the growth of phytoplankton is used by ModuleWaterQuality and ModuleMacroalgae, and its functional response is exemplified in Figure 5 1. For a detailed deduction of Equation 8 please see Section 7.

*Figure 5 1 Light limitation calculated with Equation 8 for different values of optimal light intensity (Sopt). Values used: ks = 0.04; z = 2.*

### ModuleMacroAlgae

Light limitation to algal growth inside ModuleMacroAlgae is calculate with the same function (PhytoLightLimitationFactor inside ModuleFunction) used by ModuleWaterQuality for phytoplankton (Equation 8). The differences between macroalgae and phytoplankton can be addressed using different parameter values for Sopt. More information on this module can be found in 

$$
5
$$

.

### ModuleLife

In this module the phytoplankton growth rates (C-fixation) are determined by available light and nutrients using a modified form of a growth model taken from the literature 

$$
11
$$

. The influence of light depends not just on available radiation but also on the instant C:Chla ratio of each producer group. For a comprehensive description of this model and its parameterization please consult 

$$
4
$$

.

## References

1. Parsons, T., M. Takahashi, and G. Hargrave (1984). Biological Oceanographic Processes. Pergamon Press, New York, 330p.
2. Chapra, S. C. (1997). Surface Water-Quality Modeling. USA: McGraw-Hill Companies.
3. Portela, L.I. (1996). Mathematical modelling of hydrodynamic processes and water quality in Tagus estuary. PhD Thesis, Universidade Técnica de Lisboa, Instituto Superior Técnico, Lisboa, Portugal.
4. Mateus, M. (2006). A process-oriented biogeochemical model for marine ecosystems: development, numerical study, and application. PhD Thesis, Universidade Técnica de Lisboa, Instituto Superior Técnico, Lisboa, Portugal. 252p.
5. Trancoso, A. R., S. Saraiva, L. Fernandes, P. Pina, P. Leitao and R. Neves (2005). Modelling macroalgae using a 3D hydrodynamic-ecological model in a shallow, temperate estuary. Ecological Modelling 187(2-3): 232-246.
6. EPA. Rates, constants, and kinetics formulations in surface water-quality modeling. Technical Report EPA/600/3-85/040, US Environmental Protection Agency, 1985.
7. Baretta, J. W., W. Ebenhoh and P. Ruardij (1995). The European-Regional-Seas-Ecosystem-Model, a Complex Marine Ecosystem Model. Netherlands Journal of Sea Research 33(3-4): 233-246.
8. Baretta-Bekker, J. G., J. W. Baretta and E. K. Rasmussen (1995). The Microbial Food-Web in the European-Regional-Seas-Ecosystem-Model. Netherlands Journal of Sea Research 33(3-4): 363-379.
9. Baretta-Bekker, J. G., J. W. Baretta and W. Ebenhoh (1997). Microbial dynamics in the marine ecosystem model ERSEM II with decoupled carbon assimilation and nutrient uptake. Journal of Sea Research 38(3-4): 195-211.
10. Baretta-Bekker, J. G., J. W. Baretta, A. S. Hansen and B. Riemann (1998). An improved model of carbon and nutrient dynamics in the microbial food web in marine enclosures. Aquatic Microbial Ecology 14(1): 91-108.
11. Geider, R. J., H. L. MacIntyre and T. M. Kana (1998). A dynamic regulatory model of phytoplanktonic acclimation to light, nutrients, and temperature. Limnology and Oceanography 43(4): 679-694.
12. Mateus, M, I. A. Kenov, MOHID Light Extinction Manual (2011)

## Annex I

As mentioned before, the shortwave radiation Sz at the depth z is calculated with the Lambert-Beer law (Equation 2). By replacing in Equation 7 we get:

```


$$
\frac{P}{{{P}_{\max }}}=\frac{{{S}_{0}}\,{{e}^{-{{k}_{s}}z}}}{{{S}_{opt}}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}{{e}^{-{{k}_{s}}z}} \right)}}
$$


       Equation 9
```

The average value of the limiting factor at a depth z is given by the ratio between the integral of P/Pmax and z:

```
LightLim=\frac{1}{z}\int\limits_{0}^{z}{\frac{P}{{{P}_{\max }}}dz
=\frac{1}{z}\int\limits_{0}^{z}{\frac{{{S}_{0}}}{{{S}_{opt}}}{{e}^{-{{k}_{s}}z}}{{e}^{\left( 1-\frac{{{S}_{0}}}{{{S}_{opt}}}{{e}^{-{{k}_{s}}z}} \right)}}dz}}
       Equation 10
```

The solution of the integral leads to the expression in Equation 8.

### Integration

To solve the integral in Equation 10 we start by calculating the solution of the indefinite integral,

```


$$
\int{\frac{P}{{{P}_{\max }}}dz}=\int\limits_{{}}^{{}}{\frac{{{S}_{0}}}{{{S}_{opt}}}{{e}^{-{{k}_{s}}z}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}{{e}^{-{{k}_{s}}z}} \right)}}dz}
$$


       Equation 11
```

and then multiply and divide the expression inside the integral by –ks,

```


$$
\int{\frac{P}{{{P}_{\max }}}dz}
$$




$$
=\int\limits_{{}}^{{}}{\frac{-{{k}_{s}}}{-{{k}_{s}}}\frac{{{S}_{0}}}{{{S}_{opt}}}{{e}^{-{{k}_{s}}z}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}{{e}^{-{{k}_{s}}z}} \right)}}dz}
$$


       Equation 12
```

By applying the rules for integration by parts, we can replace:

```


$$
u={{e}^{-{{k}_{s}}z}}
$$


; 

$$
du=-{{k}_{s}}{{e}^{-{{k}_{s}}z}}dz
$$


; 

$$
\int{\frac{P}{{{P}_{\max }}}dz}=\int\limits_{{}}^{{}}{-\frac{{{S}_{0}}}{{{k}_{s}}{{S}_{opt}}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}u \right)}}du}
$$


        Equation 13
```

Being a constant ks can be brought outside the integral,

```


$$
\int{\frac{P}{{{P}_{\max }}}dz}=\frac{1}{{{k}_{s}}}\int\limits_{{}}^{{}}{-\frac{{{S}_{0}}}{{{S}_{opt}}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}u \right)}}du}
$$


       Equation 14
```

and by applying the rules for integration by parts, it is possible to replace:

```


$$
s=1-\frac{{{S}_{0}}}{Sopt}u
$$


 and 

$$
ds=-\frac{{{S}_{0}}}{Sopt}du
$$


; 

$$
\int{\frac{P}{{{P}_{\max }}}dz}=\frac{1}{{{k}_{s}}}\int\limits_{{}}^{{}}{{{e}^{s}}ds}
$$


      Equation 15
```

Considering that the integral of es is es :

```


$$
\int{\frac{P}{{{P}_{\max }}}dz}=\frac{1}{{{k}_{s}}}{{e}^{s}}
$$


       Equation 16
```

By substituting and we obtain:

```


$$
\int{\frac{P}{{{P}_{\max }}}dz}=\frac{1}{{{k}_{s}}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}{{e}^{-{{k}_{s}}z}} \right)}}
$$


       Equation 17
```

Finally, the definite integral between the surface (z=0) and the depth z can be calculated as:

```


$$
\int\limits_{z}^{0}{\frac{P}{{{P}_{\max }}}dz}
$$




$$
={{\left[ \frac{1}{{{k}_{s}}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}{{e}^{-{{k}_{s}}z}} \right)}} \right]}_{z=z}}
$$




$$
-{{\left[ \frac{1}{{{k}_{s}}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}{{e}^{-{{k}_{s}}z}} \right)}} \right]}_{z=0}}
$$





$$
=\frac{1}{{{k}_{s}}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}{{e}^{-{{k}_{s}}z}} \right)}}-\frac{1}{{{k}_{s}}}{{e}^{\left( 1-\frac{{{S}_{0}}}{Sopt}{{e}^{(0)}} \right)}}
$$




$$
=\frac{{{e}^{1}}}{{{k}_{s}}}\left[ {{e}^{\left( -\frac{{{S}_{0}}}{Sopt}{{e}^{-{{k}_{s}}z}} \right)}}-{{e}^{\left( -\frac{{{S}_{0}}}{Sopt} \right)}} \right]
$$


       Equation 18
```
