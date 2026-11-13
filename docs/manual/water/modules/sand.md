# Module Sand

## Overview

Module Sand was the result of the accumulated experience gain by Hidromod in the framework of many engineering projects [

$$
1
$$

](http://www.hidromod.pt). Many of the projects were done with a old MOHID version programmed in F77. More recently Hidromod develop a new Module using the present [Mohid Programming](../../development/programming/mohid-programming.md) philosophy. The concepts of this module were described an tested in detail by Carmo (2005) (see [

$$
2
$$

](http://www.mohid.com/PublicData/products/Thesis/TFC_miguelcarmo.pdf) in Portuguese). However, if the user wants to know exactly what this module is doing it can take a look at the source code (see - [

$$
3
$$

](http://mohid.codeplex.com/SourceControl/changeset/view/63225#1208126)).

### Algorithm

From a conceptual point of view this module is very simple. The bed load fluxes are computed in a Arakawa-C grid [

$$
4
$$

](http://en.wikipedia.org/wiki/Arakawa_grids#Arakawa_C-grid). In each centre cell (Z point) the transport flow is computes. In order to compute the transport flows the user can choose from a set of formulations presented below. In a second step the Zonal (or X) component flows are interpolated to the U points and the meridional (or Y) components are interpolated to the V points. Finally to estimate the evolution of the sand thickness in the centre cell (Z point) a mass balance is done. The availability of sand in a cell is limited by the depth of the bed rock.

### Transport formulas implemented

In this section the transport formulas implemented in the MOHID system are enumerated. Based in the characteristics of the formulas and Hidromod experience the formulas can be divided by area of application:

**Rivers**

Meyer-Peter, E; Müller, R. (1948)

**Estuaries**

Ackers and White (1973)

**Coastal areas**

Bijker-Battachraya (1968) Van Rijn (1984, 1993) Bailard (1981,1984) Dibajnia (1992)

## Connect the bed load transport

To connect the bed load transport the user need to define the follow keyword in the InterfaceSedimentWater\_x.dat input file

```
SAND_TRANSPORT : 1
```

If the user choose bed load transport formulations dependent of wave parameters needs also to define the follow keywrods:

```
WAVETENSION : 1
```

```
<begin_waverugosity>
DEFAULTVALUE : 0.002
<end_waverugosity>
```

The user needs also to add the follow keywords in the bathymetry file (see - [Bathymetry](../../io/formats/bathymetry.md)). Keyword EVOLUTION used to activate the ability of the [Module GridData](README.md#module-griddata#module-griddata) to write the bathymetry evolution in a [HDF5](../../tools/hdf5.md) file and read also from a [HDF5](../../tools/hdf5.md) (important for the hot start runs). The keyword EVOLUTION\_FILE is used to define the [HDF5](../../tools/hdf5.md) filename. The keyword PROPERTY\_NAME is used to associate a name to the bathymetry array to be stored in the [HDF5](../../tools/hdf5.md) file.

```
EVOLUTION      : 1
EVOLUTION_FILE : MyModelBathymetryEvolution.hdf5
PROPERTY_NAME  : Bathymetry
```

## Main options

In this section the available options in Module Sand are briefly described. The different options provided by the module can be defined through an input data file, similarly to other MOHID modules. The name of the Module Sand for run number x is Sand\_x.dat (see [Mohid GUI](../../tools/gui/mohid-gui.md))

Keyword TRANSPORT\_METHOD is used to choose the transport formula the user wants to test from the list presented above. The options are: no transport, MeyerPeter, Ackers, VanRijn1, VanRijn2, Bailard, Dibajnia, Bijker. Keyword SAND\_DT is the time step used to compute the transport formula and by default it is equal to the global model time step (see - [Choose the model time step](../../reference/README.md#choose-the-model-time-step#choose-the-model-time-step)). The keyword OLD (0 - OFF; 1 - ON) is ON when the user wants to do a hot start (OFF - cold start).

```
TRANSPORT_METHOD : MeyerPeter
SAND_DT          : 60.
OLD              : 0
```

The user can choose if the compute bed loads change or not the bathymetry and consequently the hydrodynamics. If the user wants to take in consideration the effect of bathymetry changes in the other Mohid modules then the keyword BATHYM\_EVOLUTION : 1 must be defined. The keyword BATIM\_DT is used to define the update frequency of the model bathymetry and by default is equal to the global model time step (see - [Choose the model time step](../../reference/README.md#choose-the-model-time-step#choose-the-model-time-step)).

```
BATHYM_EVOLUTION : 1
BATIM_DT         : 60.
```

Regarding keywords used to define boundary conditions, the keyword BOUNDARY is used to define the open boundary condition. The options are 1 (null gradient) and 2 (cyclic boundary - used in academic studies mainly). The keyword DISCHARGES is used to activate a set of sinks or sources of sediments (0 - OFF, 1 - ON). These sinks/sources of sediments are defined in the discharge\_x.dat input file see [Module Discharges](discharges.md). To define a set of sinks/sources of sediments see [How to create discharges in MOHID](../../tools/gui/how-to-create-discharges-in-mohid.md).

```
BOUNDARY : 1
```

```
DISCHARGES : 0
```

The following keywords are used to define the sediments diameter. These are one of the main inputs of transport formulas.

```
<beginD90>
NAME                  : D90
UNITS                 : m
DESCRIPTION           : Diameter below which 90 percent of the particles are finer
FILE_IN_TIME          : ASCII_FILE
FILENAME              : D90Field.dat
<endD90>
```

```
<beginD50>
NAME                  : D50
UNITS                 : m
DESCRIPTION           : Diameter below which 50 percent of the particles are finer
FILE_IN_TIME          : ASCII_FILE
FILENAME              : D50Field.dat
<endD50>
```

```
<beginD35>
NAME                  : D35
UNITS                 : m
DESCRIPTION           : Diameter below which 35 percent of the particles are finer
FILE_IN_TIME          : ASCII_FILE
FILENAME              : D35Field.dat
<endD35>
```

The following keywords are used to define the sediments availability (bed rock concept). The SAND\_MIN keyword is to define the limit (in meters) below which the transport stops to avoid sediment negative thicknesses.

```
<beginrock>
NAME                  : bed rock
UNITS                 : m
DESCRIPTION           : Depth from sediment surface below which there is no sediment to be transported
FILE_IN_TIME          : ASCII_FILE
FILENAME              : BedRock.dat
<endrock>
```

```
SAND_MIN : 0.01
```

The field properties D90, D50, D35 and bed rock are defined using the options available in [Module FillMatrix](../../base/module-fillmatrix.md).

## Advanced options

The following keywords are used to convert sediments from mass per square meter in meters (POROSITY - porosity and DENS\_SAND - sand density).

```
POROSITY : 0.1
DENS_SAND : 2650.
```

The following keywords are used to smooth the solution. The transport formulas are highly non-linear and tend to generate very noisy solutions. Keyword - FILTER\_SCHEME, options: No Filter (filter OFF), Modify Lax (filter ON). This filter distributes 50% of the bathymetry evolution in the adjacent cells. The keyword FILTER\_RADIUS is used to define the radius (in number of cells) that is used to define the adjacent areas under the effect of sediment redistribution.

```
FILTER_SCHEME : No Filter
FILTER_RADIUS : 4
```

Another option is to use a factor to speed up morphodynamic processes. This factor must be used with caution. This factor multiplies by the transport loads

```
TRANSPORT_FACTOR : 1.
```

The maximum limit for the bottom shear stress used in the Meyer Peter formula can be defined with the following keyword:

```
TAU_MAX = 10.
```

To help improve the model performance it the concept of lateral transport was implemented (SMOOTH\_SLOPE : 0 - OFF, 1 - ON). In this case, it is assumed that when the bathymetry gradient perpendicular to the main bed load flow is above a critical slope (CRITICAL\_SLOPE) there is a flow in this direction. This bed load flow is equal to a percentage (FLUX\_SLOPE) of the main flow. This perpendicular bed load flow tend to smooth the slopes. This is important in river flows where the bathymetry gradients perpendicular to the main flow can be very steep.

```
SMOOTH_SLOPE : 0
CRITICAL_SLOPE : 0.1
FLUX_SLOPE : 0.1
```

## Outputs

The user in this module can output: time series of fields ([HDF5](../../tools/hdf5.md) files), time series in a point ([Time Series](../../io/formats/time-series.md) ASCII files) and time series of integral results for an area ([Boxes](../../io/formats/boxes.md) ASCII files).

### Boxes

The user can compute the sediment fluxes between boxes (or areas).

```
BOXFLUXES : 1
```

If the user wants to compute fluxes between areas needs to define a network of boxes following the ASCII MOHID format [Boxes](../../io/formats/boxes.md).

```
BOX_FILENAME : Boxes.dat
```

### Time series

The user can also make outputs of the most significant properties in a point related with bed load following the standard methodology used in the other MOHID modules (see [Module TimeSerie](timeserie.md)).

For example, if a user wants to activate the output time series option needs to add to the input file (Sand\_x.dat) the following keyword:

```
 TIME_SERIES : 1
```

### Maps (HDF5 format)

The user can also make outputs of the most significant properties fields related with bed load following the standard methodology used in the other MOHID modules (see [OUTPUT TIME](../../io/keywords/output-time.md)).

For example, if a user wants field results every hour starting from the beginning of the run it needs to define the following (time is seconds).

```
 OUTPUT_TIME : 0. 3600.
```

## References

AHMED, S. M., SATO, S. (2003); A sheetflow transport model for asymmetric oscillatory flows. Part I: Uniform grain size sediments; Coastal Engineering Journal 45, 321-337.

ACKERS, P.; WHITE, W.R. (1973); "Sediment Transport: New Approach and Analysis". Journal of the Hydraulics Division (ASCE) 99 (11): 2041–2060.

AL SALEM, A. (1993); Sediment transport in oscilatory boundary layers under sheet flow conditions; PhD thesis, Delft Hydraulics, The Netherlands.

BAGNOLD, R. (1966); An approach of sediment transport model from general physics; US Geol. Survey Prof. Paper 422-I.

BAILARD, J. A. (1984); A simplified model for longshore sediment transport. Proceedings of the 19th Coastal Engineering Conference, pp. 1454– 1470.

BAILARD, J. A., INMAN, D. L. (1981); An energetics bedload model for plane sloping beach: local transport. Journal of Geophysical Research 86 (C3), 2035– 2043.

BAYRAM, A., LARSON, M., MILLER, H., KRAUS, N. (2001); Cross-shore distribution of longshore sediment transport: comparison between predictive formulas and field measurements; Coastal Engineering Journal 44, 79– 99.

BIJKER, E. (1968); Littoral drift as function of waves and current; 11th Coastal Eng. Conf. Proc. ASCE; London, UK; pp. 415–435.

CAMENEN, B., LARROUDÉ, P. (2003); Comparison of sediment transport formulae for the coastal environment; Coastal Engineering 48, 111– 132.

DIBAJNIA, M., WATANABE, A. (1992); Sheet flow under nonlinear waves and currents. Coastal Engineering Journal, 2015– 2029.

DU BOYS, P. (1879); Le rhône et les rivières à lit affouillable; Ann; Ponts Chaussées 18 (5), 171– 195.

FERNANDES, L. (2001); Transporte de Poluentes em Estuários; Trabalho Final de Curso da Licenciatura em Engenharia do Ambiente; Instituto Superior Técnico, Universidade Técnica de Lisboa.

FRIJLINK, H. (1952); Discussion des formules de débit solide de Kalinske, Einstein et Meyer-Peter and Muller compte tenue des mesures récentes de transport dans les rivières néerlandaises; 2nd Journal Hydraulique; Société Hydraulique de France, pp. 98– 103.

KOMAR, P. D. (1998); Beach processes and sedimentation; 2nd Ed.; Pearson Education, New Jersey.

LEITÃO, P.C. (2002); Integração de Escalas e de Processos na Modelação do Ambiente Marinho; Dissertação para a obtenção do grau de Doutor em Engenharia do Ambiente; Instituto Superior Técnico, Universidade Técnica de Lisboa.

LIU, Z. (2001); Sediment Transport; Instituttet for Vand, Jord og Miljøteknik; Aalborg Universitet.

MEYER-PETER, E; MULLER, R. (1948); Formulas for bed-load transport. Proceedings of the 2nd Meeting of the International Association for Hydraulic Structures Research. pp. 39–64.

SANCHO, F. (2002); Apontamentos da disciplina de Processos Fluviais e Costeiros, Mestrado em Hidráulica, Recursos Hídricos e Ambiente; Faculdade de Ciências e Tecnologia da Universidade de Coimbra.

SILVA\*, A., NEVES\*\*, R., LEITÃO, J.C. (1997); Modelação de Processos de Transporte por Acção Combinada de Ondas e Correntes; \*HIDROMOD - Modelação em Engª, Ldª; \*\*Instituto Superior Técnico; Lisboa.

SMITH, J., SHERLOCK, A., RESIO, D. (2001); STWAVE: Steady-State Spectral Wave Model. User’s Manual for STWAVE, Version 3.0; ERDC/CHL, US Army Corps of Engineers; Washington, DC.

TRANCOSO, A. R. (2002); Modelling Macroalgae in Estuaries; Trabalho Final de Curso da Licenciatura em Engenharia do Ambiente; Instituto Superior Técnico, Universidade Técnica de Lisboa.

VAN RIJN, L.C. (1984); Sediment transport: Part I: Bed load transport; Part II: Suspended load transport; Part III: Bed forms and alluvial roughness. Journal of Hydraulic Division 110 (10), 1431– 1456; 110 (11) 1613– 1641; 110 (12) 1733-1754.

VAN RIJN, L.C. (1993); Principles of sediment transport in rivers, estuaries and coastal seas. Aqua Publication, The Netherlands, Amsterdam.

WANG, P., EBERSOLE B., SMITH E. (2002); Longshore Sand Transport – Initial Results from Large-Scale Sediment Transport Facility; ERDC/CHL, US Army Corps of Engineers, Washington, DC.

WINTER, C. (2004); Perfomance of sediment transport models in tidal environments, Workshop HWK, Delmenhorst.
