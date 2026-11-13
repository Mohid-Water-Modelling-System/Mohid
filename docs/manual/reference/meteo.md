# Meteo

This chapter explains which meteorological data we have available. Currently there are 2 meteorological models running operationally for Portugal ([MM5](http://www.mmm.ucar.edu/mm5/) and [WRF](http://www.wrf-model.org/)). The higher operational spatial resolution is 9 km but 3km are already being tested. MM5 3 day forecasts for Portugal are published on [http://meteo.ist.utl.pt/](http://meteo.ist.utl.pt/). There is also available forecasts and analysis from the global model [GFS](http://en.wikipedia.org/wiki/Global_Forecast_System) which are used for initial and boundary conditions for MM5, WRF and recently [WWIII](http://www.maretec.mohid.com/ww3/). There is also a database of meteorological observations (surface and upperair).

## MM5

MM5 meteorological surface forecasts are automatically being delivered to MARETEC in MOHID-HDF5 format. This section explains where are the HDF5 files. For other time periods and/or variables, or if you see that there is some file missing or has problems, please email me! (arosa@ist.utl.pt)

### System Operational Configuration

[MM5](http://www.mmm.ucar.edu/mm5/) is running 4 time per day (00,06,12,18 h UTC) producing 3 days forecasts with 3 spatial nested domains: D1 with 81 km resolution, D2 with 27 km and D3 with 9 km ([http://meteo.ist.utl.pt/](http://meteo.ist.utl.pt/)). These simulation were using [GFS](http://en.wikipedia.org/wiki/Global_Forecast_System) grib1 (1º resolution, which is ~80 km at our latitudes). Since August 2007, only 2 domains are simulated (D2 with 27 km and D3 with 9 km) because MM5 uses GFS grib2 output (0.5º resolution, which is ~40 km at our latitudes). Also, since XXX 7 days forecasts are produced four times per day.

### Where are the HDF5 files?

Recent files can be found in [ftp://ftp.mohid.com/](ftp://ftp.mohid.com/) (meteoIST:mm5) where:

1. **[ftp://ftp.mohid.com/mm5\_6h](ftp://ftp.mohid.com/mm5_6h)** stored in **\\\\alexandria\\public\\modelos\\Meteo\_IST\\mm5\_6h**
    -   Each file has the 2nd 6 hours of simulation. When glued, you can have a combination of the most recent available forecast.
    -   simulation 00h : files have outputs from 7:00 to 12:00
    -   simulation 06h : files have outputs from 13:00 to 18:00
    -   simulation 12h : files have outputs from 19:00 to 00:00
    -   simulation 18h : files have outputs from 1:00 to 6:00
2. **[ftp://ftp.mohid.com/mm5\_7d](ftp://ftp.mohid.com/mm5_7d)** stored in **\\\\alexandria\\public\\modelos\\Meteo\_IST\\mm5\_7d**
    -   Simulation 00h with 7 days forecasts. This is only delivered at Friday's to run MOHID Portuguese Coast Operational Model.
3. **[ftp://ftp.mohid.com/SIGEL/mm5\_3d](ftp://ftp.mohid.com/SIGEL/mm5_3d)** (not being stored)
    -   Four daily simulations with 3 days forecasts.

### References to MM5-IST configuration

- Sousa, T., Domingos, J. (2002) Previsão Meteorológica em Portugal Continental utilizando o modelo operacional e de investigação MM5. Dissertação para obtenção do Grau de Mestre em Ecologia, Gestão e Modelação dos Recursos Marinhos, Instituto Superior Técnico, 150pp. [

$$
1
$$

](http://redibericamm5.uib.es/publicacions/any2002/g30_tesi_02.pdf)

- Rodrigo's Msc thesis has a chapter on validation of MM5 forecasts for Guia

- Riflet, G.; Leitão, P.C.; Trancoso, A.R.; Canas,A.; Fernandes, L.; Fernandes, R.; Garcia, A.C.; Mateus, M.; Neves, R.J. (2007). Assessing the quality of a pre-operational model for the portuguese coast. Geophysical Research Abstracts, Vol. 9, 09979, European Geosciences Union, 15 – 20 April, Vienna, Austria

## Maretec archive

### MM5 forecasts

Here's the unexhaustive list of MM5 data outputs in hdf5 at Maretec data servers. There are three type of domains: D1 (~90 km), D2 (~30 km) and D3 (~10 km). However the domains changed during 2006/2007 and so there are the new D1, D2, D3 and the old D1, D2, D3.

- 20060601 - 20061020:
    -   D1, D2, D3 \\\\alexandria\\public\\modelos\\Meteo\_IST\\mm5\_7d\\7d\_old

- 20061020 - 20061124:
    -   D1, D2, D3 \\\\alexandria\\public\\modelos\\Meteo\_IST\\mm5\_3d

- 20061124 - 20061201:
    -   D2 \\\\alexandria\\public\\modelos\\Meteo\_IST\\mm5\_7d

- 20061201 - 20070103:
    -   D2 \\\\alexandria\\public\\modelos\\Meteo\_IST\\mm5\_7d
    -   D3 \\\\alexandria\\public\\modelos\\Meteo\_IST\\Dezembro2006

- 20070103 - 20070706:
    -   D2 \\\\alexandria\\public\\modelos\\Meteo\_IST\\mm5\_7d

- 20070706 - 20080509:
    -   D2, D3 \\\\alexandria\\public\\modelos\\Meteo\_IST\\mm5\_7d

### GFS wind files

1. **\\\\alexandria\\public\\modelos\\Meteo\_IST\\GFS\\**
    -   Daily 1 day wind forecasts(from 0h00 to 21h00) from GFS model covering the North Atlantic
2. **\\\\alexandria\\public\\modelos\\Meteo\_IST\\GFS\\latest\\**
    -   Daily 7 days wind forecasts from GFS model covering the North Atlantic

## Related Pages

- [Module\_InterfaceWaterAir](../water/modules/interfacewaterair.md)
- [ConvertToHDF5](../tools/converttohdf5.md)
