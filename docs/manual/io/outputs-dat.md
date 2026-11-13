# Outputs.dat

!IF THIS IS ZERO (0) SWAT RUNS EXACTLY LIKE THE ORIGINAL ONE. IF IT IS ONE (1) THE MODIFED VERSION IS RUNNED

SWAT\_WITH\_CHANGES  : 1

DENITRIFICATION\_THRESHOLD  : 0.99

!OUTPUT\_TIME 

$$
0 86400
$$

  ! if this is on it will create two hdf files. One with LAI, Soil moisture and sediment information per HRU and another with flow per reach. Wkt geometry of hru poligons and wkt with river line, must be provided

!Subbasin level outputs

SUBBASIN\_OUT  : 1  ! connects the output by subbasin to the sub-basins that are in the file under keyword SUBBASIN\_TIME\_SERIE\_LOCATION

REACHES\_OUT  : 1  ! connects the output to the outlet of the sub-basins that are in the file under keyword SUBBASIN\_TIME\_SERIE\_LOCATION

!HRU level outputs METEO\_OUT  : 1 ! connects the HRU output from meteo to the HRUs that are in the file HRU\_TIME\_SERIE\_LOCATION

PLANTS\_OUT  : 1 ! connects the HRU output from plants to the HRUs that are in the file HRU\_TIME\_SERIE\_LOCATION

NITRATE\_OUT  : 1 ! connects the HRU output from nitrate to the HRUs that are in the file HRU\_TIME\_SERIE\_LOCATION

PHOSPHORUS\_OUT  : 1 ! connects the HRU output from P to the HRUs that are in the file HRU\_TIME\_SERIE\_LOCATION

USLE\_OUT  : 1 ! connects the HRU output from usle to the HRUs that are in the file HRU\_TIME\_SERIE\_LOCATION

HRU\_WATER\_BALANCE\_OUT  : 1 ! connects the HRU output from Water Balance to the HRUs that are in the file HRU\_TIME\_SERIE\_LOCATION

NUTRIENT\_PER\_HRU  : 0 ! connects the HRU output from nutrients to the HRUs that are in the file HRU\_TIME\_SERIE\_LOCATION

!DEFINE THE OUTPUT OF THE GLOBAL WATERSHED

GLOBAL\_WATER\_OUT  : 1 ! agrega os resultados do modelo por fluxo de H20 de toda a bacia

GLOBAL\_WATER\_CUMULATIVE  : 1 ! Valor diário se for 0 ou valor acumulado ao longo do periodo de corrida se for 1

GLOBAL\_WATER\_M3  : 0 ! Muda o resultado para m3 em x de mm

GLOBAL\_NUTRIENTS\_OUT  : 0 ! agrega os resultados do modelo por fluxo de nutriente de toda a bacia (diário)

LAND\_USE  : 0 ! agrega os resultados do modelo por uso de solo (diário)

BASIN\_AVERAGE  : 0 ! agrega os resultados de nutrientes por cada fluxo principal de H20

!BETTER NOT CHANGE THIS

HRU\_TIME\_SERIE\_LOCATION  : ..\\res\\TS\_HRU.dat

SUBBASIN\_TIME\_SERIE\_LOCATION : ..\\res\\TS\_SubBasin.dat

GLOBAL\_TIME\_SERIE\_LOCATION  : ..\\res\\TS\_Global.dat

SUMMARY\_TIME\_SERIE\_LOCATION  : ..\\res\\TS\_Summary.dat

LAND\_USE\_TIME\_SERIE\_LOCATION : ..\\res\\TS\_LandUse.dat
