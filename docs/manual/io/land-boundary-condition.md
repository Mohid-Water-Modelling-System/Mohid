# Land boundary condition;

The land boundary presently in the MOHID system is imposed in the form of water and mass fluxes.

## Hydrodynamic\_x.dat

In the hydrodynamic input data it's necessary to turn on the water discharge option. By default no water discharge is imposed (WATER\_DISCHARGES : 0).

```
WATER_DISCHARGES : 1
```

## WaterProperties\_x.dat

For each property is possible to impose a mass discharge. By default no mass discharge is imposed (DISCHARGES : 0).

```
    <beginproperty>
     DISCHARGES     : 1
    <endproperty>
```

## Discharges\_x.dat

Here's an example on how to define a river discharge. More than one river discharges can be defined. For each discharge corresponds a <begindischarge> / <enddischarge> block.

```
   <begindischarge>
   NAME                    : Quilombo
   DESCRIPTION             : Discharge of Quilombo River (É necessário confirmar a localização)
   I_CELL                  : 142
   J_CELL                  : 92
   K_CELL                  : 1
   ALTERNATIVE_LOCATIONS   : 0
   DATA_BASE_FILE          : ..\..\GeneralData\Discharges\Quilombo.txt
   FLOW_COLUMN             : 2
       <<beginproperty>>
       NAME                    : salinity
       UNITS                   : psu
       DESCRIPTION             : salinity in santos
       DEFAULTVALUE            : 0.5
       <<endproperty>>
       <<beginproperty>>
       NAME                    : cohesive sediment
       UNITS                   : mg/l
       DESCRIPTION             : Coehsive Sediment discharge from Quilombo river
       TIME_SERIE_COLUMN       : 3
       <<endproperty>>
       <<beginproperty>>
       NAME                    : temperature
       UNITS                   : mg/l
       DESCRIPTION             : temperature discharge
       TIME_SERIE_COLUMN       : 4
       <<endproperty>>
   <enddischarge>
```
