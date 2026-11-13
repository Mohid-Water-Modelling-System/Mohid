# Settling velocity

The FreeVerticalMovement is used to add to the vertical flow a vertical movement. Presently the module is prepared to impose a constant settling velocity or a settling velocity function of salinity concentration and fine sediments concentration. In the future other types of vertical movement will be computed like phytoplankton up and down movements to optimize the nutrients uptake and/or to reduce mortality due to grazing.

## WaterProperties\_x.dat

For each particulate property it's possible to define a settling velocity condition. By default no settling velocity is assumed.

```
    <beginproperty>
     VERTICAL_MOVEMENT : 1
    <endproperty>
```

## FreeVerticalMovement\_x.dat

An example presented below defines for the cohesive sediment property a constant settling velocity of 1e-5 m/s or 0.01 mm/s. For salinities below 3 psu the settling velocity is zero. The advection transport associated with settling is computed implicitly.

```
<beginproperty>
NAME                           : cohesive sediment
WS_TYPE                        : 1
SALTINT                        : 1
SALTINTVALUE                   : 3
FREEVERT_IMPEXP_ADV            : 0
WS_VALUE                       : 1e-5
DEPOSITION		       : 1
<endproperty>
```
