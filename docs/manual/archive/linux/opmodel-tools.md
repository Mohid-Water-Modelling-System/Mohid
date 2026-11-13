# OpModel Tools

## WeeklyMonitor

- Location: **\\\\guillaume**
- Usage: **$ weeklymonitor.bat**
- Description: It returns a list of all the locations, filesizes, filenames and filestatus required to properly execute a run with the operational model for the current RunID.
- Example output:

```
D:\Data\Mercator\perl>weeklymonitor.bat
        machine    filesize                                filename
      guillaume   248608068               MercatorR20070905ANA.hdf5
      guillaume   464610694               MercatorR20070905PRE.hdf5
         kepler   248608068               MercatorR20070905ANA.hdf5
         kepler   464610694               MercatorR20070905PRE.hdf5
      guillaume   679252007                  MercatorR20070905.hdf5
      guillaume   118451163          MercatorPortugalR20070905.hdf5
      guillaume    99323992  MercatorPortugalAnalyzedR20070905.hdf5
      guillaume   108965606       MercatorEstremaduraR20070905.hdf5
      guillaume    93430374 MercatorEstremaduraAnalyzedR20070905.hd
       einstein    56530312        D2_7d_2007083100_2007090700.hdf5
       einstein    56380232        D2_7d_2007090700_2007091400.hdf5
       einstein   107762393        D3_7d_2007083100_2007090700.hdf5
       einstein   107641321        D3_7d_2007090700_2007091400.hdf5
         kepler    56530312        D2_7d_2007083100_2007090700.hdf5
         kepler         N/A        D2_7d_2007090700_2007091400.hdf5
         kepler   107762393        D3_7d_2007083100_2007090700.hdf5
         kepler         N/A        D3_7d_2007090700_2007091400.hdf5
      guillaume         N/A        D2_7d_2007083100_2007091400.hdf5
      guillaume         N/A        D3_7d_2007083100_2007091400.hdf5
      guillaume         N/A              MM5_20070905_Portugal.hdf5
      guillaume         N/A           MM5_20070905_Estremadura.hdf5
 192.168.23.151   118451163          MercatorPortugalR20070905.hdf5
 192.168.23.151    99323992  MercatorPortugalAnalyzedR20070905.hdf5
 192.168.23.151   108965606       MercatorEstremaduraR20070905.hdf5
 192.168.23.151    93430374 MercatorEstremaduraAnalyzedR20070905.hd
 192.168.23.151         N/A           WaterProperties_20070905.hdf5
 192.168.23.151         N/A           WaterProperties_20070905.fin5
 192.168.23.151         N/A              Hydrodynamic_20070905.hdf5
 192.168.23.151         N/A               Hydrodynamic_20070905.fin
 192.168.23.151         N/A           WaterProperties_20070905.hdf5
 192.168.23.151         N/A           WaterProperties_20070905.fin5
 192.168.23.151         N/A              Hydrodynamic_20070905.hdf5
 192.168.23.151         N/A               Hydrodynamic_20070905.fin
      guillaume         N/A           WaterProperties_20070905.hdf5
      guillaume         N/A              Hydrodynamic_20070905.hdf5
      guillaume         N/A    20070905_Portugal_WaterProperties.nc
      guillaume         N/A       20070905_Portugal_Hydrodynamic.nc
      guillaume         N/A           WaterProperties_20070905.hdf5
      guillaume         N/A              Hydrodynamic_20070905.hdf5
      guillaume         N/A 20070905_Estremadura_WaterProperties.nc
      guillaume         N/A    20070905_Estremadura_Hydrodynamic.nc
 data.mohid.com         N/A    20070905_Portugal_WaterProperties.nc
 data.mohid.com         N/A       20070905_Portugal_Hydrodynamic.nc
 data.mohid.com         N/A 20070905_Estremadura_WaterProperties.nc
 data.mohid.com         N/A    20070905_Estremadura_Hydrodynamic.nc

Files completed: 19 / 45
```
