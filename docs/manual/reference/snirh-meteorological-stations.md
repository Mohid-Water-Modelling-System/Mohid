# SNIRH Meteorological Stations

One useful source of meteorological data (to Portugal) to use as input on Mohid is 

$$
SNIRH
$$

 website. There are many stations with historical data available, at least for precipitation and wind speed, but sometimes for solar radiation, humidity and air temperature can also be found.

The most common data will have a daily interval. This means that you will find a single value for the property for each day. This values, in the case of precipitation, is the accumulated rain between the 00:00 hours and 23:59:59 hours of each day (indicated in the "Data" column).

To use this information as input to mohid, as a TIMESERIES, it's necessary to change the start date by 1 day, to reflect the way the MOHID interprets information on timeseries (from the date to backwards)
