# Tide Preview

**Tide Preview** allows to generate a tide timeserie from the harmonic consituents of a given tidal gauge. The harmonic components are given by the Fes95 or the [Fes2004](fes2004.md) solution, and tidal gauges can be extracted for points near the boundary.

Configuration of the **TidePrevInput.dat** file:

```
  START                       : YYYY MM DD HH MM SS -         !Start time to compute water level
  END                         : YYYY MM DD HH MM SS -         !End time to compute water level
  DT                          : real              -           !Time step to compute water level
  EXPORT_TO_XYZ               : 0/1               0           !Create a XYZ file with gauge locations
  XYZ_FILE                    : char              -           !Name of XYZ file to be created

<begintideprev>
  IN_TIDES                    : char              -           !Path to gauge file
  OUT_FILE                    : char              -           !Path to output water level time serie file
<endtideprev>
```

Sample IN\_TIDES file:

```
<begingauge>
NAME        : nth test
LONGITUDE   : -16.0000 33.0000 15.0000
LATITUDE    : 15.0000 44.0000 35.0000
METRIC_X    : -16.5543
METRIC_Y    : 15.7432
REF_LEVEL   :         2.08000
TIME_REF    :        0.000000
M2          :         0.432953       -88.0402
S2          :         0.166823       -51.6513
K1          :        0.0566607       -4.53947
K2          :        0.0460439       -54.6424
N2          :        0.0820980       -108.615
2N2         :        0.0104858       -129.940
O1          :        0.0396148       -104.233
Q1          :        0.0186481       -4.40105
P1          :        0.0121636       -161.633
M4          :         0.000000       0.000000
Mf          :        0.0140184       -1.40625
Mm          :       0.00688141       -4.94457
Mtm         :       0.00290994        2.06745
MSqm        :      0.000401051        4.38925
<endgauge>
```
