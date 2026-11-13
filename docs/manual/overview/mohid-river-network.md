# Mohid River Network

## Overview

Mohid River Network (MRN) is a standalone version of Mohid Land that computes flow and property transport only the [Drainage Network](../water/modules/drainagenetwork.md). MRN capabilities are described in [Trancoso, et al (2009)](http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6V78-4VKMW36-4&_user=2459750&_coverDate=04%2F01%2F2009&_alid=1203175416&_rdoc=1&_fmt=high&_orig=search&_cdi=5836&_docanchor=&view=c&_ct=2&_acct=C000057394&_version=1&_urlVersion=0&_userid=2459750&md5=37512e3e9fe68a06a987760bc7ecbf62)

## Keywords

When running properties with T90 decay, solar radiation is needed. Put it in the Discharges format and add the following keywords to the [Model\_X.dat](../land/how-to-create-a-mohid-river-project-step-by-step.md) input file:

```
SOLARRADIATION_FILE     : ..\..\GeneralData\SolarRadiation.dat
SOLARRADIATION_COLUMN   : 2
```

If runnig properties with SURFACE\_FLUXES, other atmospheric variables should also be given with the following keywords:

```
AIR_TEMPERATURE_FILE     : ..\..\GeneralData\MeteoData.dat
AIR_TEMPERATURE_COLUMN   : 2

CLOUD_COVER_FILE         : ..\..\GeneralData\MeteoData.dat
CLOUD_COVER_COLUMN       : 3

RELATIVE_HUMIDITY_FILE   : ..\..\GeneralData\MeteoData.dat
RELATIVE_HUMIDITY_COLUMN : 4

WIND_SPEED_FILE          : ..\..\GeneralData\MeteoData.dat
WIND_SPEED_COLUMN        : 5
```

## UserManual

[How to create a MOHID River Project Step-by-Step](../land/how-to-create-a-mohid-river-project-step-by-step.md)
