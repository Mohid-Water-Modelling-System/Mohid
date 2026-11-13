# Software to run Mohid

In order to run Mohid you need the following software

## System requirements

- Windows XP or Windows Vista or Windows 7 ([not your operating system of choice?](../reference/developers.md))
- 1 GB ram (2 GB recommended)
- Hardware graphical accelerator
- 250 MB of free hard-disk space is the bare minimum. >10 GB, or more, of free disk space is recommended. (Depending on the modeler's goals, results could require 50 MB of free disk space, but results could also require 500 GB of free disk space).
- Internet access (for documentation and support)

**WORD OF CAUTION: For MOHID GIS to work correctly, you need to configure the windows settings so that the decimal separator be a point. For example: 4.1516 is correct but 4,1516 is incorrect. If the separator isn't a point, then importing data files into MOHID GIS may yield errors.**

**WORD OF CAUTION2: Google earth must show coordinates in decimal degrees.**

## Pre-requisites

- DOT Net framework
- Office web components

## Download MOHID package

- Register, login, download and install the latest version of the MOHID package at the [MOHID main webpage](http://www.mohid.com/). It includes:
    -   [MohidWater](../reference/mohidwater.md), [MOHID Land](../land/mohid-land.md), [GIS](../tools/gui/mohid-gis.md), [GUI](../tools/gui/mohid-gui.md), MOHID PostProcessor, MOHID Time Series Editor
- Requirements for installing Mohid Package:
    -   The programs to be installed in your computer were designed to work on Microsoft Windows 2000 and XP;
    -   To install these package be sure that the .NET Platform v1.1 is installed in your computer. If you do not have the .NET Platform v1.1, please go to: <[http://www.microsoft.com/downloads/details.aspx?FamilyID=262d25e3-f589-4842-8157-034d1e7cf3a3&displaylang=en](http://www.microsoft.com/downloads/details.aspx?FamilyID=262d25e3-f589-4842-8157-034d1e7cf3a3&displaylang=en)\>;
    -   You also need Office Web Components to be installed. They can be found here: [http://www.microsoft.com/downloads/details.aspx?FamilyID=982b0359-0a86-4fb2-a7ee-5f3a499515dd&DisplayLang=en](http://www.microsoft.com/downloads/details.aspx?FamilyID=982b0359-0a86-4fb2-a7ee-5f3a499515dd&DisplayLang=en);
    -   It is strongly advised to install this package in the default Windows Program Files folder;
    -   Colour quality in your screen display settings must be set to 24 bits or more;
    -   OpenGL must be supported by your graphics card;
    -   Be sure to grant Read/Write access to the Temporary folder created in your installation folder (../Mohid/Temp);

## MohidWater primer

- Download and install [Google earth](http://earth.google.com/).
- Download and store the [kml2mohid.exe](http://www.mohid.com/hydrogroup/data_sources/kml2mohid.exe) program.

- Download and uncompress the [ConvertToXyz](http://www.mohid.com/hydrogroup/data_sources/convert.zip) program.
- Download and store the [ETOPO2V2](http://www.ngdc.noaa.gov/mgg/global/etopo2.html) bathymetry data in **raw binary format**, **cell centered**, **little-endian**.

- Download and uncompress the [Mohid-Tide](http://maretec.mohid.com/PublicData/products/Software/Mohid-Tide_v0.4.zip) program.
- Download and store the [Fes2004](http://www.legos.obs-mip.fr/en/share/soa/cgi/getarc/v0.0/index.pl.cgi?contexte=SOA&donnees=maree&produit=modele_fes) tide solution. **Data access**, **fes2004**, **tide**, **tide.fes2004.nce**.

### [Tutorials](tutorials.md)

Once you downloaded and installed everything above, you can read the following [tutorials](tutorials.md).

## Internal links

- [Main Page](../overview/main-page.md)
- [Developers](../reference/developers.md)
