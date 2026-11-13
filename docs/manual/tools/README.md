# Tools

Pre- and post-processing programs, graphical interfaces and helper utilities.

## Sections

- [Gui](gui/README.md) (20 pages)

## Pages

- [2gif](2gif.md) — Please put the executable directory in the PATH environment variable.
- [AssimilationZones](assimilationzones.md) — The AssimilationZones is an application used for generating assimilation decay coefficients (also called relaxation times) files, usable in MOHID Wate…
- [Box2HDF5](box2hdf5.md) — Box2HDF5 is a little tool that creates a HDF5 file from a Boxes definition file, a geometry file and a bathymetry file. This hdf5 file serve for check…
- [Convert2netcdf](convert2netcdf.md) — The input data file must be placed in the working folder and must be named Convert2netcdf.dat.
- [ConvertHDF5ToGridData](converthdf5togriddata.md) — The application ConvertHDF5ToGridData converts an HDF5 file into one or more grid data files (it extracts or one instant or all the intants for the de…
- [ConvertToHDF5](converttohdf5.md) — The ConvertToHDF5 is an application which allows the making of several operations, called actions, involving HDF5 files: conversion of data in other f…
- [ConvertToHDF5 WOA2005](converttohdf5-woa2005.md) — The World Ocean Atlas (WOA) 2005 climatology provides results for water temperature, salinity and several water quality and biology properties. The da…
- [ConvertToXYZ](converttoxyz.md) — 1. Map a network drive to the data file (ETOPO, GEBCO, NOAA etc...).
- [CopyMohidProject](copymohidproject.md) — Eventually it is useful to copy a full aplication of Mohid but without the output (the hdf5 files).
- [Digital Terrain Creator](digital-terrain-creator.md) — The Digital Terrain Creator is a small executable used to generate bathymetries from baseline data. Its inputs are a grid file(.grd) and an ascii xyz…
- [DrawFluxesInHDF5](drawfluxesinhdf5.md) — The latest source code is in the codeplex inside the "SmallTools\\DrawFluxesInHDF5" folder. An older source code is in the Mohid\V4 repository, inside…
- [Fes2004](fes2004.md) — Fes2004 is the tide level and harmonic predictor for tide of POC-LEGOS. Development credits go to Thierry Letellier and Florent Lyard.
- [FillMatrix](fillmatrix.md) — The FillMatrix application can create a grid data file or HDF5 file from two possible origins:
- [FilterBathymetry](filterbathymetry.md) — The application FilterBathymetry performs filtering of MOHID bathymetries.
- [GenerateGrid](generategrid.md) — GenerateGrid is a MOHID support tool used to create regular grids with a non-constant grid size. It reads an input data file named GridGenerator.dat.
- [Generating curvilinear grids](generating-curvilinear-grids.md) — Add to your MOHID GIS project the polygon that will be the boundary of your meshing domain
- [Geometry diagnostic](geometry-diagnostic.md) — Once the vertical discretization or geometry is imposed and the bathymetry is chosen, the bottom layer can yield stability problems when using shaved…
- [GLUES HDF5 FILES](glues-hdf5-files.md)
- [H5DUMP](h5dump.md) — H5dump is a small tool that comes with the hdf5 libraries (windows or unix). It allows to dump in ASCII the content from hdf5 files.
- [H5utils](h5utils.md) — H5utils are open-source command-line tools to manipulate HDF5 files. They can be downloadaded from h5utils-1.12.1.tar.gz.
- [HDF5](hdf5.md) — HDF5 is a general purpose library and file format for storing scientific data. It's the standard input/output format of spatial and/or temporal data s…
- [HDF5Exporter](hdf5exporter.md) — The application HDF5Exporter produces time series from data contained in HDF5 files in MOHID format .
- [HDF5Extractor](hdf5extractor.md) — The application HDF5Extractor extracts from an HDF5 file data to produce one or more smaller HDF5 files.
- [HDF5Statistics](hdf5statistics.md) — The application HDF5Statistics produces time statistics in HDF5 format from data existent in HDF5 files.
- [Hudson](hudson.md) — Hudson allows automatic daily builds and daily tests with a reporting system on top. Good to use with MOHID.
- [InvertedBarometer](invertedbarometer.md) — This application creates water level time series for tide gauge locations according with the inverted (inverse) barometer effect.
- [JoinTimeSeries](jointimeseries.md) — JoinTimeSeries is a MOHID support tool designed to join several time series files in a single file.
- [Libfproj4](libfproj4.md) — libfproj4 is a Fortran wrapper for PROJ4 library. It was developed as of the GLIMMER project and can be found in libfproj4-1.0.tar.gz. This wrapper is…
- [Mohid Support Tools](mohid-support-tools.md) — MOHID Water Modelling System has a large variety of support tools used to help creating simulations, preparing and processing input and output data fi…
- [Mohid Swan](mohid-swan.md) — In this wiki it would be commented how to implement the Swan model into the Mohid modelling software.
- [Ncdump](ncdump.md) — ASCII dumps the datafile header only:
- [NOAAnetcdfTOhdf5](noaanetcdftohdf5.md) — The NOAAnetcdfTOhdf5 is an application which allows the conversion of NOAA meteorological reanalysis files in NetCDF format to HDF5 format readable by…
- [PROJ4](proj4.md) — This page explains what is PROJ4 library and how to use it from the command line, in Windows and in Linux.
- [RapidSvn](rapidsvn.md) — RapidSVN is a cross-platform GUI front-end for the Subversion revision system written in C++ using the wxWidgets framework. This project also includes…
- [SmoothBatimNesting](smoothbatimnesting.md) — Animation consisting of three frames showing respectively i) the father bathymetry, ii) the son bathym, iii) and the interpolated bathym.
- [Tide Preview](tide-preview.md) — Configuration of the TidePrevInput.dat file:
- [TortoiseSvn](tortoisesvn.md) — Simply install presto the tortoise svn Tortoise svn. You can also see video screencasts instead of the picture-only tutorial below...
- [WOAnetcdfTOhdf5](woanetcdftohdf5.md) — This program converts the NetCDF files of the World Ocean Atlas 2005 (WOA05) into HDF5, it allows to create HDF5 files for Mohid using the Module LIFE…

## Short entries

### Box2GridData

Box2GridData is a [MOHID support tool](mohid-support-tools.md) used to convert a [boxes](../io/formats/boxes.md) file into [Grid Data](../io/formats/grid-data.md) format. This tool enables the preview of the boxes file before running MOHID, in order to check for possible errors when defining the integration boxes.

### H5diff

**h5diff** is a command-line tool that looks for differences in results in two given hdf5 files.

#### Syntax

```
> h5diff file1 file2
```

#### External references

- [HDF5 tools home](http://www.hdfgroup.org/products/hdf5_tools/)
- [HDF5tools tutorial](http://www.speedup.ch/workshops/w37_2008/HDF5-Tutorial-PDF/HDF5-Tools.pdf)

### HDF5Operator

HDF5Operator is a tool able to compute some operations regarding MOHID results files namely:

- compute differences (absolute and relative) between parametric simulations
- sum different data sets within a single results file multiplied by a factor

### HDFView

HDFView is java based (cross-platform) HDF4 and HDF5 viewer and editor. It allows the browsing of hdf data in numerical tables format. It also allows the manual editing of the values in the tables; [sometimes this is a useful feature](mohid-swan.md).

#### External Links

- [Hdfview download page](http://www.hdfgroup.org/hdf-java-html/hdfview/)

### MohidReservoirOptimization

MohidReservoirOptimization is part of Mohid system. It is used to create bathymetries to be used with the CEQUAL-W2 EPA's program. It is writen in FORTRAN95.

- [User manual](../reference/user-manual.md).
- [Technical manual](../reference/README.md#technical-manual#technical-manual).

### Search

<table style="padding:1em"><tbody><tr><td width="15%"></td><td width="17%"><a href="http://wiki.mohid.com/wikisearch.html"><img id="mohid" src="images/Mohid.gif" alt="google"></a></td><td width="33%"><input type="text" name="q" size="40"></td><td width="10%"><input type="submit" name="sa" value="Search"></td><td width="15%"></td></tr></tbody></table>
