# Mohid Graphical User Interfaces

The graphical user interfaces of the MOHID Water Modeling System are a set of programs to pre and post process the input and output data needed by the different numerical programs. This chapter describes how the user interfaces integrated in the MOHID Water Modeling System work and how they interact with the numerical code.

The core graphical user interfaces are:

- [Mohid GUI](mohid-gui.md) - a graphical user interface to handle file structure organization of all files which are needed to set up a MOHID numerical simulation;
- [Mohid GIS](mohid-gis.md) - a geographical information system which handles spatial and temporal variable data required or produced by MOHID numerical programs;
- [Mohid PostProcessor](mohid-postprocessor.md) - a graphical user interface which displays data stored in HDF files as animation on the screen;
- [Mohid Time Series Editor](mohid-time-series-editor.md) – a graphical user interface which allows the user to visualize in a quick way time series data required or produced by MOHID numerical programs.
- [Mohid Statistics Analyser](mohid-statistics-analyser.md) - an application to perform statistical analysis of [time series](../../io/formats/time-series.md).

Data exchange between MOHID Graphical User Interfaces and MOHID numerical programs is done by files. There are two types of files which are used: (i) ASCII text files which follow a proper formatting and [HDF files](../../io/formats/hdf-file.md). See also [Input Data Formats](../../io/formats/input-data-formats.md).
