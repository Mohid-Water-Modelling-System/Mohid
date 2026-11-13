# Paraview for MOHID

This page will describe how to open and visualize MOHID HDF5 and NetCDF results using the ParaView software. ParaView is an open-source, multi-platform data analysis and visualization application. ParaView users can quickly build visualizations to analyze their data using qualitative and quantitative techniques.

HDF5 files from MOHID need to be obtained by using executable files compiled since February 2019 that have included the field Corners3D. If earlier, HDF5 files need to be converted to NetCDF file using XXXXXXXXX.

## Required Software

Please install the latest version of Paraview that can be found at [https://www.paraview.org/download/](https://www.paraview.org/download/)

You will also need the Phyton V3 that can be obtained form the Miniconda distribution [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html). It will be needed the h5py package, please open a Anaconda prompt and install it using the following instruction:

```
pip install h5py
```

You will also need to download the package [https://github.com/maretec/MOHID\_python\_tools](https://github.com/maretec/MOHID_python_tools)

## Opening HDF5 files

To be able to open MOHID HDF5 in paraview it needs a descriptor file with the extension .xdmf that is generated using a python tool that is located at [https://github.com/maretec/MOHID\_python\_tools/tree/master/MXDMF](https://github.com/maretec/MOHID_python_tools/tree/master/MXDMF)

It is needed to run the file RunWin.bat ([https://github.com/maretec/MOHID\_python\_tools/blob/master/MXDMF/run/RunWin.bat](https://github.com/maretec/MOHID_python_tools/blob/master/MXDMF/run/RunWin.bat)) or RunLinux.sh ([https://github.com/maretec/MOHID\_python\_tools/blob/master/MXDMF/run/RunLinux.sh](https://github.com/maretec/MOHID_python_tools/blob/master/MXDMF/run/RunLinux.sh)) according to the operational system where the location of the TestFiles need to be edited. The tool also is able to glue files.

Open using XDMF reader

Scale Z should be negative to visualize 3D files due to the bathymetries being defined as positive values.

```
Visualization
Surface - fill surface
Points - show particules
```

To filter land, you need to use the threshold filter, then you select the variable and then the maximum value

## Visualizing variables
