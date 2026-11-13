# Module HDF5

## Overview

Module HDF5 is an interface between [Mohid Water Modelling System](../overview/README.md#mohid-water-modelling-system#mohid-water-modelling-system) and [HDF5](../tools/hdf5.md) library that reads and writes matricial data into HDF5 format.

The advantages of HDF5 are that the files can store much more data, the file access is much faster and the files can be read even if the run does not reach the end.

## Programming

**Constructing phase**

In the constructing phase the following routines must be called: - GetHDF5FileAccess - ConstructHDF5

In the /Grid group the following data should be stored: - WaterPoints3D - [Bathymetry](../io/formats/bathymetry.md) - ConnectionX - ConnectionY

**Working Cycle**

In the working cycle the following routines should be called: - HDF5SetSize - HDF5WriteData (for Time, ConnectionZ (SZZ), OpenPoints3D and for each property) - HDF5FlushMemory

**Destructor phase**

In the destructing phase the files must be closed. To do so one should call: - KillHDF5

## Errors

Frequent errors are found as:

```
HDF5-DIAG: Error detected in HDF5 library version: 1.6.4 thread 0.  Back trace f
ollows.
  #000: D:\muqun\hdf5164-release-fortran\hdf5\src\H5G.c line 297 in H5Gopen(): g
roup not found
    major(10): Symbol table layer
    minor(48): Object not found
```

## See also

[HDF5 file](../io/formats/hdf-file.md)
