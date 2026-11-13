# Mohid Input/Output

In MOHID, input and output is made via ASCII and [HDF5](../tools/hdf5.md) files. ASCII files are handled by [Module EnterData](../base/module-enterdata.md) and HDF5 files by [Module HDF5](../base/module-hdf5.md). Both the ASCII files and HDF5 files have [specific format conventions for MOHID](../io/formats/input-data-formats.md). These conventions were defined in order to reduce [input data errors](../io/input-data-errors.md) and standardize the output formats so that [Mohid Graphical User Interfaces](../tools/gui/mohid-graphical-user-interfaces.md) can use them straightforwardly.

## Input

To input data into MOHID please see [Module FillMatrix](../base/module-fillmatrix.md) and [Module EnterData](../base/module-enterdata.md).

## Output

### HDF outputs

- [OUTPUT\_TIME](../io/keywords/output-time.md)
- [SURFACE\_OUTPUT\_TIME](../io/keywords/surface-output-time.md)

### Time series outputs

[Module TimeSerie](../water/modules/timeserie.md)
