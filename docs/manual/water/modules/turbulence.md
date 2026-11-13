# Module Turbulence

## Overview

## Concepts

## Horizontal turbulence

### Constant

### Smagorinsky

### Proportional to depth and to the square of velocity

## Vertical turbulence

### Constant

### Nihoul

### Leendertse and Liu

### Backhaus and Hainbucher

### Pacanowski and Philander

### GOTM

Burchard et al, 1999 – ([http://www.gotm.net](http://www.gotm.net)), a turbulence models library coupled with MOHID, including a k-ε model and Mellor-Yamada second order turbulent closure model (Mellor and Yamada, 1982).

## Other features

## Outputs

### Time series

### Box integration

### Maps (HDF5 format)

To write 3D results use the keyword [OUTPUT\_TIME](../../io/keywords/output-time.md) and define keyword:

```
OUTPUT_HDF           : 1
```

### Statistics

Assign the following keywords:

```
STATISTICS               : 1
STATISTICS_FILE          : .. location of the statistics configuration file

STATISTICS_2D            : 1
STATISTICS_FILE_2D       : .. location of the 2D statistics configuration file
```

Configuration of the statistics file is described in [Module Statistics](statistics.md)

## References
