# Module Functions

## Overview

Module Functions is a set of scientific mathematical functions or routines that are used by various others classes but that did not fit as specific methods of a class. This way this class can be seen as a run-time mathematical library included in the model.

## List of functions

### MOHID Utilities

- SetMatrixValue
- ReadTimeKeyWords
- ConstructPropertyID
- Check\_Hydrodynamic\_Property
- Check\_Water\_Property
- GetDataOnlineString
- IsOdd
- Chunk\_K
- Chunk\_J
- TimeToString
- ChangeSuffix
- SortNumerically\_3D

### Math & Solvers

- InvSingularDiagMatrix2D
- CholeskyFactorization
- CholLinSystemSolver
- TriangLinSystemSolver
- THOMAS\_2D
- THOMAS\_3D
- THOMASZ
- tridag
- tridag\_cyclic
- CGS2D
- LISOLVE
- SIPSOL
- BICGSTAB2D
- Secant
- Normcrossprod
- SLPMIN
- SLPMIN2

### Advection & Diffusion & Boundary Conditions

- ComputeAdvectionFace
- ComputeAdvection1D\_V2
- ComputeAdvection1D
- ComputeAdvection3D
- ComputeDiffusion3D
- OrlanskiCelerity2D

### Heat & Mass transfer

- SaturatedVaporPressure
- LongWaveDownward
- LongWaveUpward
- LatentHeat
- SensibleHeat
- AerationFlux

### Oxygen Saturation

- OxygenSaturation
- OxygenSaturationHenry
- OxygenSaturationCeQualW2

### Density & Specific heat

- Density
- Sigma
- SigmaWang
- SigmaUNESCO
- SigmaLeendertse
- SigmaUNESCOPressureCorrection
- SigmaJMD95PressureCorrection
- SigmaMel96PressureCorrection
- **ConvertTemperature** - converts in-situ temperature into potential temperature
- SpecificHeatUNESCO
- CalcPotentialEnergy

### Grids & Coordinates

- RODAXY
- FromCartesianToGrid
- FromGridToCartesian
- ComputeGridZone
- GetLambda0
- GetEllipsoid
- LatLonToUTM
- UTMToLatLon
- LatLonToLambertSP2

### Interpolation

- InterpolateValueInTime
- InterpolateMatrix2DInTime
- InterpolateMatrix3DInTime
- LinearInterpolation
- InterpolateLinearyMatrix2D
- InterpolateLinearyMatrix3D
- InterpolateProfile
- InterpolateProfileR8
- QuadraticInterpolation
- QuadraticInterpolProfile
- PolIntProfile
- polint

### Biogeochemistry & Water quality

- PhytoLightLimitationFactor
- ComputeT90\_Canteras
- ComputeT90\_Chapra
- SettlingVelocity
