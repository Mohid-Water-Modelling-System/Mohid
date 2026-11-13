# AssimilationPreProcessor

The application **AssimilationPreProcessor** performs the calculation of the covariance structure of data in MOHID format HDF5 files while performing EOF Analysis.

Running options for this application are specified in an [input file](assimilationpreprocessor.md) whose path is indicated in a [nomfich.dat file](../io/nomfich-dat-file.md).

## Introduction

AssimilationPreProcessor calculates from fields in HDF5 files the covariance structure to be used as initial state error covariance in sequential data assimilation operations in MOHID Water.

The fields to be analysed are specified through the definition of the state. The state is a set of variables which characterizes the system to be studied at a particular time instant.

These variables can refer to several hydrodynamic and water properties and represent also several spatial locations (or cells) according to the horizontal and vertical grids defined in the HDF5 files. The spatial locations are specific for each property and therefore can be different according to property.

The time instants contained in the input HDF5 files can be selected using a decimation process. This can be used to reduce computational costs of covariance generation and processing and also if subsequent time instants provide very similar data.

Following the state definition the covariance matrix of the state is calculated. This matrix is then processed to be in a form (the covariance structure) suitable to be used in MOHID Water for the initialization of sequential assimilation processes.

Currently, the covariance matrix is processed only only for SEEK filter sequential data assimilation operations. This calculation involves the performance of Empirical Orthogonal Functions (EOF) analysis, as suggested by Pham et al. (1998). Hence, AssimilationPreProcessor can also be used to analyse fields without further use in data assimilation.

If state is composed of several properties with different magnitudes then EOFs can be dominated by the properties with larger magnitude, causing a lack of state variability representation. In this multivariate analysis it is convenient to normalize the states previous to covariance calculation. This is accomplished considering a normalization factor for each property (multiplied to the state value) calculated as the inverse of the average standard deviation of all variables of this property. This is a similar approach to the one of Hoteit (2001), which uses the inverse of the square root of the average of variances of each property.

To reduce the computational costs involved in the EOF analysis the covariance matrix is calculated in the space of the time instants instead of the space of state variables. This implies a lower dimension of the matrix if the number of state variables is larger than the number of time instants considered for covariance calculation.

The EOF Analysis is performed by eigenvalues and eigenvectors decomposition of the covariance matrix, using the power method. The obtained eigenvectors are multiplied by the square root of the respective eigenvalue to obtain the EOF in the time instants space. These are then used together with covariance matrix and eigenvalue to obtain the EOF expansion coefficient.

Finally, the EOF in the time instants space is translated to the state variables space.

Optionally, can be made the reconstruction of the state variables, at the several time instants, using the calculated EOFs and the average state.

*References:*

Hoteit, I., 2001, *Filtres de Kalman Reduits et Efficaces pour l'Assimilation de Données en Oceanographie*, Ph.D. thesis, Université de Joseph Fourrier - Grenoble I.

Pham, D., J. Verron and M. Roubaud, 1998, "A singular evolutive extended Kalman filter for data assimilation in oceanography", *Journal of Marine Systems*, 16, pp. 323-340.

**Typical use:**

Perform EOF analysis of spatial property fields usable for sequential data assimilation in MOHID Water.

**Data input requirements:**

One or more HDF5 files containing spatial fields of hydrodynamic and water properties.

**Ouput:**

One HDF5 file with results of EOF analysis (EOFs, eigenvalues, inertia) and the covariance structure, together with general statistics calculated in the processing (e.g. average, standard deviation).

One TimeSeries file with the expansion coeficient is produced for each EOF calculated.

When the state reconstruction is commanded a HDF5 file is produced with the fields used for EOF analysis reconstructed using only the calculated EOFs.

## Input file

The name of the input file must be provided in the [nomfich.dat file](../io/nomfich-dat-file.md) in use.

```
(block for each HDF5 file containing data to be analysed; may be several, block order is
irrelevant)
<BeginHDF5File>
NAME                    : ... (path/name of HDF5 file with data to extract time series)
<EndHDF5File>

START_TIME              : ... (start time for data analysis: yyyy mm dd hh mm ss)
END_TIME                : ... (end time for data analysis: yyyy mm dd hh mm ss)

HDF5_MAP_ITEM           : ... (time independent map item name in HDF5 file)

3D_HDF5                 : 0/1 (0 = 2D HDF5 files, 1 = 3D HDF5 files; 0 = default)

METHOD                  : 1/... (sequential data assimilation method to which the covariance
                                 structure is intended; 1=SEEK; if the application is intended
                                 for EOF analysis then 1 should be chosen)

NORMALIZATION           : 0/1 (0 = no normalization of state; 1 = normalization of state; 0 =
                               default)

DECIMATION_FACTOR       : ... (interval in number of instants not to be considered in the
                               decimation process; 0 = no decimation; 0 = default)

STATE_RECONSTRUCTION    : 0/1 (0 = no reconstruction of state with EOFs estimated; 1 =
                               reconstruction of state with EOFs estimated; 0 = default)

MAX_BUFFER_SIZE         : ... (size in bytes of buffer for expansion coefficient time series;
                               100000 = default)

OUTPUTFILENAME          : ... (path/name of output file (HDF5) with analysis results)

(if METHOD : 1:)
STATECOV_RANK           : ... (dimension of the covariance subspace; number of EOFs estimated
                               if METHOD : 1)

(if STATE_RECONSTRUCTION : 1:)
STATE_OUTPUTFILENAME    : ... (path/name of output file (HDF5) with reconstructed state)

(block for each property/parameter which belongs to state definition, may be several)
<beginproperty>
NAME                    : ... (property name, according with MOHID V4)

UNITS                   : ... (property units)

DIMENSION               : 2D/3D (3D = default)

HDF_GROUP               : ... (complete path in HDF5 file to parameter data, according with
                               MOHID V4)

STATE_WINDOW            : ... ... ... ... ... ... (state spatial window limits: ILB (i lower-
                                                   left cell), JLB (j lower-left cell), IUB
                                                   (i upper-right cell), JUB (j upper-right
                                                   cell), KLB (lower layer), KUB (upper layer))

TYPE_ZUV                : Z/U/V (type of horizontal grid: Z = cell center, U = cell U faces,
                                 V = cell V faces; Z = default)

(if TYPE_ZUV : Z and (NAME : velocity U or NAME : velocity V)
CONVERT_TO_FACES        : 0/1 (0 = not convert to horizontal faces grid; 1 = convert to faces
                               grid; 0 = default)
<endproperty>
```

```
Remarks:

- only METHOD : 1 (SEEK) is currently implemented; this should be the option to use in
  preprocessing for MOHID Water applications using SEEK and SFEK filters as sequential data
  assimilation methods;

- if the application is used to perform EOF analysis then METHOD : 1 must be used;

- in DECIMATION_FACTOR the value should be inserted according with the following example: if 1
  out of 6 time instants is to be considered then DECIMATION_FACTOR : 5 must be used;

- if STATE_RECONSTRUCTION : 1 then state is reconstructed for the time instants considered in
  the analysis using the number of EOFs defined in STATECOV_RANK;

- the order in which the property blocks appear in the input file is considered to the order
  in which the properties appear in state definition;

- in TYPE_ZUV should be indicated the type of grid in which the property fields are available
  in input HDF5 files, which is usually Z; if input HDF5 files containing hydrodynamic
  properties are generated in MOHID Water simulations with the OUTPUT_FACES keyword selected
  in Hydrodynamic_#.dat file then U or V types can be selected;

- if EOFs are intended to be produced to use in SEEK/SFEK filters initialization and velocity
  U/velocity V properties are considered for state definition then the properties must be
  specified in state definition as having type U or V grid: this can be made indicating in
  TYPE_ZUV if velocities are present in input HDF5 files in faces grids or if otherwise by
  selecting CONVERT_TO_FACES : 1.
```

## Sample

```
<BeginHDF5File>
NAME                    : Hydrodynamic_3.hdf5
<EndHDF5File>

<BeginHDF5File>
NAME                    : Hydrodynamic_4.hdf5
<EndHDF5File>

<BeginHDF5File>
NAME                    : Hydrodynamic_5.hdf5
<EndHDF5File>

<BeginHDF5File>
NAME                    : Hydrodynamic_6.hdf5
<EndHDF5File>

START_TIME              : 1972 11 1 0 0 0
END_TIME                : 1972 11 10 0 0 0

HDF5_MAP_ITEM           : WaterPoints3D

3D_HDF5                 : 1

METHOD                  : 1

STATECOV_RANK           : 50

OUTPUTFILENAME          : InitialCovWrongModel_Nov72Dez72_Sampling1h.hdf5

NORMALIZATION           : 1

DECIMATION_FACTOR       : 0

STATE_RECONSTRUCTION    : 1

STATE_OUTPUTFILENAME    : RecStateWrongModel_Nov72Dez72_Sampling1h.hdf5

<beginproperty>
NAME                    : velocity U
UNITS                   : m/s
DIMENSION               : 3D

HDF_GROUP               : /Results/FacesVelocityU

STATE_WINDOW            : 1 162 1 162 1 1

TYPE_ZUV                : U
CONVERT_TO_FACES        : 0
<endproperty>

<beginproperty>
NAME                    : velocity V
UNITS                   : m/s
DIMENSION               : 3D

HDF_GROUP               : /Results/FacesVelocityV

STATE_WINDOW            : 1 162 1 162 1 1

TYPE_ZUV                : V
CONVERT_TO_FACES        : 0
<endproperty>

<beginproperty>
NAME                    : water level
UNITS                   : m
DIMENSION               : 2D

HDF_GROUP               : /Results/water level

STATE_WINDOW            : 1 162 1 162

TYPE_ZUV                : Z
<endproperty>
```
