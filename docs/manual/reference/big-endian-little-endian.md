# Big-endian little-endian

The Intel Fortran Compiler can write unformatted sequential files in big-endian format and also can read files produced in big-endian format by using the little-endian-to-big-endian conversion feature.

Both on IA-32-based processors and on Itanium®-based processors, Intel Fortran handles internal data in little-endian format. The little-endian-to-big-endian conversion feature is intended for Fortran unformatted input/output operations in unformatted sequential files. The feature enables:

In order to use the little-endian-to-big-endian conversion feature, specify the numbers of the units to be used for conversion purposes by setting the F\_UFMTENDIAN environment variable. Then, the READ/WRITE statements that use these unit numbers, will perform relevant conversions. Other READ/WRITE statements will work in the usual way.

In the general case, the variable consists of two parts divided by a semicolon. No spaces are allowed inside the F\_UFMTENDIAN value. The variable has the following syntax:

```
F_UFMTENDIAN=MODE | [MODE;] EXCEPTION
```

where:

```
MODE = big | little
EXCEPTION = big:ULIST | little:ULIST | ULIST
ULIST = U | ULIST,U
U = decimal | decimal -decimal
```

```
     MODE defines current format of data, represented in the files; it can be omitted.
     The keyword little means that the data have little endian format and will not be converted. This keyword is a default.
     The keyword big means that the data have big endian format and will be converted. This keyword may be omitted together with the colon.
```

```
     EXCEPTION is intended to define the list of exclusions for MODE; it can be omitted. EXCEPTION keyword (little or big) defines data format in the files that are connected to the units from the EXCEPTION list. This value overrides MODE value for the units listed.
```

```
     Each list member U is a simple unit number or a number of units. The number of list members is limited to 64.
     decimal is a non-negative decimal number less than 232.
```

Converted data should have basic data types, or arrays of basic data types. Derived data types are disabled.

Command lines for variable setting with different shells:
