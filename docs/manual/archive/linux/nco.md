# NCO

They are [netcdf](netcdf.md) files averagers, concateners and renamer operators that act on a netcdf dataset [hyperslab](../../reference/README.md#hyperslab#hyperslab) or an ensemble of netcdf datasets [hyperslabs](../../reference/README.md#hyperslab#hyperslab).

## Installing

### Linux

#### From distros repository

```
> sudo yum install nco
```

or

```
> sudo apt-get install nco
```

#### Manually (recommended)

A manual installation allows to add opendap support besides hdf5 support.

Prerequisites:

- hdf5 libraries with fortran support,
- nc-dap libraries,

```
> wget nco
> tar -xvf nco.tar.gz
> cd nco
> ./configure --enable-hdf5 --enable-fortran --enable-dods
```

## Examples

### ncrename

Rename the variable p to pressure and t to temperature in netCDF ‘in.nc’. In this case p must exist in the input file (or ncrename will abort), but the presence of t is optional:

```
> ncrename -v p,pressure -v .t,temperature in.nc
```

Rename the attribute long\_name to largo\_nombre in the variable u, and no other variables in netCDF ‘in.nc’.

```
> ncrename -a u:long_name,largo_nombre in.nc
```

### ncrcat concatenate

Concatenates along the record dimension (usually the time dimension)

```
> ncrcat month[1].nc month[2].nc month[3].nc month[1-3].nc
```

### ncks extracts

Extract variables time and pressure from netCDF ‘in.nc’. If ‘out.nc’ does not exist it will be created. Otherwise the you will be prompted whether to append to or to overwrite ‘out.nc’:

```
> ncks -v time,pressure in.nc out.nc
```

```
> ncks -C -v time,pressure in.nc out.nc
```

The first version of the command creates an ‘out.nc’ which contains time, pressure, and any coordinate variables associated with pressure. The ‘out.nc’ from the second version is guaranteed to contain only two variables time and pressure.

### ncks append

Appends variables from input file into output file. Will not work if dimension variables are different or have different sizes.

```
> ncks -A fl_1.nc fl_2.nc
```

### nca averages

The following command performs the time average over an ensemble of files (they don't need to be the same size):

```
> ncra file1.nc file2.nc ... filen.nc mean.nc
```

## Applied examples

### Renaming variables names

Renames the Vel.\_X variables to Vel\_X. This is because the "." character isn't a very good idea to be present in a netcdf variable name, or in any variable at all. The matlab scripts don't read such variables.

```
>  for i in `find . | grep 'Portugal_Hydrodynamic'`;
>  do ncrename -v .Vel._X,Vel_X -v .Vel._Y,Vel_Y -v .Vel._Z,Vel_Z -v .Vel._Z,Vel_Z $i;
>  done;
```

### Biscay exercise. Formatting data correctly

G. Reffray helped us a lot to get the correct netcdf files format for the Biscay intercomparison. Here are the lines he provided:

```
> ncpdq -U -h -O -a lat,lon MOHID_GRIDT_12.nc MOHID_GRIDT_12.nc
> ncpdq -h -U -F -O -d deptht,3,43 MOHID_GRIDT_12.nc MOHID_GRIDT_12.nc
> ncpdq -h -O -U -a -deptht MOHID_GRIDT_12.nc MOHID_GRIDT_12.nc
> ncks -h -O -v temperature,salinity MOHID_GRIDT_12.nc MOHID_GRIDT_12.nc
```

The first line transposes the lat and lon dimensions. The second line trims the deptht dimension. The third line inverts the deptht dimension, and finally, the last line, excludes all variables, except the temperature and the salinity.

### PCOMS subset netcdf extraction and gluing

L. Fernandes produces daily netcdf datasets for the portuguese coast available on a thredds catalog [here](http://data.mohid.com:8080/thredds/catalog/MohidAll/pcoms/catalog.html). Our goal is to extract a subset and glue the files together.

```
> ncks
```

### PCOMS merging

```
> time for i in `ls -1 | grep PCOMSV1_T_S`;
for> do j=${i/_T_S/_U_V_SSH}; k=${i/_T_S/};
for> echo merging into $k ...;
for> ncks -A $i $j;
for> mv $j $k;
for> rm -f $i;
for> echo [OK];
for> done
```

### PCOMS renaming a variable in all files

Rename variable *vm* into CF compliant variable name *sea\_water\_speed*:

```
> time for i in `ls -1 | grep U_V`; do sudo ncrename -v vm,sea_water_speed $i; echo . ; done
```

Rename attribute long\_name into CF compliant attribute

```
> time for i in `ls -1 | grep U_V`;
for> do sudo ncatted -O -a standard_name,vm,m,c,"sea_water_speed" $i;
for> echo .;
for> done
```

## External references

- [NCO homepage examples](http://jisao.washington.edu/data_sets/nco/)
