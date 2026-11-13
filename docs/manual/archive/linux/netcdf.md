# Netcdf

NetCDF (network Common Data Form) is an interface for array-oriented data access and a library that provides an implementation of the interface. The netCDF library also defines a machine-independent format for representing scientific data. Together, the interface, library, and format support the creation, access, and sharing of scientific data.

## Netcdf library

The netcdf library provides some useful binary tools suchs as [ncdump](../../tools/ncdump.md). It is required to link against [ConvertToHdf5](../../tools/converttohdf5.md) and [Convert2netcdf](../../tools/convert2netcdf.md) tools.

### Netcdf 4

Netcdf 4 allows to read hdf5 files as well and is fully backwards compatible.

#### How to compile in windows

Install the cygwin platform with the gcc, g++, g77 and gfortran packages. Download and unpack netcdf 4 latest releases from the cygwin bash console. The windows ifort doesn't work quite the same as the linux ifort, so we must use gFortran to compil...

```
> CPPFLAGS=-DgFortran ./configure --enable-dll --disable-dap
> make
> make install
```

Unfortunately, the library will not be compatible with the intel fortran compiler ...

### Netcdf 3.6.x

### Netcdf 3.5.1

#### How to compile for intel fortran in windows

This how to will only explain how to compile the static libraries.

- Get the netcdf source distribution version [3.5.1](http://www.unidata.ucar.edu/downloads/netcdf/ftp/old/netcdf-3.5.1.zip).
- Unpack somwhere
- From the netcdf-3.5.1/src/libsrc

```
nmake /f msofts.mak
```

- From the netcdf-3.5.1/src/fortran

```
nmake /f msofts.mak
```

- From the netcdf-3.5.1/src/f90

```
ifort -c  -Gm typeSizes.f90
ifort -c  -Gm netcdf.f90
ar cru ../libsrc/netcdfs.lib netcdf.obj typeSizes.obj
ranlib ../libsrc/netcdfs.lib
```

- Copy the netcdfs.lib and the mod files to IntelLibs and that's it!

## Netcdf tools

### Nc bins

- ncgen
- [ncdump](../../tools/ncdump.md)

### NcBrowse

### [Matlab](matlab.md)

### [NCO](nco.md) NetCdf Operators

### [Perl](perl.md)

To compile the NetCDF module in perl in **linux**:

- 1 Make sure you know where the **netcdf** libraries are installed in your system:

```
> sudo find / | grep libnetcdf.a
> sudo find / | grep netcdf.h
```

- 2 Go to this [site](http://home.badc.rl.ac.uk/iwi/netcdf-perl/) and extract the [module source code](http://home.badc.rl.ac.uk/iwi/netcdf-perl/downloads/netcdf-perl-1.2.3.tar.gz);
- 3 Unpack the module

```
> tar -xvf netcdf-perl-x.tar.gz
> cd netcdf-perl-x/src
> ./configure
> make
> sudo make install
```

- 4 Locate where it was installed

```
> sudo find / | grep NetCDF.pm
```

- 5 Copy the [ncmanipulator file](http://home.badc.rl.ac.uk/iwi/netcdf-perl/downloads/ncmanipulate.pm) in the same directory as *NetCDF.pm*
- 6 That should do it!

*NOTE: Couldn't do it in windows :(*

[Check out](perl.md) a sample perl program to manipulate a netcdf file.

## Problems while compiling netcdf for fortran in windows or in linux

The binaries distributed for windows are currently only for the c and c++ interface. This is a problem since we want the netcdf libraries with the fortran interface.

### What the compilation thread should be

Here's the log of a linux compilation with make:

First we build the c library:

```
make[2]: Entering directory `/home/guillaume/Transferências/netcdf-3.6.1/src/libsrc
'
cc -c -g -O2 -I.. -I.   -DpgiFortran attr.c
cc -c -g -O2 -I.. -I.   -DpgiFortran dim.c
cc -c -g -O2 -I.. -I.   -DpgiFortran error.c
cc -c -g -O2 -I.. -I.   -DpgiFortran -DVERSION=`cat ../VERSION` libvers.c
cc -c -g -O2 -I.. -I.   -DpgiFortran nc.c
cc -c -g -O2 -I.. -I.   -DpgiFortran ncio.c
cc -c -g -O2 -I.. -I.   -DpgiFortran ncx.c
cc -c -g -O2 -I.. -I.   -DpgiFortran putget.c
cc -c -g -O2 -I.. -I.   -DpgiFortran string.c
cc -c -g -O2 -I.. -I.   -DpgiFortran v1hpg.c
cc -c -g -O2 -I.. -I.   -DpgiFortran v2i.c
cc -c -g -O2 -I.. -I.   -DpgiFortran var.c
ar cru libnetcdf.a attr.o dim.o error.o libvers.o nc.o ncio.o ncx.o putget.o string
.o v1hpg.o v2i.o var.o
ranlib libnetcdf.a
make[2]: Leaving directory `/home/guillaume/Transferências/netcdf-3.6.1/src/libsrc'
```

Then we build the C interface with fortran:

```
make[2]: Entering directory `/home/guillaume/Transferências/netcdf-3.6.1/src/fortra
n'
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-attio.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-control.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-dim.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-genatt.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-geninq.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-genvar.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-lib.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-misc.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-v2compat.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-vario.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-var1io.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-varaio.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-varmio.c
cc -c -g -O2 -I.. -I../libsrc   -DpgiFortran fort-varsio.c
ar cru ../libsrc/libnetcdf.a fort-attio.o fort-control.o fort-dim.o fort-genatt.o f
ort-geninq.o fort-genvar.o fort-lib.o fort-misc.o fort-v2compat.o fort-vario.o fort
-var1io.o fort-varaio.o fort-varmio.o fort-varsio.o
ranlib ../libsrc/libnetcdf.a
make[2]: Leaving directory `/home/guillaume/Transferências/netcdf-3.6.1/src/fortran
```

Finally we build the fortran library:

```
make[2]: Entering directory `/home/guillaume/Transferências/netcdf-3.6.1/src/f90'
ifort -c  typeSizes.f90
ifort -c  netcdf.f90
ar cru ../libsrc/libnetcdf.a netcdf.o typeSizes.o
ranlib ../libsrc/libnetcdf.a
make[2]: Leaving directory `/home/guillaume/Transferências/netcdf-3.6.1/src/f90'
```

Follow these steps and you should be fine.

### Problems arose from [name mangling](http://en.wikipedia.org/wiki/Name_decoration)

However different compilers will produce different [name decorations](http://en.wikipedia.org/wiki/Name_decoration) and then [linkers](http://en.wikipedia.org/wiki/Linker) will have a hell of a time trying to resolve the references from different [object files](http://en.wikipedia.org/wiki/Object_file) and/or [libraries](http://en.wikipedia.org/wiki/Library_%28computer_science%29). In netcdf, it is likely for this to be the case.

#### in windows

For example, here's an extract of the external references of ftest.obj after compiling with df (compaq Fortran):

```
df /c /nologo /Zi /Od ftest.for
dumpbin /SYMBOLS ftest.obj
```

```
094 00000000 UNDEF  notype ()    External     | _NCDREN@20
095 00000000 UNDEF  notype ()    External     | _NCVID@16
096 00000000 UNDEF  notype ()    External     | _NCVREN@20
097 00000000 UNDEF  notype ()    External     | _NCAREN@28
098 00000000 UNDEF  notype ()    External     | _NCVDEF@28
099 00000000 UNDEF  notype ()    External     | _NCVGT@24
09A 00000000 UNDEF  notype ()    External     | _NCVGTC@32
09B 00000000 UNDEF  notype ()    External     | _NCVGT1@20
09C 00000000 UNDEF  notype ()    External     | _NCVG1C@24
09D 00000000 UNDEF  notype ()    External     | _NCVPT@24
09E 00000000 UNDEF  notype ()    External     | _NCVPTC@32
09F 00000000 UNDEF  notype ()    External     | _NCVPT1@20
0A0 00000000 UNDEF  notype ()    External     | _NCVP1C@24
```

and here's an extract after compiling ftest.for with [ifort](../../building/linux/ifort.md) (Intel Fortran):

```
ifort /c /nologo /Zi /Od ftest.for
dumpbin /SYMBOLS ftest.obj
```

```
1D8 00000000 UNDEF  notype ()    External     | _NCDREN
1D9 00000000 UNDEF  notype ()    External     | _NCVID
1DA 00000000 UNDEF  notype ()    External     | _NCVREN
1DB 00000000 UNDEF  notype ()    External     | _NCAREN
...
1E5 00000000 UNDEF  notype ()    External     | _NCVDEF
...
205 00000000 UNDEF  notype ()    External     | _NCVGT
206 00000000 UNDEF  notype ()    External     | _NCVGTC
```

Hence the same object will have problems linking to other objects or libraries if the [name mangling](http://en.wikipedia.org/wiki/Name_decoration) convention is diferent! How do you solve that?

Well, let's take a look at the static version of the netcdf library, netcdfs.lib:

```
071 00000DE0 SECT5  notype ()    External     | _NCDREN@20
072 00000ED0 SECT5  notype ()    Static       | _c_ncdren
073 00000000 UNDEF  notype ()    External     | _ncdimrename
074 00000F00 SECT5  notype ()    External     | _NCVINQ@36
075 00000000 UNDEF  notype ()    External     | _c2f_dimids
076 000010B0 SECT5  notype ()    Static       | _c_ncvinq
077 00000000 UNDEF  notype ()    External     | _ncvarinq
078 000010F0 SECT5  notype ()    External     | _NCVPT1@20
079 00000000 UNDEF  notype ()    External     | _f2c_coords
07A 00001150 SECT5  notype ()    Static       | _c_ncvpt1
```

As you can see the name mangling convention is compatible with the df compiler default settings and not with the [ifort](../../building/linux/ifort.md) compiler default settings. Thus it will link correctly with the df-obj whereas it will return unreferenced errors with the ifort-obj.

So, what's the solution? Well, you have to tell the compiler to use the correct [name mangling convention](http://en.wikipedia.org/wiki/Name_decoration):

```
ifort /Gm /c /nologo /Zi /Od ftest.for
```

And if you try link it against the netcds.lib, it'll work! Thanks for the /Gm option defined in [ifort](../../building/linux/ifort.md).

#### in linux

In linux machines, sometimes, such name mangling may occur. Let's check which name convention is used with the [ifort](../../building/linux/ifort.md) compiler:

```
> ifort -c -fpp ModuleMERCATORFormat.F90
> nm ModuleMERCATORFormat.F90 | grep netcdf
                U netcdf_mp_nf90_close_
                U netcdf_mp_nf90_get_att_text_
                U netcdf_mp_nf90_get_var_1d_eightbytereal_
                U netcdf_mp_nf90_get_var_1d_fourbytereal_
                U netcdf_mp_nf90_get_var_2d_fourbytereal_
                U netcdf_mp_nf90_get_var_3d_fourbytereal_
                U netcdf_mp_nf90_get_var_4d_fourbytereal_
                U netcdf_mp_nf90_inq_dimid_
                U netcdf_mp_nf90_inquire_
                U netcdf_mp_nf90_inquire_dimension_
                U netcdf_mp_nf90_inquire_variable_
                U netcdf_mp_nf90_inq_varid_
                U netcdf_mp_nf90_open_
```

Now let's verify what the *libnetcdf.a* library uses the same name decoration:

```
> nm /usr/lib/libnetcdf.a | grep netcdf | grep '_$'
...
0000000000001796 T netcdf_mp_nf90_put_var_text_
000000000000563a T netcdf_mp_nf90_put_var_twobyteint_
0000000000000190 T netcdf_mp_nf90_redef_
0000000000000452 T netcdf_mp_nf90_rename_att_
00000000000003b0 T netcdf_mp_nf90_rename_dim_
000000000000178c T netcdf_mp_nf90_rename_var_
...
```

In this case, it turns out that both the netcdf library and the calling module use the same name decoration, thus there will be no error at link time arose from calling conventions for netcdf functions.

## How to Use Netcdf5 files

The netcdf libraries should be updated in MOHID projects. While that is not done there is a workaround that is to convert netcdf v5 files in v4:

1. Go to [Netcdf libraries](http://www.unidata.ucar.edu/software/netcdf/docs/winbin.html) and download netCDF4 libraries (32 or 64bit)
2. Run the executable to install
3. Go to your folder where it was installed (e.g. C:\\Program Files\\netCDF 4.3.0\\bin) and copy those files to any work folder where to convert the netcf file (e.g D:\\ConvertNetcdf)
4. Copy to the work folder the netcdf file to convert
5. Start a command line and go to the work folder (e.g. cd D:\\ConvertNetcdf)
6. Run in command prompt: nccopy -k classic 

$$
name of final file in v4 (with extension)
$$

 

$$
name of original file in v5 (with extension)
$$


7. Done
8. You can test the final file to check that is v4 run in the same folder command prompt: ncdump -k 

$$
name of final file (with extension)
$$


9. It should report "classic" if it is a v4 file

## External references

- [Netcdf at unidata](http://www.unidata.ucar.edu/software/netcdf/)
- [NCO homepage](http://nco.sourceforge.net/)
- [Matlab interface for netcdf](http://www.marine.csiro.au/sw/matlab-netcdf.html)
