# *NIX platforms

This is a tutorial on how to install MOHID in linux and derived other platforms. To complement this tutorial, you can also check out some [video tutorials](../web/codeplex.md) in the [Codeplex](../web/codeplex.md) article.

## Installing

### Prerequisites

1. gcc, g++ and make-tools for linux
2. Subversion
3. Intel Fortran 9 or greater

*NOTE*: this article was written assuming the fortran compiler is *Intel*. If this is not the case then change where there is written "ifort" with the fortran compiler in use (ex: openf95, pgf90, g90, etc...). In particular, the MOHID code is found to build well with the [open64](../../development/build/open64.md) fortran.

#### HDF5 libraries

To compile Mohid it is mandatory to have the HDF5 [

$$
1
$$

](http://hdf.ncsa.uiuc.edu/HDF5/) libraries compiled.

##### Extract from Codeplex

HDF5 libraries are available at Codeplex.

##### Compile from source

In HDF5's web page it is possible to obtain both binary distribuitions and the source code. Most of the times it is necessary to build from source; so download the source code, e. g. hdf5-1.6.5.tar.gz;

```
> tar -zxvf hdf5-1.6.5.tar.gz</code>
> cd hdf5-1.6.5</code>
> F9X=ifort ./configure --enable-fortran (--prefix=$HDF5ROOTDIR --enable-parallel)</code>
> make</code>
> make check</code>
> make install</code>
```

Detailed information may be found in file ./release\_docs/INSTALL.

#### Zlib libraries

To use the HDF5 libraries it is mandatory to have installed the zlib libraries. So install the zlib package for your linux distro. Also, keep a record of where the libz.a library is installed (it might come in handy later on).

##### Extract from Codeplex

Zlib libraries are available at Codeplex.

#### Netcdf libraries

The netcdf libraries are optional since they are not required to run [MohidWater](../../reference/mohidwater.md). However they are useful for other modules of the Mohid Suite such as [ConvertToHDF5](../../tools/converttohdf5.md). You can get them from the [official site](http://www.unidata.ucar.edu/software/netcdf/). Here's the [report](netcdf.md) on how to build them for windows.

### Checking out the MOHID code

Use subversion to checkout the latest version of the code stored at [codeplex](../web/codeplex.md).

```
> mkdir [mohidpath]
> svn co https://mohid.svn.codeplex.com/svn [mohidpath]
```

Or update to the latest version of the code

```
> svn update
```

### Building MohidWater

Make sure that you have Intel Fortran installed (any version later than 9.x should do). Follow instructions in the README file that is in the Codeplex distribution, under the Linux folder (you are required to create the Editme.mk file from a template and edit your hdf5 and Zlib libraries paths).

To build the Linux makefile solution, just type

```
> make Nix
```

To clean the distribution, type:

```
> make Nix.clean
```

Make sure you also check our [videos](../web/codeplex.md).

## Running Mohid projects in linux

### General guidelines that change from windows

1. There are no graphical interfaces available. No GUI, no PostProcessor and no GIS.
2. Linux is case-sensitive. This is very important.
3. Linux reads paths to files and folders using the slash '/' instead of the backslash '\\'.

Required files for mohidwater to run (preserve letters capitalization!). You must copy (and edit) the files manually:

```
exe/nomfich.dat
exe/tree.dat
```

## TroubleShooting

**Q: the makefile keeps giving me errors :(**

A:Here's a list of what could be causing these errors:

1. Open the makefiles and check that fields are properly filled.
2. Make sure that **mod** and **bin** directories exist in the module's folder. Otherwise simply make them with the mkdir command.
3. Make sure that "ifort" is a valid command. If it isn't then you probably didn't loaded its environment variables. Type "source /opt/intel/fc/9.0/bin/ifortvars.sh" for that.
4. Make sure that the intel fortran compiler was properly installed.

**Q: the makefile complains when compiling ModuleGOTM.f90**

A: Probably you have to rename

GOTM**V**ariables\_in.**f90** ---> GOTM**v**ariables\_in.**f90**

GOTM**V**ariables\_out.**f90** ---> GOTM**v**ariables\_out.**f90**.

This error shows up in linux because it's case-sensitive whereas windows isn't. Also, make sure that the extension is written in lower case.

**Q: By the time i'm linking the final MohidWater binary I get some error about** unreferenced functions "compress2" and others in HDF5 files...

A: Probably you need to specify to link against the libz.a library. Identify the exactly the location of libz.a (ex: in my system it can be found at /opt/zlib/lib/libz.a) then add edit accordingly the *Solutions/Linux/Editme.mk* file.

## See also

- [Makefile](../../building/linux/makefile.md)
- [Mohid Makefile](../../building/linux/mohid-makefile.md)
- [ifort](../../building/linux/ifort.md)
- [open64](../../development/build/open64.md)
