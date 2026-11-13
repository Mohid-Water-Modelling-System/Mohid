# HDF5

HDF5 is a general purpose library and file format for storing scientific data. It's the standard input/output format of spatial and/or temporal data sets in MOHID.

## HDF5 libraries

### Windows

The windows binaries are served directly from our [codeplex](../archive/web/codeplex.md) repository, so there is no need to install them.

### Linux

To compile Mohid it is mandatory to have the HDF5 [

$$
1
$$

](http://hdf.ncsa.uiuc.edu/HDF5/) libraries compiled. In HDF5's web page it is possible to obtain both binary distribuitions and the source code. Most of the times it is necessary to build from source which takse the following steps:

1. download the source code, e. g. hdf5-1.6.5.tar.gz;
2. `> tar -zxvf hdf5-1.6.5.tar.gz`
3. `> cd hdf5-1.6.5`
4. `> F9X=ifort ./configure --enable-fortran (--prefix=$HDF5ROOTDIR --enable-parallel)`
5. `> make`
6. `> make check`
7. `> make install`

Libraries may be found in ./hdf5/lib/ and, typically, in /opt/hdf5/hdf5/lib. Detailed information may be found in file ./release\_docs/INSTALL.

#### Zlib libraries

To use the HDF5 libraries it is mandatory to have installed the zlib libraries. So install the zlib package for your linux distro. Also, keep a record of where the libz.a library is installed (it might come in handy later on).#

### Back to Mohid developers page

- [Developers](../reference/developers.md)

## References

HDF5 Homepage - [

$$
2
$$

](http://hdf.ncsa.uiuc.edu/HDF5/)
