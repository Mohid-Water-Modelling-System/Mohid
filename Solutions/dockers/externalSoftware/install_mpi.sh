#!/usr/bin/bash
export CXX=icpc
export FC=ifort
export F77=ifort

## intel flags compiler
export CXXFLAGS='-O0 -ip'
export FCFLAGS='-O0 -ip'
export FFLAGS='-O0 -ip'

source /opt/intel/oneapi/setvars.sh

DIRINSTALL=$HOME/apps_intel
zlib='zlib-1.2.11'
hdf5='hdf5-1.8.17'
netcdf='netcdf-4.4.1.1'
netcdff='netcdf-fortran-4.4.4'
mpiv="mpich-3.2"

## trata do mpi 
cd /externalSoftware
tar -xf mpich-3.2.tar.gz
cd $mpiv || exit 1
CC=$CC FC=$FC ./configure --prefix=$DIRINSTALL/$mpiv --enable-fast=O0 --disable-error-checking --without-timing --without-mpit-pvars --with-device=ch4:ofi --disable-cxx || exit 1
make || exit 1
make install || exit 1

echo "########## $netcdf #######" >> ~/.bashrc
echo "export NETCDF=$DIRINSTALL/${netcdf}" >> ~/.bashrc
echo 'export PATH=$PATH:$NETCDF/bin:$NETCDF/lib:$NETCDF/include' >> ~/.bashrc
echo >> ~/.bashrc
echo 'export NETCDF_ROOT=$NETCDF' >> ~/.bashrc
echo 'export NETCDF4_ROOT=$NETCDF' >> ~/.bashrc
echo 'export NETCDF_LIB=$NETCDF/lib' >> ~/.bashrc
echo 'export NETCDF_INC=$NETCDF/include' >> ~/.bashrc
echo >> ~/.bashrc
echo 'export NETCDF_GF_ROOT=$NETCDF' >> ~/.bashrc
echo 'export NETCDF4_GF_ROOT=$NETCDF' >> ~/.bashrc
echo 'export NETCDF_GF_LIB=$NETCDF/lib' >> ~/.bashrc
echo 'export NETCDF_GF_INC=$NETCDF/include' >> ~/.bashrc
echo >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_LIB' >> ~/.bashrc
echo >> ~/.bashrc
echo 'export CPPFLAGS="$CPPFLAGS -I$NETCDF_INC"' >> ~/.bashrc
echo 'export LDFLAGS="$LDFLAGS -L$NETCDF_LIB"' >> ~/.bashrc
echo >> ~/.bashrc

## trata do netcdf Fortran 
cd ..
tar -xf "$netcdff.tar.gz"
cd $netcdff
export CPPFLAGS="-I$DIRINSTALL/$hdf5/include -I$DIRINSTALL/$zlib/include"
export LDFLAGS="-L$DIRINSTALL/$hdf5/lib -L$DIRINSTALL/$zlib/lib"
export NETCDF=$DIRINSTALL/${netcdf}
export PATH=$PATH:$NETCDF/bin:$NETCDF/lib:$NETCDF/include
export NETCDF_ROOT=$NETCDF
export NETCDF4_ROOT=$NETCDF
export NETCDF_LIB=$NETCDF/lib
export NETCDF_INC=$NETCDF/include
export NETCDF_GF_ROOT=$NETCDF
export NETCDF4_GF_ROOT=$NETCDF
export NETCDF_GF_LIB=$NETCDF/lib
export NETCDF_GF_INC=$NETCDF/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_LIB
export CPPFLAGS="$CPPFLAGS -I$NETCDF_INC"
export LDFLAGS="$LDFLAGS -L$NETCDF_LIB"
CC=$CC FC=$FC ./configure --prefix=$DIRINSTALL/$netcdf
make
make install
