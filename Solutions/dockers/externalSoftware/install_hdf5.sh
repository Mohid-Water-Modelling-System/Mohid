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

echo "###### $zlib #########" >> ~/.bashrc
echo "ZLIB=$DIRINSTALL/$zlib" >> ~/.bashrc
echo 'export PATH=$PATH:$ZLIB/lib:$ZLIB/include' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ZLIB/lib' >> ~/.bashrc



    echo
    echo " #### Install $hdf5 ####"
    cd /externalSoftware
    tar -xf "$hdf5.tar.gz"
    cd $hdf5 || exit 1
    CC=$CC FC=$FC  ./configure --with-zlib=$DIRINSTALL/$zlib --prefix=$DIRINSTALL/$hdf5 --enable-fortran --enable-fortran2003 || exit 1
    make || exit 1
    #make check || exit 1
    $sudo make install || exit 1
    echo "####### $hdf5 #######" >> ~/.bashrc
    echo "HDF5=$DIRINSTALL/$hdf5" >> ~/.bashrc
    echo 'export PATH=$PATH:$HDF5/bin:$HDF5/lib:$HDF5/include' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5/lib' >> ~/.bashrc
    echo >> ~/.bashrc
