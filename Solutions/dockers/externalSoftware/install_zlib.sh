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

echo
echo " #### Install $zlib ####"
cd /externalSoftware
tar -xf "$zlib.tar.gz"
cd $zlib
./configure --prefix=$DIRINSTALL/$zlib || exit 1
make || exit 1
$sudo make install || exit 1
echo "###### $zlib #########" >> ~/.bashrc
echo "ZLIB=$DIRINSTALL/$zlib" >> ~/.bashrc
echo 'export PATH=$PATH:$ZLIB/lib:$ZLIB/include' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ZLIB/lib' >> ~/.bashrc
echo >> ~/.bashrc
