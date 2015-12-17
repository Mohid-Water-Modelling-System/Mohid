export CC=icc
export CXX=icpc
export F77=ifort
export CFLAGS='-O3 -xHost -ip'
export CXXFLAGS='-O3 -xHost -ip'
export FFLAGS='-O3 -xHost -ip'
tar -zxvf szip-2.1.tar.gz 
cd szip-2.1
./configure --prefix=/home/jauch/Development/Libraries
make
make check
make install
