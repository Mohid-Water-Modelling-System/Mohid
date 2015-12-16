export CC=icc
export FC=ifort
export CXX=icpc
export CFLAGS='-O3 -xHost -ip'
export CXXFLAGS='-O3 -xHost -ip'
export FFLAGS='-O3 -xHost -ip'
#export LDFLAGS="-Wl,--rpath,/usr/local/lib"
tar -xvf hdf5-1.8.15-patch1.tar
cd hdf5-1.8.15-patch1
./configure --prefix=/home/jauch/Development/Libraries --enable-fortran \
--with-szlib=/home/jauch/Development/Libraries --with-zlib=/home/jauch/Development/Libraries
make
make check
make install

