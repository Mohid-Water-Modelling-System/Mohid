export CC=icc
export CFLAGS='-O3 -xHost -ip'
tar -zxvf zlib-1.2.8.tar.gz
cd zlib-1.2.8 
./configure --prefix=~/Libraries
make
make check
make install
