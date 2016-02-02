#tar -xvf hdf5-1.8.15-patch1.tar
cd mpich-3.2
./configure CC=icc F77=ifort FC=ifort CXX=icpc --prefix=~/Libraries --enable-fast=all,O3 \
--enable-shared --with-pm=hydra
make
make install

