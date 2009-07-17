#PLEASE EDIT UNCOMMENTED LINES
#BEFORE RUNNING THE MOHID MAKEFILE

#Edit the prefix of the installed binaries
# > make nix.install
export VER = x64_single

#Where do you want to install the binary files?
# > make nix.install
export DESTDIR = /usr/bin/mohid

#Where are the hdf5 libraries (with --enable-fortran) in your system?
export HDF5 = /home/Projects/hdf5/hdf5-1.6.5/hdf5/lib

#Where is the libz.a in your system?
export ZLIBINC = /usr/lib64

#Where is the libnetcdf.a (with --enable-fortran) in your system?
export NETCDFINC = /usr/lib64
