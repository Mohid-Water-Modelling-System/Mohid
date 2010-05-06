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

#Uncomment the desired pre-processing definitions
#If you don't know what these mean, then leave them commented ;)
FPP_DEFINES := ""
#FPP_DEFINES := ${FPP_DEFINES} -D_INCREASE_MAXINSTANCES
##FPP_DEFINES := ${FPP_DEFINES} -D_SHORT_LINE_LENGTH
##FPP_DEFINES := ${FPP_DEFINES} -D_LONG_LINE_LENGTH
#FPP_DEFINES := ${FPP_DEFINES} -D_EXTRA_LONG_LINE_LENGTH
#FPP_DEFINES := ${FPP_DEFINES} -D_EXTRA_SHORT_LINE_LENGTH
#FPP_DEFINES := ${FPP_DEFINES} -D_USE_MPI
#FPP_DEFINES := ${FPP_DEFINES} -D_LAGRANGIAN_
#FPP_DEFINES := ${FPP_DEFINES} -D_WAVES_
#FPP_DEFINES := ${FPP_DEFINES} -DOVERLAP
#FPP_DEFINES := ${FPP_DEFINES} -D_SEDIMENT_
#FPP_DEFINES := ${FPP_DEFINES} -D_AIR_
#FPP_DEFINES := ${FPP_DEFINES} -D_LAGRANGIAN_GLOBAL_

export FPP_DEFINES

