#PLEASE EDIT UNCOMMENTED LINES
#BEFORE RUNNING THE MOHID MAKEFILE

#Edit the prefix of the installed binaries
# > make nix.install
export VER = x64_single

#Where do you want to install the binary files?
# > make nix.install
export DESTDIR = ~/

#Where are the hdf5 libraries (with --enable-fortran) in your system?
# > sudo updatedb; locate hdf5.mod
export HDF5INC = ../../../../../linux/hdf5/include
# > sudo updatedb; locate libhdf5
export HDF5LIB = ../../../../../linux/hdf5/lib

#Where is the libz.a in your system?
export ZLIBINC = ../../../../../linux/zlib/lib

#Activate the extra modules that require the netcdf libraries in the ConvertToHdf5 tool
#Two valid options
# true
# false
#Un-comment your choice
#export IS_NETCDF = true
export IS_NETCDF = false

#export IS_PROJ4F = true
export IS_PROJ4F = false

#Uncomment the desired pre-processing definitiions
#_NO_NETCDF is activated by default.
FPP_DEFINES :=

ifeq ($(IS_NETCDF),false)
    FPP_DEFINES := ${FPP_DEFINES} -D_NO_NETCDF
endif

ifeq ($(IS_PROJ4F),true)
    FPP_DEFINES := ${FPP_DEFINES} -D_USE_PROJ4
endif


FPP_DEFINES := ${FPP_DEFINES} -D_USE_MPI


export FPP_DEFINES
