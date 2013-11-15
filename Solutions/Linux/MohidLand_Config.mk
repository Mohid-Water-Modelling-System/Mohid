#PLEASE EDIT UNCOMMENTED LINES
#BEFORE RUNNING THE MOHID MAKEFILE

#Edit the prefix of the installed binaries
# > make nix.install
export VER = r8_omp_x64

#Where do you want to install the binary files?
# > make nix.install
export DESTDIR = /usr/bin/mohid

#Where are your intel libraries?
export INTEL_PATH = /opt/intel/lib/intel64/

#Where are the hdf5 libraries (with --enable-fortran) in your system?
# > sudo updatedb; locate hdf5.mod
export HDF5MODINC = /home/user/Work/Mohid_External_Libs/include
# > sudo updatedb; locate libhdf5
export HDF5LIBINC = /home/user/Work/Mohid_External_Libs/lib

#Where is the libz.a in your system?
export ZLIBINC = /home/user/Work/Mohid_External_Libs/lib
export SZLIBINC = /home/user/Work/Mohid_External_Libs/lib

#Find CODEPLEX Respository revision and if there is local changes
LOCAL_CHANGES := $(shell ./test_local_changes.sh)
CODEPLEX_VERSION := $(shell ./get_last_version.sh)
COMPILING_DATE := '$(shell date)'

#Activate the extra modules that require the netcdf libraries in the 
export IS_NETCDF = false
export IS_PROJ4F = false

#Uncomment the desired pre-processing definitiions
#_NO_NETCDF is activated by default.
FPP_DEFINES := 

#FPP_DEFINES := ${FPP_DEFINES} -D_INCREASE_MAXINSTANCES
#FPP_DEFINES := ${FPP_DEFINES} -D_SHORT_LINE_LENGTH
#FPP_DEFINES := ${FPP_DEFINES} -D_LONG_LINE_LENGTH
#FPP_DEFINES := ${FPP_DEFINES} -D_EXTRA_LONG_LINE_LENGTH
#FPP_DEFINES := ${FPP_DEFINES} -D_EXTRA_SHORT_LINE_LENGTH
#FPP_DEFINES := ${FPP_DEFINES} -D_USE_MPI
#FPP_DEFINES := ${FPP_DEFINES} -D_LAGRANGIAN_
#FPP_DEFINES := ${FPP_DEFINES} -D_WAVES_
#FPP_DEFINES := ${FPP_DEFINES} -DOVERLAP
#FPP_DEFINES := ${FPP_DEFINES} -D_SEDIMENT_
#FPP_DEFINES := ${FPP_DEFINES} -D_AIR_
#FPP_DEFINES := ${FPP_DEFINES} -D_LAGRANGIAN_GLOBAL_
#FPP_DEFINES := ${FPP_DEFINES} -D_PHREEQC_      #Not available yet
#FPP_DEFINES := ${FPP_DEFINES} -D_PHREEQC_X64_  #Not available yet
FPP_DEFINES := ${FPP_DEFINES} -D_COMPILED_BY=\"MohidLinux\"
FPP_DEFINES := ${FPP_DEFINES} -D_COMPILING_DATE=\"${COMPILING_DATE}\"
FPP_DEFINES := ${FPP_DEFINES} -D_CODEPLEX_VERSION=\"${CODEPLEX_VERSION}${LOCAL_CHANGES}\"
FPP_DEFINES := ${FPP_DEFINES} -D_NO_NETCDF
FPP_DEFINES := ${FPP_DEFINES} -D_STACK_LIMITS_
FPP_DEFINES := ${FPP_DEFINES} -D_BIG_LINE_LENGTH

export FPP_DEFINES

#rpath paths for intel, hdf5, szlib and zlib libraries


