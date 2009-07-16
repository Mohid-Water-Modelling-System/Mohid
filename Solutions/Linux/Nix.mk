SHELL = /bin/sh

###### NIX Exported variables #########
export CP  = sudo cp
export DEL = rm
export O = o
export F = F90
export MOD = mod
export CC= ifort
export CCFLAGS  = -c -fpp -warn all -nologo -convert big_endian -D_USE_NIX# Debug: -g Profiling: -p
export LFLAGS   = -fpp -nologo -warn all -i-static -convert big_endian -D_USE_NIX# Profiling: -p
export LLFLAGS  =
export AR = ar rc
export SRCBASE1 = ../../../Shared/MOHID.Base.1
export SRCBASE2 = ../../../Shared/MOHID.Base.2
export SRCWATER = ../../../Modulus.Software/MOHID.Water
export BASE1INC = ../Mohid_Base_1
export BASE1 = Mohid_Base_1.lib
export BASE2INC = ../Mohid_Base_2
export BASE2 = Mohid_Base_2.lib
export WATER = MohidWater
export LHDF5FORTRAN = libhdf5_fortran.a
export LHDF5 = libhdf5.a
export LHDF5HL = libhdf5_hl.a
export ZLIB = libz.a
export BASELIBS = \
       $(BASE1INC)/$(BASE1) \
       $(BASE2INC)/$(BASE2) \
       $(HDF5)/$(LHDF5FORTRAN) \
       $(HDF5)/$(LHDF5) \
       $(HDF5)/$(LHDF5HL) \
       $(ZLIBINC)/$(ZLIB)
export LNETCDF  = netcdf.a
export NETCDFLIBS := \
                     $(BASELIBS) \
                     $(NETCDFINC)/$(LNETCDF)
#### End of exported variables ####

include Editme.mk
include Files.mk
include Actions.mk
include Dependencies.mk
