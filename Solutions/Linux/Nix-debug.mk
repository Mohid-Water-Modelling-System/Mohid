SHELL = /bin/sh

#------User configuration file---------------

# guillaume

include Editme.mk
francisco
#-- NIX platform specific global variables --

export CP  = sudo cp
export DEL = rm
export O = o
export F = F90
export MOD = mod
export CC= ifort
export CCFLAGS  = -c -fpp -warn all -nologo -convert big_endian -fpe0 -traceback -g -p -D_USE_NIX $(FPP_DEFINES)# Debug: -g Profiling: -p
export LFLAGS   = -fpp -nologo -warn all -i-static -convert big_endian -p -D_USE_NIX# Profiling: -p
export LLFLAGS  =
export MKFLAGS =
export AR = ar rc
export SUFFLIB = .lib
export SUFFPROG = 

# MohidBase1
export SRCBASE1 = ../../../Software/MOHIDBase1
export BASE1INC = ../Mohid_Base_1
export BASE1 = Mohid_Base_1$(SUFFLIB)

# MohidBase2
export SRCBASE2 = ../../../Software/MOHIDBase2
export BASE2INC = ../Mohid_Base_2
export BASE2 = Mohid_Base_2$(SUFFLIB)

# MohidWater
export SRCWATER = ../../../Software/MOHIDWater
export WATER = MohidWater$(SUFFPROG)

# MohidLand
export SRCLAND = ../../../Software/MOHIDLand
export LAND = MohidLand$(SUFFPROG)

export SRCCONVERT = ../../../Software/ConvertToHDF5
export CONVERT = ConvertToHDF5$(SUFFPROG)

# HDF5 lib
export LHDF5FORTRAN = libhdf5_fortran.a
export LHDF5 = libhdf5.a
export LHDF5HL = libhdf5_hl.a

# Z lib
#export ZLIB = libz.a
export ZLIB = libz.so


# All libs folders
export BASELIBS := \
       $(BASE1INC)/$(BASE1) \
       $(BASE2INC)/$(BASE2) \
       $(HDF5)/$(LHDF5FORTRAN) \
       $(HDF5)/$(LHDF5) \
       $(HDF5)/$(LHDF5HL) \
       $(ZLIBINC)/$(ZLIB)

# Netcdf lib
export LNETCDF  = libnetcdf.a

# libfproj4 static lib
export LPROJ4F = libfproj4.a

# All libs folders (including netcdf)
ifeq ($(IS_NETCDF),true)
    export BASELIBS := \
            $(BASELIBS) \
            $(NETCDFINC)/$(LNETCDF)
endif

ifeq ($(IS_PROJ4F),true)
    export BASELIBS := \
            $(BASELIBS) \
            $(PROJ4FLIB)/$(LPROJ4F)
endif

#------Files and modules lists------

METAFILES = \
        README \
        Editme_template.smk \
        Nix.smk
MODULES = \
          Makefiles \
          MohidLand \
          MohidWater \
          Mohid_Base_2 \
          Mohid_Base_1 \
          ConvertToHDF5

#------Makefile rules---------------

include Makefiles/Makemodules.mk

#------Modules dependencies----------

Mohid_Base_2.all : Mohid_Base_1.all
MohidLand.all : Mohid_Base_2.all
MohidWater.all : Mohid_Base_2.all
ConvertToHDF5.all : Mohid_Base_2.all

