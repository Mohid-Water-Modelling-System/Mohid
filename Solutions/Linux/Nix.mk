SHELL = /bin/sh

#------User configuration file---------------

include Editme.mk

#-- NIX platform specific global variables --

export CP  = sudo cp
export DEL = rm
export Obj = o
export F = F90
export MOD = mod
export CC= ifort
export CCFLAGS  = -c -r8 -inline-level=0 -fpp -warn all -nologo -convert big_endian -fpe0 -D_USE_NIX -traceback -mcmodel=large -heap-arrays 64 -openmp -check bounds $(FPP_DEFINES) # Debug: -g; Profiling: -p; Openmp: -openmp; Endianness: -convert big_endian
export LFLAGS   = -openmp -lpthread -fpp -nologo -warn all -i-static -convert big_endian -traceback -D_USE_NIX -mcmodel=large # Profiling: -p; Openmp: -openmp -lpthread;endianness: -convert big_endian
export LLFLAGS  =
export MKFLAGS =
export AR = ar rc
export SUFFLIB = .lib
export SUFFPROG =

#Software repository
export SRCREP = ../../../Software

# MohidBase1
export BASE1INC = ../Mohid_Base_1
export BASE1 = Mohid_Base_1$(SUFFLIB)

# MohidBase2
export BASE2INC = ../Mohid_Base_2
export BASE2 = Mohid_Base_2$(SUFFLIB)

#ConvertToNetcdf
export CONVERT2NETCDFINC = ../SmallTools/Convert2netcdf

# HDF5 lib
export LHDF5FORTRAN = libhdf5_fortran.a
export LHDF5 = libhdf5.a
export LHDF5HL = libhdf5_hl.a

# Z lib
#export ZLIB = libz.a
export ZLIB = libz.so

# Netcdf lib
export LNETCDF  = libnetcdf.a
export LNETCDFF = libnetcdff.a

# Curl lib
export LCURL = libcurl.so

# libfproj4 static lib
export LPROJ4F = libfproj4.a

# All libs folders
export BASELIBS := \
       $(BASE2INC)/$(BASE2) \
       $(BASE1INC)/$(BASE1)

export HDFLIBS := \
       $(HDF5LIB)/$(LHDF5FORTRAN) \
       $(HDF5LIB)/$(LHDF5) \
       $(HDF5LIB)/$(LHDF5HL) \
       $(ZLIBINC)/$(ZLIB)

# All libs folders (including netcdf)
# Order of libraries *is relevant* at link-time
export NETCDFLIBS := \
	$(NETCDFLIB)/$(LNETCDFF) \
	$(NETCDFLIB)/$(LNETCDF) \
	$(CURLLIB)/$(LCURL)

export PROJ4FLIBS := \
    $(PROJ4FLIB)/$(LPROJ4F)

export MODULENETCDFOBJ := \
	$(CONVERT2NETCDFINC)/ModuleNETCDF.$(O)

#------Files and modules lists------

METAFILES = \
        README \
        Editme_template.smk \
        Nix.smk

MODULES := 
MODULES := $(MODULES) Makefiles
MODULES := $(MODULES) MohidWater
MODULES := $(MODULES) Mohid_Base_2
MODULES := $(MODULES) Mohid_Base_1
MODULES := $(MODULES) MohidLand
MODULES := $(MODULES) ConvertToHDF5
MODULES := $(MODULES) ConvertToXYZ
# MODULES := $(MODULES) MohidRiver
# MODULES := $(MODULES) HDF5Extrator

NCMODULES :=
NCMODULES := $(NCMODULES) SmallTools/Convert2netcdf

ifeq ($(IS_NETCDF),true)
MODULES := $(MODULES) $(NCMODULES)
endif

#------Makefile rules---------------

include Makefiles/Makemodules.mk

#------Modules dependencies----------

MohidLand.all : Mohid_Base_2.all
MohidWater.all : Mohid_Base_2.all
ifeq ($(IS_NETCDF),true)
SmallTools/Convert2netcdf.all : Mohid_Base_1.all
ConvertToHDF5.all : Mohid_Base_2.all \
                    SmallTools/Convert2netcdf.all
Mohid_Base_2.all : Mohid_Base_1.all \
		SmallTools/Convert2netcdf.all
else
ConvertToHDF5.all : Mohid_Base_2.all
Mohid_Base_2.all : Mohid_Base_1.all
endif
ConvertToXYZ.all : Mohid_Base_2.all
MohidRiver.all : Mohid_Base_2.all
HDF5Extrator.all : Mohid_Base_2.all

