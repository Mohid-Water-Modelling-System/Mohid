SHELL = /bin/sh

#------User configuration file---------------

include MohidLand_Config.mk

#-- NIX platform specific global variables --

export ORIGIN = ORIGIN
export CP  = sudo cp
export DEL = rm
export F = F90
export Obj = o
export MOD = mod
export CC= ifort
export CCFLAGS  = -openmp -fpp -r8 -std03 -m64 -ip -zero -xHost -O3 -parallel -inline-level=1 -assume buffered_io -warn all -nologo -convert big_endian -fpe0 -D_USE_NIX -traceback -mcmodel=small -static-intel -openmp-link=static -static-libgcc -heap-arrays 64 -c $(FPP_DEFINES) 
export LFLAGS   = -openmp -fpp -r8 -std03 -m64 -ip -zero -xHost -O3 -parallel -inline-level=1 -assume buffered_io -warn all -nologo -convert big_endian -fpe0 -D_USE_NIX -traceback -mcmodel=small -static-intel -openmp-link=static -static-libgcc -heap-arrays 64 -Wl,-rpath,'$$$$ORIGIN/.',-rpath,$(INTEL_PATH) $(FPP_DEFINES)
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

# HDF5 lib
export LHDF5FORTRANLIB = libhdf5_fortran.a
export LHDF5LIB = libhdf5.a
export LHDF5HLLIB = libhdf5_hl.a

# Z lib
export ZLIB = libz.so
export SZLIB = libsz.so

# All libs folders
export BASELIBS := \
       $(BASE2INC)/$(BASE2) \
       $(BASE1INC)/$(BASE1)

export HDFLIBS := \
       $(HDF5LIBINC)/$(LHDF5FORTRANLIB) \
       $(HDF5LIBINC)/$(LHDF5LIB) \
       $(HDF5LIBINC)/$(LHDF5HLLIB) \
       $(ZLIBINC)/$(ZLIB) \
       $(ZLIBINC)/$(SZLIB)

# All libs folders (including netcdf)
# Order of libraries *is relevant* at link-time

#------Files and modules lists------

METAFILES = \
        README \
        Editme_template.smk \
        Nix.smk

MODULES := 
MODULES := $(MODULES) Makefiles
MODULES := $(MODULES) Mohid_Base_2
MODULES := $(MODULES) Mohid_Base_1
MODULES := $(MODULES) MohidLand

#------Makefile rules---------------

include Makefiles/Makemodules.mk

#------Modules dependencies----------

MohidLand.all : Mohid_Base_2.all
Mohid_Base_2.all : Mohid_Base_1.all
