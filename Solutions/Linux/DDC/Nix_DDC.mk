SHELL = /bin/sh

#------User configuration file---------------

include Editme.mk

#-- NIX platform specific global variables --

export CP  = sudo cp
export DEL = rm
export Obj = o
export F = F90
export MOD = mod
#export CC= ifort
export CC= mpif90
#export CCFLAGS  = -c      -r8 -inline-level=0 -fpp -warn all -nologo -convert big_endian -fpe0 -D_USE_NIX -traceback -mcmodel=large -heap-arrays 64 -openmp -check bounds $(FPP_DEFINES) # Debug: -g; Profiling: -p; Openmp: -openmp; Endianness: -convert big_endian
export CCFLAGS  = -c  -O3  -r8 -inline-level=2 -fpp -warn all -nologo -convert big_endian -fpe0 -D_USE_NIX -traceback  -mcmodel=large -heap-arrays 64 -check bounds $(FPP_DEFINES) -I/opt/mpich/include # Debug: -g; Profiling: -p; Openmp: -openmp; Endianness: -convert big_endian
export LFLAGS   = -O3 -lpthread -fpp -nologo -warn all -convert big_endian  -D_USE_NIX -mcmodel=large  -shared-intel  # Profiling: -p; Openmp: -openmp -lpthread;endianness: -convert big_endian; -shared-intel -i-static
export LLFLAGS  =
export MKFLAGS =
export AR = ar rc
export SUFFLIB = .lib
export SUFFPROG =

#Software repository
export SRCREP = ../../../Software

# MohidBase1
export BASE1INC = ../../Mohid_Base_1
export BASE1 = Mohid_Base_1$(SUFFLIB)

# HDF5 lib
export HDF5MODINC = ../../../../../hdf5/linux/include
export LHDF5FORTRAN = libhdf5_fortran.a
export LHDF5 = libhdf5.a
export LHDF5HL = libhdf5_hl.a

# Z lib
#export ZLIB = libz.a
export ZLIB = libz.so
export SZLIB = libsz.so


# All libs folders
export BASELIBS := \
       $(BASE1INC)/$(BASE1)

export HDFLIBS := \
       $(HDF5LIB)/$(LHDF5FORTRAN) \
       $(HDF5LIB)/$(LHDF5) \
       $(HDF5LIB)/$(LHDF5HL) \
       $(ZLIBINC)/$(ZLIB) \
       $(ZLIBINC)/$(SZLIB)


#------Files and modules lists------

METAFILES = \
        README \
        Editme_template.smk \
        Nix.smk

MODULES :=
MODULES := $(MODULES) Makefiles
MODULES := $(MODULES) DDCParser
MODULES := $(MODULES) DDCWorker
MODULES := $(MODULES) Mohid_Base_1

NCMODULES :=

#------Makefile rules---------------

include Makefiles/Makemodules.mk

#------Modules dependencies----------

DDCParser.all : Mohid_Base_1.all
DDCWorker.all : Mohid_Base_1.all

