# Makefile CCCP

This section describes a Cross-Compiler, Cross-Platform (CCCP) [makefile](makefile.md) development methodology using the [make](makefile.md) tool from GNU, also known as gmake. The syntax to use would be:

```
make platform.action
```

*(ex: make Win, make Unix.clean, etc...).*

## Description

A CCCP makefile project a bunch of library, modules and program folders within a single root folder.

### Tree description

#### Example tree

```
> tree cproj/
cproj/
|-- Editme.mk
|-- Makefile
|-- Makefiles
|   |-- Makefiles.mk
|   |-- Makemodules.mk
|   `-- Makefile
|-- Nix.mk
|-- lib1
|   `-- Makefile
|-- lib2
|   `-- Makefile
`-- prog1
    `-- Makefile
```

#### Regular folders

Every folder must contain the following files:

Makefile

encapsulates a list of all the folder's source-files and meta-files, a bunch of includes containing the makefile rules and, last but not least, the target dependencies. Several local variables are defined within:

INCS

LIBS

SRCDIR

contains the path to the source files

FILES

contains the list of source files (using the $(S) suffix)

METAFILES

contains the list of files required by the Makefile to run in this project (using the *.smk* suffix) besides the *Makefile* itself (of course).

TARGET

is the target name that the local makefile is supposed to build.

MODULES

contains the list of modules (special root folder only).

Dependencies

is a list of the dependent files and the respective dependencies.

#### Special folders and special files

There are two special folders, mandatory, and each containing special files, beyond the mandatory regular files listed above:

1. **the *root* folder**

    Makefile

    it's completely different from the regular folder *Makefile*. It contains the rules to build the Solution for a series of platform.

    Nix.mk

    it contains the regular folder *Makefile* structure and content. It also contains all the necessary configuration global variables specific to the CCCP solution. It must be configured for a specific platform (Linux, windows,etc...) by the administrator.

    Editme.mk

    Each user must edit this configuration file before building the solution. It contains the paths to required external libraries such as [HDF5](../../tools/hdf5.md) or [Netcdf](../../archive/linux/netcdf.md).

2. **the *Makefiles* folder**

    Makemodules.mk

    contains the makefile rules to build (all), clean, install, etc, the target modules. It is invoked only from the root file *Nix.mk*.

    Makefiles.mk

    contains the makefile rules to build, clean, install, etc, the target files.

### The makefile's global variables

The idea is to use **global variables** to pass arguments to the different modules makefiles. Thus, each specific platform makefile in the root folder defines these global variables. Then, it calls the modules makefiles who will use the global variables in their syntax. If a systematic approach of using global variables is implemented, then it is fairly easy to have a cross-platform makefile project implemented; and easily extended.

## A template CCCP makefile solution

### Special root folder's *Makefile*

```
 SHELL = /bin/sh

 export MAKE = make

 .PHONY: default
 default:
 	@echo Please choose your platform from:
 	@echo Type "make (win nix).(' '  clean install sos)".

 #make all (add other platforms "Nix Nix Nix_double Win" ...)
 PLAT := Nix
 .PHONY : $(PLAT)
 $(PLAT):
 	$(MAKE) -f $@.mk all

 #make clean
 PLATCLEAN := $(addsuffix .clean, $(PLAT))
 .PHONY : $(PLATCLEAN)
 $(PLATCLEAN) :
 	$(MAKE) -f $(@:.clean=).mk clean

 #make install
 PLATINST := $(addsuffix .install, $(PLAT))
 .PHONY : $(PLATINST)
 $(PLATINST) :
 	$(MAKE) -f $(@:.install=).mk install
```

### Special root folder's *Editme.mk*

```
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
```

### Special root folder's *Nix.mk*

```
SHELL = /bin/sh

#------User configuration file---------------

include Editme.mk

#-- NIX platform specific global variables --

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
export SUFFLIB = .lib
export SUFFPROG =
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
export NETCDFLIBS := \
                     $(BASELIBS) \
                     $(NETCDFINC)/$(LNETCDF)

#------Files and modules lists------

METAFILES = \
       README \
       Editme_template.smk \
       Nix.smk
MODULES = \
         Makefiles \
         MohidWater \
         Mohid_Base_2 \
         Mohid_Base_1

#------Makefile rules---------------

include Makefiles/Makemodules.mk

#------Modules dependencies----------

Mohid_Base_2.all : Mohid_Base_1.all
MohidWater.all : Mohid_Base_2.all
```

### Special Makefiles folder's *Makefile*

```
METAFILES = \
        Makefiles.smk \
        Makemodules.smk

SHELL = /bin/sh

.PHONY: all
all:
        @echo No such action in this module
```

### Special Makefiles folder's *Makefiles.mk*

```
OBJS = $(FILES:.$(S)=.$(O))
SRCF = $(SRCDIR)
SUFF := $(suffix $(TARGET))

#Make all
.PHONY: all
all: $(TARGET)

$(TARGET) : $(OBJS) $(LIBS)
ifeq ($(SUFF),$(SUFFLIB))
        @$(AR) $@ $^
else
        @$(CC) $(LFLAGS) -o $@ $^ $(LLFLAGS)
endif
        @echo Finished building $@.

#Fortran compilation rule
%.$(O) : $(SRCF)/%.$(F)
        @$(CC) $(CCFLAGS) $(INCS) $<
        @echo $* .................. [OK]

#make clean
.PHONY: clean
clean:
        @-$(DEL) *.$(O) *.$(MOD) $(TARGET)
        @echo erased $(TARGET) files.

#make install
SOURCE := $(TARGET)
.PHONY: install
install: $(SOURCE)
        @-$(CP) $< $(DESTDIR)/`date +%G%m%d`_$(addprefix $(VER), $<)
        @echo Installed $<.
```

### Special Makefiles folder's *Makemodules.mk*

```
#make all
MODALL = $(addsuffix .all, $(MODULES))
.PHONY: all $(MODALL)
all: $(MODALL)
$(MODALL):
	-$(MAKE) -C $(@:.all=) all

#make clean
MODCLEAN = $(addsuffix .clean, $(MODULES))
.PHONY: clean $(MODCLEAN)
clean: $(MODCLEAN)
$(MODCLEAN):
	-$(MAKE) -C $(@:.clean=) clean

#make install
MODINSTALL = $(addsuffix .install, $(MODULES))
.PHONY: install $(MODINSTALL)
install: $(MODINSTALL)
$(MODINSTALL) :
	-$(MAKE) -C $(@:.install=) install
```

### Regular folder typical *Makefile*

```
INCS = -I$(HDF5)
LIBS =
SRCDIR = $(SRCBASE1)
FILES = \
       ModuleGlobalData.$(S) \
       ModuleDrainageNetwork.$(S) \
       ModuleProfile.$(S)
TARGET = $(BASE1)

SHELL = /bin/sh

include ../Makefiles/Makefiles.mk

ModuleMacroAlgae.${O} : ModuleFunctions.${O}
ModuleTimeSerie.${O} : ModuleEnterData.${O}
ModuleDischarges.${O} : ModuleFunctions.${O} \
                     ModuleTimeSerie.${O} \
                     ModuleDrawing.${O}
```

## See also

- [Makefile](makefile.md)
- [Mohid Makefile](mohid-makefile.md)
