# Makefile

A makefile is a script run by make (or nmake under windows) that contains information about the dependencies of a series of targets and, additionally, makes them. For example, make can build the binaries for a large software project.

## Syntax

Here's the syntax to use the make command:

```
make Target
make -C Subdir Target
make -f filename Target
```

*(ex: make -f mymakefile all).*

In most projects in linux, the most common targets are:

- *all*: builds the full project.
- *install*: copies the executable files in their final destination on the system.
- *clean*: deletes the object and the executable files created by target *all*.

*NOTE: Here's how to ignore built-in pattern rules!!*

```
> make -r target
```

This is very important for the following rules:

```
%.o : %.c
%.dvi : %.aux
```

## [A large project cross-platform, cross-compiler makefile methodology](makefile-cccp.md)

[Here](makefile-cccp.md) is described a CCCP makefile methodology particularly suitable for [MOHID](../../overview/mohid.md).

## Sample makefiles

In particular, the gnumake program is very well documented [here](http://www.gnu.org/software/make/manual/make.html).

### [Mohid\_Base\_1](../../water/README.md#mohid-base-1#mohid-base-1)

The source code of a sample **Makefile** is given bellow: `[makefile,Y]`

```
#
#  Use GNUmake!
#  Win32: rename gmake.exe to make.exe and load ifortvars.bat
#  *nix: load ifortvars.sh
#

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .f90 .o

DEL = rm -f

BIN = ./bin
MOD = ./mod
SRC = .
OBJ = .

#-----Users: edit these lines!!!----------
#*nix vs windows
#Uncomment one of these
#HDF5     = ../IntelLibs
#MODOPT   = -module:$(MOD)
#OBJOPT   = -object:$*.o
HDF5     = /opt/hdf5/hdf5/lib
MODOPT   = -module $(MOD)
OBJOPT   = -o $*.o
#-----Users: end of block ----------------

TARGET = $(BIN)/mohidbase1.lib

CC       = ifort
CCFLAGS  = -c -fpp -nologo
AR       = ar rc
#--------Users: fill these lines with your sourcesafe data-------------
SERVER   =
USER     =
PASS     =
SOSHOME  =
PROJECT  = Mohid_Base_1
GET      = soscmd -command GetFile \
           -server $(SERVER):8890 \
           -name $(USER) -password $(PASS) \
           -database "W:\SourceSafe\Mohid_v4\srcsafe.ini" \
           -project $(PROJECT) \
           -soshome $(SOSHOME) \
           -file
#--------End of sourcesafe data-------------

#Rules to apply to a standard .o file
.f90.o:
         @$(CC) $(CCFLAGS) $< $(OBJOPT) $(MODOPT)
         @echo  $* ......... [OK]

OBJS = \
       ModuleGlobalData.o \
       ModuleTime.o \
       ModuleEnterData.o \
       ModuleFunctions.o \
       ModuleBenthos.o \
       ModuleCEQUALW2.o \
       ModuleTimeSerie.o \
       ModuleDischarges.o \
       ModuleDrawing.o \
       ModuleLUD.o \
       ModuleWaterQuality.o \
       ModuleSedimentQuality.o \
       ModuleLife.o \
       ModuleInterface.o \
       ModuleHydroIntegration.o \
       ModuleLightExtinction.o \
       ModuleStopWatch.o \
       ModuleTriangulation.o \
       ModuleHDF5.o \
       ModuleDrainageNetwork.o \
       ModuleProfile.o

.PHONY: all clean

all: $(TARGET)
	@echo
	@echo Finished building $(TARGET)
	@echo

$(TARGET) : $(OBJS)
	@$(AR) $(TARGET) $(OBJS)
	@echo  $(TARGET) ....... [OK]

clean :
	@-$(DEL) *.o $(MOD)/*.mod $(BIN)/*.lib
 	@echo $(TARGET) ......... [Erased]

#This command here gets out of sourcesafe the latest version
sos: $(OBJS:.o=)

$(OBJS:.o=)::
        @-$(GET) $@.f90
        @echo  $@.f90 ....... [Updated]

#--------Dependencies list-------------

ModuleGlobalData.o : ModuleGlobalData.f90

ModuleTime.o : ModuleGlobalData.o \
              ModuleTime.f90

ModuleEnterData.o : ModuleTime.o \
                   ModuleEnterData.f90

ModuleFunctions.o : ModuleEnterData.o \
                   ModuleFunctions.f90

ModuleBenthos.o : ModuleEnterData.o \
                 ModuleBenthos.f90

ModuleCEQUALW2.o : ModuleFunctions.o \
                  ModuleCEQUALW2.f90

ModuleTimeSerie.o : ModuleEnterData.o \
                   ModuleTimeSerie.f90

ModuleDischarges.o : ModuleFunctions.o \
                    ModuleTimeSerie.o \
                    ModuleDischarges.f90

ModuleDrawing.o : ModuleFunctions.o \
                 ModuleDrawing.f90

ModuleLUD.o : ModuleGlobalData.o \
             ModuleLUD.f90

ModuleWaterQuality.o : ModuleFunctions.o \
                      ModuleLUD.o \
                      ModuleWaterQuality.f90

ModuleSedimentQuality.o : ModuleFunctions.o \
                         ModuleLUD.o \
                         ModuleSedimentQuality.f90

ModuleLife.o : ModuleFunctions.o \
              ModuleLife.f90

ModuleInterface.o : ModuleWaterQuality.o \
                   ModuleSedimentQuality.o \
                   ModuleCEQUALW2.o \
                   ModuleLife.o \
                   ModuleBenthos.o \
                   ModuleInterface.f90

ModuleHydroIntegration.o : ModuleTime.o \
                          ModuleHydroIntegration.f90

ModuleLightExtinction.o : ModuleFunctions.o \
                         ModuleTimeSerie.o \
                         ModuleLightExtinction.f90

ModuleStopWatch.o : ModuleTime.o \
                   ModuleStopWatch.f90

ModuleTriangulation.o : ModuleGlobalData.o \
                       ModuleTriangulation.f90

#The following module requires the HDF5 mod files
ModuleHDF5.o : ModuleGlobalData.o \
              ModuleHDF5.f90
	@$(CC) $(CCFLAGS) -I$(HDF5) $*.f90 $(OBJOPT) $(MODOPT)
	@echo  $* ......... [OK]

ModuleDrainageNetwork.o : ModuleDischarges.o \
                         ModuleLightExtinction.o \
                         ModuleStopWatch.o \
                         ModuleInterface.o \
                         ModuleDrainageNetwork.f90

ModuleProfile.o : ModuleEnterData.o \
                 ModuleHDF5.o \
                 ModuleProfile.f90

#--------End of dependencies list---------
```

### ADCP to AVU

```
SHELL = /bin/sh

#CCCPEV (Cross-Compiler, Cross-Platform Environment Variables)
CC=gcc
CCFLAGS=-c
LCFLAGS=
OCFLAG=-o
RM=del
H=h
C=c
O=o
A=exe
D=tds
L=lib

#Files list
TARGET=adcp_avu.$(A)
SRC=adcp_avu_this.$(C) readstring.$(C)
OBJ=$(SRC:.$(C)=.$(O))
LIB=
HDR=$(SRC:.$(C)=.$(H))

#Action "make all"
.PHONY=all
all : $(TARGET)
	-$(RM) $(OBJ) $(^:.$(A)=.$(D))

$(TARGET) : $(LIB) $(OBJ)
	$(CC) $(LCFLAGS) $(OCFLAG)$@ $^

$(LIB) :
	@echo Warning: $@ is missing.

%.$(O) : %.$(C)
	$(CC) $(CCFLAGS) $<

#Action "make clean"
.PHONY=clean
clean :
	-$(RM) $(OBJ) $(TARGET) $(TARGET:.$(A)=.$(D))

#Dependencies list
%.$(C) : %.$(H)  #Fortran is %.$(H) : %.$(O)
adcp_avu_this.$(O) : readstring.$(H)
```

### Mururoa

Here's another sample [makefile](http://en.wikipedia.org/wiki/Make) made for \\\\guillaume\\projects\\mururoa based on [\*NIX\_platforms](../../archive/linux/nix-platforms.md) and on [this](http://www.eng.hawaii.edu/Tutor/Make/index.html) tutorial.

```
SHELL = /bin/sh

#Suffixes which we'll be working on
.SUFFIXES: .cpp .o

#windows/linux bash commands
COPY = xcopy /-Y
DEL = del /Q

BIN = ../bin
INC = ../inc
SRC = .
OBJ = .
LIB = ../lib

CC       = gcc
CCFLAGS  = -c -I$(INC)
LFLAGS   = -L$(LIB)
CPP      = g++

#Rules to apply to a standard .o file
.cpp.o:
	$(CC) $(CCFLAGS) ./$<

OBJS = \
       prog_3d.o      \
       prog_prm.o     \
       prog_io.o      \
       prog_cdf_out.o \
       prog_cdf.o     \
       prog_blc.o     \
       prog_all.o     \
       msnemo.o

INCS = \
       $(INC)/StdAfx.h     \
       $(INC)/prog_3d.h    \
       $(INC)/prog_all.h   \
       $(INC)/prog_blc.h   \
       $(INC)/prog_boo.h   \
       $(INC)/prog_cdf.h   \
       $(INC)/prog_cdf_out.h   \
       $(INC)/prog_ind.h   \
       $(INC)/prog_io.h    \
       $(INC)/prog_prm.h   \
       $(INC)/prog_mai.h   \
       $(INC)/prog_switches.h  \

LIBS = \
       $(LIB)/netcdf.lib

#Compiles the full solution and cleans the intermediate files
all: msnemo postclean

#Makes the main executable
#links the object files against the libraries
msnemo: $(OBJS)
	$(CPP) $(LFLAGS) -o $(BIN)/$@.exe $(OBJS) $(LIBS)

#Cleans to the original distro
clean:
	-$(DEL) *.o

#Cleans everything except binaries, source and object files
postclean:

#Installs the binaries in the right places of the system
install:

#Dependencies
prog_3d.o: $(INCS)
prog_prm.o: $(INCS)
prog_io.o: $(INCS)
prog_cdf_out.o: $(INCS)
prog_cdf.o: $(INCS)
prog_blc.o: $(INCS)
prog_all.o: $(INCS)
msnemo.o: $(INCS)
```

## Troubleshoots

### use DFWIN

This is a windows platform specific library. It cannot be ported to Linux. To build the code in linux you need to change it first. This is currently(20061201) the case in *ModuleGlueHDF5*. If you get an error attempting to make it in windows, first make sure you're calling **make** inside an Intel Fortran Development Console, which is where all the environment variables are correct.

### **ifconsol**

If you get an error about an "ifconsol", then it's probably because you have incorrect environment variables. Under windows, make sure you're inside an Intel Fortran Development Console.

## Other references

1. The Wikipedia [article](http://en.wikipedia.org/wiki/Make).
2. The [GNUMake](http://www.gnu.org/software/make/manual/make.html) documentation. A classic standard.
3. [Here's](http://www.eng.hawaii.edu/Tutor/Make/index.html) a popular tutorial.
4. The 20060629 workshop [powerpoint](http://www.mohid.com/hydrogroup/presentations/HG_060629_Formacao_Makefiles.ppt).
