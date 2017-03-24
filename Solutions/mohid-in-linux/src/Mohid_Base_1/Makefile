
NAME            = MohidBase_1

FC = ifort
AR = ar rcs

LIB_DIR         = $(HOME)/apps_intel

DEBUG_FLAGS     = -g -traceback
LANG_FLAGS      = -cpp -real_size 64
WARNINGS_FLAGS  = -w
OPT_FLAGS       = -convert little_endian -fPIC
OTH_FLAGS       = -xHost -ip
MODOUT          = -module

FPP_DEFINES     = -D_LAGRANGIAN_GLOBAL_ -D_USE_NIX -D_STACK_LIMITS_ -D_USE_PROJ4

ZLIB            = ${LIB_DIR}/zlib-1.2.11
INC_ZLIB        = $(ZLIB)/include
LIB_ZLIB        = $(ZLIB)/lib

HDF5            = ${LIB_DIR}/hdf5-1.8.15
INC_HDF5        = $(HDF5)/include
LIB_HDF5        = $(HDF5)/lib

NETCDF          = ${LIB_DIR}/netcdf-4.4.1.1
INC_NETCDF      = $(NETCDF)/include
LIB_NETCDF      = $(NETCDF)/lib

PROJ4           = ${LIB_DIR}/proj-4.9.3
INC_PROJ4       = $(PROJ4)/include
LIB_PROJ4       = $(PROJ4)/lib

PROJ4F          = ${LIB_DIR}/proj4-fortran
INC_PROJ4F      = $(PROJ4F)/include
LIB_PROJ4F      = $(PROJ4F)/lib

CCFLAGS         = ${WARNINGS_FLAGS} ${DEBUG_FLAGS} ${LANG_FLAGS} ${OPT_FLAGS} ${OTH_FLAGS} ${FPP_DEFINES}
OBJ_CCFLAGS     = $(CCFLAGS) -I$(INC_HDF5) -I$(INC_NETCDF) -I$(INC_PROJ4) -I$(INC_PROJ4F) -Iinclude

LIBS            = \
                -L$(LIB_HDF5) -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
                -L$(LIB_NETCDF) -lnetcdf -lnetcdff \
                -L$(LIB_PROJ4F) -lfproj4 \
                -L$(LIB_PROJ4) -lproj \
                -L$(LIB_ZLIB) -lz \
                -lm

MODOUT          = -module 

OBJ             = build/ModuleGlobalData.o build/ModuleLUD.o build/ModuleTriangulation.o build/ModuleTime.o build/ModuleEnterData.o build/ModuleWWTPQ.o \
                  build/ModuleStopWatch.o build/ModuleFunctions.o build/ModuleMacroAlgae.o build/ModuleWaterQuality.o build/ModuleSedimentQuality.o \
                  build/ModuleHydroIntegration.o build/ModuleSeagrassWaterInteraction.o build/ModuleHDF5.o build/ModuleHDF5_OO.o \
                  build/ModuleSeagrassSedimInteraction.o build/ModuleLife.o build/ModuleCEQUALW2.o build/ModuleBenthos.o build/ModuleDrawing.o build/ModuleProfile.o \
                  build/ModuleBivalve.o build/ModuleBenthicEcology.o build/ModuleInterface.o build/ModuleTimeSerie.o build/ModuleDischarges.o \
                  build/ModuleLightExtinction.o build/ModuleDrainageNetwork.o build/ModuleTask2000.o

.PHONY: all clean

all: lib/libmohidbase1.a

lib/libmohidbase1.a: ${OBJ}
	${AR} $@  ${OBJ}

build/ModuleGlobalData.o: src/ModuleGlobalData.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleLUD.o: src/ModuleLUD.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleTriangulation.o: src/ModuleTriangulation.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleTime.o: src/ModuleTime.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleEnterData.o: src/ModuleEnterData.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleWWTPQ.o: src/ModuleWWTPQ.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleStopWatch.o: src/ModuleStopWatch.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleFunctions.o: src/ModuleFunctions.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleMacroAlgae.o: src/ModuleMacroAlgae.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleWaterQuality.o: src/ModuleWaterQuality.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleSedimentQuality.o: src/ModuleSedimentQuality.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleHydroIntegration.o: src/ModuleHydroIntegration.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleSeagrassWaterInteraction.o: src/ModuleSeagrassWaterInteraction.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleHDF5.o: src/ModuleHDF5.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleHDF5_OO.o: src/ModuleHDF5_OO.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleSeagrassSedimInteraction.o: src/ModuleSeagrassSedimInteraction.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleLife.o: src/ModuleLife.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleCEQUALW2.o: src/ModuleCEQUALW2.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleBenthos.o: src/ModuleBenthos.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleDrawing.o: src/ModuleDrawing.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleProfile.o: src/ModuleProfile.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleBivalve.o: src/ModuleBivalve.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleBenthicEcology.o: src/ModuleBenthicEcology.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include 

build/ModuleInterface.o: src/ModuleInterface.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleTimeSerie.o: src/ModuleTimeSerie.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleDischarges.o: src/ModuleDischarges.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleLightExtinction.o: src/ModuleLightExtinction.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleDrainageNetwork.o: src/ModuleDrainageNetwork.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleTask2000.o: src/ModuleTask2000.F
	${FC} -c $? ${OBJ_CCFLAGS} -o $@ ${MODOUT}include

clean:
	-rm -f build/*.o
	-rm -f include/*.mod
	-rm -f lib/*.so lib/*.a 
