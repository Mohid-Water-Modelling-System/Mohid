
NAME            = MohidWater

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
OBJ_CCFLAGS     = $(CCFLAGS) -I$(INC_HDF5) -I$(INC_NETCDF) -I$(INC_PROJ4) -I$(INC_PROJ4F) -Iinclude -I../Mohid_Base_1/include -I../Mohid_Base_2/include

LIBS            = \
                -L$(LIB_HDF5) -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
                -L$(LIB_NETCDF) -lnetcdf -lnetcdff \
                -L$(LIB_PROJ4F) -lfproj4 \
                -L$(LIB_PROJ4) -lproj \
                -L$(LIB_ZLIB) -lz \
                -lm

MODOUT          = -module 

OBJ             = build/ModuleGOTM.o build/ModuleTurbGOTM.o build/ModuleFreeVerticalMovement.o \
                  build/ModuleToga.o build/ModuleGauge.o build/ModuleOil.o \
                  build/ModuleOil_0D.o build/ModuleHNS.o build/ModuleOpenBoundary.o \
                  build/ModuleTurbulence.o build/ModuleHydrodynamicFile.o build/ModuleAssimilation.o \
                  build/ModuleWaves.o build/ModuleJet.o build/ModuleSand.o \
                  build/ModuleConsolidation.o build/ModuleHydrodynamic.o build/ModuleWaterProperties.o \
                  build/ModuleLagrangian.o build/ModuleLagrangianGlobal.o build/ModuleSedimentProperties.o \
                  build/ModuleSediment.o build/ModuleInterfaceSedimentWater.o build/ModuleInterfaceWaterAir.o \
                  build/ModuleModel.o
#                  build/ModuleSequentialAssimilation.o \

.PHONY: clean all

all: bin/MohidWater.exe

bin/MohidWater.exe: ${OBJ}
	${FC} ${OBJ_CCFLAGS} $? ../Mohid_Base_1/build/*.o ../Mohid_Base_2/build/*.o src/MohidWater.F90 ${LIBS} -o bin/MohidWater.exe

build/ModuleGOTM.o: src/ModuleGOTM.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleTurbGOTM.o: src/ModuleTurbGOTM.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleFreeVerticalMovement.o: src/ModuleFreeVerticalMovement.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleToga.o: src/ModuleToga.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleGauge.o: src/ModuleGauge.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleOil.o: src/ModuleOil.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleOil_0D.o: src/ModuleOil_0D.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleHNS.o: src/ModuleHNS.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleOpenBoundary.o: src/ModuleOpenBoundary.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleTurbulence.o: src/ModuleTurbulence.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleHydrodynamicFile.o: src/ModuleHydrodynamicFile.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleAssimilation.o: src/ModuleAssimilation.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleWaves.o: src/ModuleWaves.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleJet.o: src/ModuleJet.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleSand.o: src/ModuleSand.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleConsolidation.o: src/ModuleConsolidation.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleHydrodynamic.o: src/ModuleHydrodynamic.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleWaterProperties.o: src/ModuleWaterProperties.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleLagrangian.o: src/ModuleLagrangian.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleLagrangianGlobal.o: src/ModuleLagrangianGlobal.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleSedimentProperties.o: src/ModuleSedimentProperties.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleSediment.o: src/ModuleSediment.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleInterfaceSedimentWater.o: src/ModuleInterfaceSedimentWater.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleInterfaceWaterAir.o: src/ModuleInterfaceWaterAir.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleSequentialAssimilation.o: src/ModuleSequentialAssimilation.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleModel.o: src/ModuleModel.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

install:
	ln -sf ../src/MohidWater/bin/MohidWater.exe ../../bin/MohidWater.exe

clean:
	-rm -f build/*.o
	-rm -f include/*.mod
	-rm -f bin/*.exe
