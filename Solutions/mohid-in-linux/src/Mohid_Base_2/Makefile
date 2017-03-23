
NAME            = MohidBase_2

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
OBJ_CCFLAGS     = $(CCFLAGS) -I$(INC_HDF5) -I$(INC_NETCDF) -I$(INC_PROJ4) -I$(INC_PROJ4F) -Iinclude -I../Mohid_Base_1/include

LIBS            = \
                -L$(LIB_HDF5) -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
                -L$(LIB_NETCDF) -lnetcdf -lnetcdff \
                -L$(LIB_PROJ4F) -lfproj4 \
                -L$(LIB_PROJ4) -lproj \
                -L$(LIB_ZLIB) -lz \
                -lm

MODOUT          = -module 

OBJ             = build/ModuleHorizontalGrid.o build/ModuleStatistic.o build/ModuleGridData.o build/ModuleBasinGeometry.o \
                  build/ModuleHorizontalMap.o build/ModuleBoxDif.o build/ModuleGeometry.o build/ModuleMap.o \
                  build/ModuleAdvectionDiffusion.o build/ModuleInterpolation.o build/ModuleNETCDF.o \
                  build/ModuleField4D.o build/ModuleFillMatrix.o build/ModuleChainReactions.o build/ModuleAtmosphere.o
                  
.PHONY: all clean

all: lib/libmohidbase2.a

lib/libmohidbase2.a: ${OBJ}
	${AR} $@  ${OBJ}

build/ModuleHorizontalGrid.o: src/ModuleHorizontalGrid.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleStatistic.o: src/ModuleStatistic.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleGridData.o: src/ModuleGridData.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleBasinGeometry.o: src/ModuleBasinGeometry.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleHorizontalMap.o: src/ModuleHorizontalMap.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleBoxDif.o: src/ModuleBoxDif.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@  ${MODOUT}include

build/ModuleGeometry.o: src/ModuleGeometry.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@  ${MODOUT}include

build/ModuleMap.o: src/ModuleMap.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@  ${MODOUT}include

build/ModuleAdvectionDiffusion.o: src/ModuleAdvectionDiffusion.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleInterpolation.o: src/ModuleInterpolation.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleNETCDF.o: src/ModuleNETCDF.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleField4D.o: src/ModuleField4D.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleFillMatrix.o: src/ModuleFillMatrix.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@  ${MODOUT}include

build/ModuleChainReactions.o: src/ModuleChainReactions.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@  ${MODOUT}include

build/ModuleAtmosphere.o: src/ModuleAtmosphere.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@  ${MODOUT}include


clean:
	-rm -f build/*.o
	-rm -f include/*.mod
	-rm -f lib/*.so lib/*.a 
