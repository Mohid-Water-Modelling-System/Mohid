
NAME            = MohidLand

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

OBJ             = build/ModuleRunOff.o build/ModuleRunoffProperties.o build/ModulePorousMedia.o \
                  build/ModulePorousMediaProperties.o build/ModuleVegetation.o build/ModuleSnow.o \
                  build/ModuleIrrigation.o build/ModuleReservoirs.o build/ModuleBasin.o

.PHONY: clean all

all: bin/MohidLand.exe

bin/MohidLand.exe: ${OBJ}
	${FC} ${OBJ_CCFLAGS} $? ../Mohid_Base_1/build/*.o ../Mohid_Base_2/build/*.o src/MohidLand.F90 ${LIBS} -o bin/MohidLand.exe

build/ModuleRunOff.o: src/ModuleRunOff.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleRunoffProperties.o: src/ModuleRunoffProperties.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModulePorousMedia.o: src/ModulePorousMedia.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModulePorousMediaProperties.o: src/ModulePorousMediaProperties.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleVegetation.o: src/ModuleVegetation.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleSnow.o: src/ModuleSnow.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleIrrigation.o: src/ModuleIrrigation.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleReservoirs.o: src/ModuleReservoirs.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

build/ModuleBasin.o: src/ModuleBasin.F90
	${FC} -c $? ${OBJ_CCFLAGS} ${LIBS} -o $@ ${MODOUT}include

install:
	ln -sf ../src/MohidLand/bin/MohidLand.exe ../../bin/MohidLand.exe

clean:
	-rm -f build/*.o
	-rm -f include/*.mod
	-rm -f bin/*.exe
