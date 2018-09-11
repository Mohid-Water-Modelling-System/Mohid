#!/usr/bin/bash
#==============================================================================
#title           : compile_mohid.sh
#description     : This script is an attempt to compile MOHID in a linux
#                  machine with Intel compiler ifort.
#author          : Jorge Palma (jorgempalma@tecnico.ulisboa.pt)
#date            : 20180712
#usage           : bash compile_mohid.sh
#notes           :
#==============================================================================

### Make changes to fit your setup ###
set -e

#### Choose your compiler ####
FC=ifort

#### Debug mode ####
IS_DEBUG=false

#### libraries path ####
DIR_REQ=$HOME/apps_intel

ZLIB=$DIR_REQ/zlib-1.2.11
HDF5=$DIR_REQ/hdf5-1.8.15
NETCDF=$DIR_REQ/netcdf-4.4.1.1
MPI=$DIR_REQ/mpich-3.2
PROJ4=$DIR_REQ/proj-4.9.3
PROJ4F=$DIR_REQ/proj4-fortran
IPHREEQC=$DIR_REQ/iphreeqc-3.3.11-12535
PHREEQCRM=$DIR_REQ/phreeqcrm-3.3.11-12535

#### Activate modules ####
USE_OPENMP=true                  # default : true
USE_MPI=false                    # default : false
USE_HDF=true                     # default : true
USE_NETCDF=true                  # default : true
USE_PROJ4=true                   # default : true
USE_IEEE_ARITHMETIC=true         # default : true
USE_CUDA=false                   # default : false
USE_PHREEQC=false                # default : false
USE_PAGELOCKED=false             # default : false
USE_PAD_MATRICES=false           # default : false
USE_BIG_MAX_COLUMNS=false        # default : false
USE_LAGRANGIAN=false             # default : false
USE_LAGRANGIAN_GLOBAL=true       # default : true
USE_GUI=false                    # default : false
USE_GGI=false                    # default : false
USE_BFM=false                    # default : false
USE_INCREASE_MAXINSTANCES=false  # default : false
USE_ONLINE=false                 # default : false
USE_NIX=true                     # default : true
USE_OUTPUT_OFF=false             # default : false
USE_STACK_LIMITS=true            # default : true
USE_SCANLINE_UNSAFE=false        # default : false
USE_GOOGLEMAPS=false             # default : false
USE_SEDIMENT=false               # default : false  # true => error
USE_WAVES=false                  # default : false  # true => error
USE_AIR=false                    # default : false
USE_MODIS=false                  # default : false
USE_SEQASSIMILATION=false        # default : false
USE_OPENMI=false                 # default : false

########################################################################
########################################################################
START_OF_COMPILE=`date`

#set -x #echo on
SRCREP=`pwd`
SRCREP_SRC=$SRCREP/src
SRCREP_BIN=$SRCREP/bin
SRCREP_TEST=$SRCREP/test
MACHINE_TYPE=`uname -m`

FPP_DEFINES=""
if [ $USE_MPI == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_USE_MPI"
fi
if [ $USE_OPENMI == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_OPENMI_"
fi
if [ $USE_NETCDF == false ]; then
  FPP_DEFINES="$FPP_DEFINES -D_NO_NETCDF"
fi
if [ $USE_HDF == false ]; then
  FPP_DEFINES="$FPP_DEFINES -D_NO_HDF5"
fi
if [ $USE_IEEE_ARITHMETIC == false ]; then
  FPP_DEFINES="$FPP_DEFINES -D_NOT_IEEE_ARITHMETIC"
fi
if [ $USE_CUDA == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_USE_CUDA -D_ENABLE_CUDA"
fi
if [ $USE_PROJ4 == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_USE_PROJ4"
fi
if [ $USE_PHREEQC == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_PHREEQC_"
  if [ $MACHINE_TYPE == 'x86_64' ]; then
    FPP_DEFINES="$FPP_DEFINES -D_PHREEQC_X64_"
  fi
fi
if [ $USE_PAGELOCKED == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_USE_PAGELOCKED"
fi
if [ $USE_PAD_MATRICES == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_PAD_MATRICES"
fi
if [ $USE_BIG_MAX_COLUMNS == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_BIG_MAX_COLUMNS"
fi
if [ $USE_LAGRANGIAN == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_LAGRANGIAN_"
fi
if [ $USE_LAGRANGIAN_GLOBAL == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_LAGRANGIAN_GLOBAL_"
fi
if [ $USE_GUI == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_GUI_"
fi
if [ $USE_GGI == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_CGI_"
fi
if [ $USE_BFM == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_BFM_"
fi
if [ $USE_INCREASE_MAXINSTANCES == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_INCREASE_MAXINSTANCES"
fi
if [ $USE_ONLINE == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_ONLINE_"
fi
if [ $USE_NIX == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_USE_NIX"
fi
if [ $USE_OUTPUT_OFF == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_OUTPUT_OFF_"
fi
if [ $USE_STACK_LIMITS == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_STACK_LIMITS_"
fi
if [ $USE_SCANLINE_UNSAFE == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_SCANLINE_UNSAFE"
fi
if [ $USE_GOOGLEMAPS == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_GOOGLEMAPS"
fi
if [ $USE_SEDIMENT == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_SEDIMENT_"
fi
if [ $USE_WAVES == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_WAVES_"
fi
if [ $USE_AIR == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_AIR_"
fi
if [ $USE_MODIS == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_USE_MODIS"
fi
if [ $USE_SEQASSIMILATION == true ]; then
  FPP_DEFINES="$FPP_DEFINES -D_USE_SEQASSIMILATION"
fi
########################
#-- specific global variables --
DEL="rm"
F90=".F90"
F77=".F"
Obj=".o"
SUFFLIB=".a"
Mod=".mod"
Exe=".exe"
AR="ar rcs"

LANG_FLAGS=
DEBUG_FLAGS=
WARNINGS_FLAGS=
OPT_OPENMP=
OPT_FLAGS=
OTH_FLAGS=
if [[ $FC == *"gfortran"* ]]; then
    LANG_FLAGS='-cpp -fdefault-real-8 -fdefault-double-8'
    if [ $IS_DEBUG == true ]; then
        DEBUG_FLAGS="-g3 -fcheck=all " #-fbacktrace -fdump-fortran-optimized -fbounds-check"
        # -fsanitize=undefined -fsanitize=address -fsanitize=leak"
        # -ffpe-trap=invalid,zero,overflow,underflow"

        WARNINGS_FLAGS="-Wall -Wno-unused-function -Wno-conversion -Wno-unused-variable -Wno-tabs -Wno-c-binding-type"
        #-Wmaybe-uninitialized -Wunused-dummy-argument  -Wcharacter-truncation
    else
        WARNINGS_FLAGS="-w -pedantic"
    fi
    if [ $USE_OPENMP == true ]; then
        OPT_OPENMP='-fopenmp'
    fi
    OPT_FLAGS="-O2 -ffast-math -march=x86-64 -fconvert=little-endian -fPIC -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans"
    MODOUT="-J"
elif [[ $FC == *"ifort"* ]]; then
    LANG_FLAGS="-cpp -real_size 64 -heap-arrays"
    if [ $IS_DEBUG == true ]; then
        DEBUG_FLAGS="-g -traceback -heap-arrays "
    fi
    WARNINGS_FLAGS="-w"
    if [ $USE_OPENMP == true ]; then
        OPT_OPENMP='-qopenmp'
    fi
    OPT_FLAGS="-O2 -convert little_endian -fPIC"
    OTH_FLAGS="-g -traceback -xHost -ip -fpe0"
    MODOUT="-module "
fi

CCFLAGS="$WARNINGS_FLAGS $DEBUG_FLAGS $LANG_FLAGS $OPT_OPENMP $OPT_FLAGS $OTH_FLAGS $FPP_DEFINES"

############################################
## settings Includes / Libraries ##
BASE1_SRC=$SRCREP_SRC/Mohid_Base_1
BASE2_SRC=$SRCREP_SRC/Mohid_Base_2

MPIINC=$MPI/include
MPILIB=$MPI/lib
MPI_LIB="-lmpi -lmpifort"

ZLIBINC=$ZLIB/include
ZLIBLIB=$ZLIB/lib
Z_LIB=lz

HDF5INC=$HDF5/include
HDF5LIB=$HDF5/lib
HDF5_LIB="-lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5"

NETCDFINC=$NETCDF/include
NETCDFLIB=$NETCDF/lib
NETCDF_LIB="-lnetcdf -lnetcdff"

PROJ4INC=$PROJ4/include
PROJ4LIB=$PROJ4/lib
PROJ4_LIB=-lproj

PROJ4FINC=$PROJ4F/include
PROJ4FLIB=$PROJ4F/lib
PROJ4F_LIB=-lfproj4

IPHREEQCINC=$IPHREEQC/include
IPHREEQCLIB=$IPHREEQC/lib
IPHREEQC_LIB=-liphreeqc

PHREEQCRMINC=$PHREEQCRM/include
PHREEQCRMLIB=$PHREEQCRM/lib
PHREEQCRM_LIB=-lphreeqcrm

############################################
############################################

#### Includes / Libraries ####
INCLUDES=
LIBS=
if [ $USE_NETCDF == true ]; then
  INCLUDES="$INCLUDES -I${NETCDFINC}"
  LIBS="$LIBS -L${NETCDFLIB} $NETCDF_LIB"
fi
if [ $USE_HDF == true ]; then
  INCLUDES="$INCLUDES -I${HDF5INC}"
  LIBS="$LIBS -L${HDF5LIB} $HDF5_LIB -ldl"
fi
if [ $USE_PROJ4 == true ]; then
  INCLUDES="$INCLUDES -I${PROJ4INC} -I${PROJ4FINC}"
  LIBS="$LIBS -L${PROJ4LIB} $PROJ4_LIB -L${PROJ4FLIB} $PROJ4F_LIB"
fi
if [ $USE_MPI == true ]; then
  INCLUDES="$INCLUDES -I${MPIINC}"
  LIBS="$LIBS -L${MPILIB} $MPI_LIB"
fi
if [ $USE_PHREEQC == true ]; then
  INCLUDES="$INCLUDES -I${IPHREEQCINC} -I${PHREEQCRMINC}"
  LIBS="$LIBS -L${IPHREEQCLIB} $IPHREEQC_LIB -L${PHREEQCRMLIB} $PHREEQCRM_LIB"
fi

INCLUDES="$INCLUDES -I${BASE1_SRC}/include -I${BASE2_SRC}/include"
LIBS="$LIBS -I${ZLIBINC} -L${ZLIBLIB} -lz  -lm"


#### -Don't touch anything below this line- ####

# echo colors ##
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color
OK=${GREEN}OK${NC}
NOK=${RED}NOK${NC}
ERROR=${RED}ERROR${NC}
WARNING=${RED}warning${NC}

mkdir -p bin


#### FUNCTIONS ####
## remove all pre-compile files ##
CLEAN(){
  find . -name \*.a -type f -delete
  find . -name \*.so -type f -delete
  find . -name \*.exe -type f -delete
  find . -name \*.o -type f -delete
  find . -name \*.mod -type f -delete
  find . -name \*.log -type f -delete
  find . -name UsedKeyWords_1.dat -type f -delete
  find src/  -name \*.exe -type l -delete
  find bin/  -name \*.exe -type l -delete
}

## program help ##
HELP(){
           echo
           echo "Compile Mohid script"
           echo "Usage: $0 [-h] [-c] [-v]   [-option1] [-option2] [-option3] ..."
           echo "    -h|-help                  : Show this help"
           echo "    -v|-verbose               : Show more text"
           echo "    -c|clean                  : Remove build files"
           echo
           echo " Options to compile"
           echo "    -mb1|mohid_base_1         : Compile Mohid Base 1"
           echo "    -mb2|mohid_base_2         : Compile Mohid Base 2"
           echo "    -ml|mohid_land            : Compile Mohid Land"
           echo "    -mw|mohid_water           : Compile Mohid Water"
           echo "    -mr|mohid_river           : Compile Mohid River (experimental)"
           echo "    -mt|mohid_tools           : Compile Mohid Tools"
           echo
}

## compile Mohid_Base ##
COMPILE_MOHID_BASE(){
  array=$1[@]
  base=$2
  modules_Mohid_Base=("${!array}")

  rm -f build/*${Obj} include/*${Mod} lib/*.a lib/*.so

  if [ $OPT_VERBOSE == 1 ]; then
    echo
    echo " $FC -c $CCFLAGS src/....$F90 $INCLUDES $LIBS  -o build/....${Obj}  ${MODOUT}include"
    echo
  fi

  for module in ${modules_Mohid_Base[*]}; do
    if [[ $USE_HDF == false && ($module == *"HDF5"* || $module == ModuleProfile || $module == ModuleDrainageNetwork) ]]; then
      continue
    fi
    if [[ $USE_HDF == false && ($module == ModuleHorizontalGrid || ModuleStatistic) ]]; then
      continue
    fi
    if [[ $USE_NETCDF == false && $module == *"NETCDF"* ]]; then
      continue
    fi
    if [[ $USE_PHREEQC == false && ($module == IPhreeqc_interface || $module == IPhreeqcRM_interface || $module == ModulePhreeqCRM || $module == ModulePhreeqC) ]]; then
      continue
    fi

    echo -ne " compiling $module ...                                                      \r"
    if [ -f "src/${module}$F90" ]; then
      $FC -c $CCFLAGS src/${module}$F90 $INCLUDES $LIBS -o build/${module}${Obj} ${MODOUT}include
    elif [ -f "src/${module}$F77" ]; then
      $FC -c $CCFLAGS src/${module}$F77 $INCLUDES $LIBS -o build/${module}${Obj} ${MODOUT}include
    else
      echo -e "${ERROR} src/${module}$F90 File not found!"
      exit 0
    fi

    if [ $module == IPhreeqc_interface ]; then
      if [ ! -f "include/iphreeqc${Mod}" ]; then
        echo -e "${ERROR} include/iphreeqc${Mod} File not created!"
        exit 0
      fi
    elif [ $module == IPhreeqcRM_interface ]; then
      if [ ! -f "include/phreeqcrm${Mod}" ]; then
        echo -e "${ERROR} include/phreeqcrm${Mod} File not created!"
        exit 0
      fi
    elif [ ! -f "include/${module,,}$Mod" ]; then
      echo -e "${ERROR} include/${module,,}$Mod File not created!"
      exit 0
    fi
  done

  #gfortran -shared -o libmohid_base_1.so *${Obj}
  $AR lib/lib${base}.a  build/*${Obj}

  echo -e " compile $base ${OK}                                                      "
  echo
}

## compile Mohid ##
COMPILE_MOHID(){
  array=$1[@]
  name=$2
  modules_Mohid=("${!array}")

  rm -f build/*${Obj} include/*${Mod}

  if [ $OPT_VERBOSE == 1 ]; then
    echo
    echo " $FC -c $CCFLAGS src/....$F90 $INCLUDES $LIBS  -o build/....${Obj}  ${MODOUT}include"
    echo
  fi

  for module in ${modules_Mohid[*]}; do
    echo -ne " compiling $module ...                                                      \r"

   if [ -f "src/${module}$F90" ] || [ -f "src/${module}$F77" ]; then
      $FC -c $CCFLAGS src/${module}$F90 $INCLUDES $LIBS -o build/${module}${Obj} ${MODOUT}include
    else
      echo -e "${ERROR} src/${module}$F90 File not found!"
      exit 0
    fi

    if [ ! -f "include/${module,,}$Mod" ];  then
      echo -e "${ERROR} include/${module,,}$Mod File not created!"
      exit 0
    fi
  done

  echo -ne " compiling ${name} ...                                                  \r"

  if [ $OPT_VERBOSE == 1 ]; then
    echo
    echo " $FC $CCFLAGS build/*${Obj}  $BASE1_SRC/build/*${Obj}  $BASE2_SRC/build/*${Obj}  src/${name}${F90}  $INCLUDES -Iinclude $LIBS  -o bin/${name}${Exe}"
    echo
  fi

  if [ $name == 'MohidWater' ]; then
    $FC $CCFLAGS build/*${Obj}  $BASE1_SRC/build/*${Obj}  $BASE2_SRC/build/*${Obj}  src/Main${F90}  $INCLUDES -Iinclude $LIBS  -o bin/${name}${Exe}
  else
    $FC $CCFLAGS build/*${Obj}  $BASE1_SRC/build/*${Obj}  $BASE2_SRC/build/*${Obj}  src/${name}${F90}  $INCLUDES -Iinclude $LIBS  -o bin/${name}${Exe}
  fi

  if [ ! -f "bin/${name}${Exe}" ]; then
    echo -e "${ERROR} ${name}${Exe} File not created!"
    exit 0
  else
    echo -e " compile ${name} ${OK}                                                      "
    echo
  fi

  cd $SRCREP_BIN; ln -sf ../src/${name}/bin/${name}${Exe} .
  
  if [ $name == 'MohidWater' ]; then
    ln -sf $SRCREP/bin/MohidWater.exe $SRCREP_TEST/mohidwater/25m_deep/exe/MohidWater.exe
  fi
  
  if [ $name == 'MohidLand' ]; then
     ln -sf $SRCREP/bin/MohidLand.exe $SRCREP_TEST/mohidland/schematicWatershed/exe/MohidLand.exe
  fi
  
}

## compile Mohid_Tools ##
COMPILE_MOHID_TOOLS(){
  array=$1[@]
  tool=$2
  module_array=("${!array}")

  if [ $OPT_VERBOSE == 1 ]; then
    echo
    echo " $FC -c $CCFLAGS src/....$F90 $INCLUDES $LIBS  -o build/....${Obj}  ${MODOUT}include"
    echo
  fi

  echo -ne "\n"
  for module in ${module_array[*]}; do
    if [[ $USE_HDF == false || $USE_NETCDF == false ]]; then
      return
    fi
    echo -ne "  compiling $module ...                                                      \r"
    if [ -f "src/${module}$F90" ]; then
      $FC -c $CCFLAGS src/${module}$F90 $INCLUDES $LIBS -o build/${module}${Obj} ${MODOUT}include
      if [ ! -f "include/${module,,}$Mod" ]; then
        echo -e "  ${ERROR} ${module,,}$Mod File not created!"
        exit 0
      fi
    else
      echo -e "  ${ERROR} ${module}$F90 File not found!"
      exit 0
    fi
  done

  echo -ne " compiling $tool ...                                                  \r"

  if [ $OPT_VERBOSE == 1 ]; then
    echo
    echo " $FC $CCFLAGS build/*${Obj}  $BASE1_SRC/build/*${Obj}  $BASE2_SRC/build/*${Obj}  src/${tool}${F90} $INCLUDES -Iinclude $LIBS  -o bin/${tool}${Exe}"
    echo
  fi

  $FC $CCFLAGS build/*${Obj}  $BASE1_SRC/build/*${Obj}  $BASE2_SRC/build/*${Obj}  src/${tool}${F90} $INCLUDES -Iinclude $LIBS  -o bin/${tool}${Exe}
  if [ ! -f "bin/${tool}${Exe}" ]; then
    echo -e "${ERROR} ${tool}${Exe} File not created!"
    exit 0
  else
      echo -e " compile $tool ${OK}                                                      "
      echo
  fi
}

## Mohid_Base_1 ##
MOHID_BASE_1(){
  echo
  echo "#### Mohid Base 1 ####"

  cd $BASE1_SRC

  modules_Mohid_Base_1=( \
    ModuleGlobalData \
    ModuleLUD \
    ModuleTriangulation \
    ModuleTime \
    ModuleEnterData \
    ModuleWWTPQ \
    ModuleStopWatch  \
    ModuleFunctions  \
    ModuleMacroAlgae  \
    ModuleWaterQuality  \
    ModuleSedimentQuality  \
    ModuleHydroIntegration  \
    ModuleSeagrassWaterInteraction  \
    ModuleHDF5  \
    ModuleHDF5_OO  \
    ModuleSeagrassSedimInteraction  \
    ModuleLife \
    ModuleCEQUALW2  \
    ModuleBenthos \
    ModuleDrawing  \
    ModuleProfile \
    ModuleBivalve  \
    ModuleBenthicEcology  \
    IPhreeqc_interface
    IPhreeqcRM_interface
    ModulePhreeqCRM
    ModulePhreeqC
    ModuleInterface  \
    ModuleTimeSerie  \
    ModuleDischarges  \
    ModuleLightExtinction  \
    ModuleDrainageNetwork  \
    ModuleMPImanagement  \
    ModuleTask2000)

  COMPILE_MOHID_BASE modules_Mohid_Base_1 "mohidbase1"
}

## Mohid_Base_2 ##
MOHID_BASE_2(){
  echo
  echo "#### Mohid Base 2 ####"

  cd $BASE2_SRC

  modules_Mohid_Base_2=( \
    ModuleHorizontalGrid \
    ModuleStatistic \
    ModuleGridData \
    ModuleBasinGeometry \
    ModuleHorizontalMap \
    ModuleBoxDif \
    ModuleGeometry  \
    ModuleMap  \
    ModuleAdvectionDiffusion  \
    ModuleInterpolation  \
    ModuleNETCDF  \
    ModuleField4D  \
    ModuleFillMatrix  \
    ModuleChainReactions  \
    ModuleAtmosphere      \
    ModuleTwoWay)

  COMPILE_MOHID_BASE modules_Mohid_Base_2 "mohidbase2"
}

## Mohid_Land ##
MOHID_LAND(){
  echo
  echo "#### Mohid Land ####"

  cd $SRCREP_SRC/MohidLand

  modules_Mohid_Land=( \
    ModuleRunOff \
    ModuleRunoffProperties \
    ModulePorousMedia \
    ModulePorousMediaProperties \
    ModuleVegetation \
    ModuleSnow \
    ModuleIrrigation  \
    ModuleReservoirs  \
    ModuleBasin)

  COMPILE_MOHID modules_Mohid_Land "MohidLand"
}

## Mohid_Water ##
MOHID_WATER(){
  echo
  echo "#### Mohid Water ####"

  cd $SRCREP_SRC/MohidWater

  modules_Mohid_Water=( \
    ModuleTurbine  \
    ModuleGOTM  \
    ModuleTurbGOTM  \
    ModuleFreeVerticalMovement  \
    ModuleToga  \
    ModuleGauge  \
    ModuleOil  \
    ModuleOil_0D  \
    ModuleHNS  \
    ModuleOpenBoundary  \
    ModuleTurbulence  \
    ModuleHydrodynamicFile  \
    ModuleAssimilation  \
    ModuleWaves  \
    ModuleJet  \
    ModuleSand  \
    ModuleConsolidation  \
    ModuleHydrodynamic  \
    ModuleWaterProperties  \
    ModuleLagrangian  \
    ModuleLagrangianGlobal  \
    ModuleSedimentProperties  \
    ModuleSediment  \
    ModuleInterfaceSedimentWater  \
    ModuleInterfaceWaterAir  \
    
    #ModuleSequentialAssimilation  \
    ModuleModel)

  COMPILE_MOHID modules_Mohid_Water "MohidWater"
}

## Mohid_River ##
MOHID_RIVER(){
  echo
  echo "#### Mohid River ####"

  cd $SRCREP_SRC/MohidRiver/src
  rm -f *${Obj} *${Mod}

  echo -ne " compiling Mohid River ...                                                  \r"
  $FC $CCFLAGS $BASE1_SRC/*${Obj}  $BASE2_SRC/*${Obj}  RiverNetwork${F90}  $INCLUDES $LIBS  -o ../bin/MohidRiver${Exe}

  if [ ! -f "../bin/MohidRiver${Exe}" ]; then
    echo -e "${ERROR} MohidRiver${Exe} File not created!"
    exit 0
  else
    echo -e " compile Mohid River ${OK}                                                      "
  fi

  cd $SRCREP_BIN; ln -sf ../src/MohidRiver/bin/MohidRiver${Exe} .
}

## Mohid_Tools ##
MOHID_TOOLS(){
  echo
  echo "#### Mohid Tools ####"

  if [[ $USE_HDF == false || $USE_NETCDF == false ]]; then
      echo -e " ${WARNING}: USE_HDF and USE_NETCDF must be true and both libraries must be installed"
      echo
      exit 1
  fi
  modules_Mohid_Tools=( \
    BasinDelimiter \
    ConvertGridDataToHDF5 \
    ConvertHDF5ToGridData \
    ConvertToHDF5 \
    ConvertToXYZ
    DigitalTerrainCreator \
    FillMatrix \
    HDF5Exporter \
    HDF5Extractor \
    HDF5Statistics \
    DomainDecompositionConsolidation \
    #Shell \   ## error
    )

  for tool in ${modules_Mohid_Tools[*]}; do
    cd $SRCREP_SRC/Tools/src/$tool
    rm -f build/*${Obj} include/*${Mod}
    rm -f bin/*.exe
  done

  for tool in ${modules_Mohid_Tools[*]}; do
    cd $SRCREP_SRC/Tools/src/$tool
    echo -ne " compiling $tool ...                                                      \r"

    if [ $tool = 'ConvertToHDF5' ]; then
      modules_ConvertToHDF5=( \
        ncdflib \
        ModuleAladinFormat \
        ModuleARPSFormat \
        ModuleARPSToWW3 \
        ModuleCFFormat \
        ModuleCFPolcomFormat \
        #ModuleConvertModisL2 \        ## obsolete (need hdf4) => USE_MODIS=false
        #ModuleConvertModisL3Mapped \  ## obsolete (need hdf4) => USE_MODIS=false
        #ModuleConvertOceanColorL2 \   ## obsolete (need hdf4) => USE_MODIS=false
        ModuleCOWAMAAsciiWind \
        ModuleERA40Format \
        ModuleHYCOMFormat \
        ModuleIHRadarFormat
        ModuleMERCATORFormat \
        ModuleMOG2DFormat \
        ModuleReadSWANNonStationary \
        ModuleWOAFormat \
        ModuleNetCDFCF_2_HDF5MOHID \
        ModuleEUCenterFormat \
        ModuleGFSasciiWind \
        ModuleGlueHDF5Files \
        ModuleHDF5ToASCIIandBIN \
        ModuleHellermanRosensteinAscii \
        ModuleInterpolateGrids \
        ModuleInterpolateTime \
        ModuleLevitusFormat \
        ModuleMM5Format \
        ModulePatchHDF5Files \
        ModuleSWAN \
        ModuleTecnoceanAscii \
        ModuleWRFFormat \
        ModuleCALMETFormat)
        COMPILE_MOHID_TOOLS modules_ConvertToHDF5 "$tool"

    elif [ $tool = 'ConvertToXYZ' ]; then
      modules_ConvertToXYZ=( \
        ModuleASCII \
        ModuleEtopo2 \
        ModuleEtopo5 \
        ModuleGEBCO \
        ModuleNASA \
        ModuleNOAA_ShoreLine \
        ModuleSRTM30)
        COMPILE_MOHID_TOOLS modules_ConvertToXYZ "$tool"

    elif [ $tool = 'HDF5Exporter' ]; then
      modules_HDF5Exporter=( \
        ModuleExportHDF5ToTimeSerie)
        COMPILE_MOHID_TOOLS modules_HDF5Exporter "$tool"

    elif [ $tool = 'HDF5Extractor' ]; then
      modules_HDF5Extractor=( \
        ModuleHDF5Extractor)
        COMPILE_MOHID_TOOLS modules_HDF5Extractor "$tool"

    elif [ $tool = 'HDF5Statistics' ]; then
      modules_HDF5Statistics=( \
        ModuleHDF5Statistics)
        COMPILE_MOHID_TOOLS modules_HDF5Statistics "$tool"

    elif [ $tool = 'DomainDecompositionConsolidation' ]; then
      modules_DDC=( \
        ModuleHashTable \
        ModuleDDC)
        COMPILE_MOHID_TOOLS modules_DDC "$tool"

    elif [ $tool = 'Shell' ]; then
      modules_Shell=( \
        ModuleShell)
        COMPILE_MOHID_TOOLS modules_Shell "$tool"

    elif [ -f "src/${tool}$F90" ];  then

      if [ $OPT_VERBOSE == 1 ]; then
        echo
        echo " $FC $CCFLAGS $BASE1_SRC/build/*${Obj} $BASE2_SRC/build/*${Obj} src/${tool}$F90  $INCLUDES $LIBS -o bin/${tool}${Exe}"
        echo
      fi

      $FC $CCFLAGS $BASE1_SRC/build/*${Obj} $BASE2_SRC/build/*${Obj} src/${tool}$F90  $INCLUDES $LIBS -o bin/${tool}${Exe}

      if [ ! -f "bin/${tool}${Exe}" ];  then
        echo -e "${ERROR} bin/${tool}${Exe} File not created!"
        exit 0
      else
        echo -e " compile ${tool} ${OK}                                                      "
        echo
      fi

    elif [ -f "src/${tool}.f90" ];  then

      if [ $OPT_VERBOSE == 1 ]; then
        echo
        echo " $FC $CCFLAGS $BASE1_SRC/build/*${Obj} $BASE2_SRC/build/*${Obj} src/${tool}.f90  $INCLUDES $LIBS -o bin/${tool}${Exe}"
        echo
      fi

      $FC $CCFLAGS $BASE1_SRC/build/*${Obj} $BASE2_SRC/build/*${Obj} src/${tool}.f90  $INCLUDES $LIBS -o bin/${tool}${Exe}

      if [ ! -f "bin/${tool}${Exe}" ];  then
        echo -e "${ERROR} bin/${tool}${Exe} File not created!"
        exit 0
      else
        echo -e " compile ${tool} ${OK}                                                      "
        echo
      fi

    else
      echo -e "${ERROR} ${tool}$F90 File not found!"
      exit 0
    fi

    cd $SRCREP_BIN; ln -sf ../src/Tools/src/${tool}/bin/${tool}${Exe} .

  done
}

# ----------- Main -------------- ##

# If no argument given, show HELP and exit
OPT_VERBOSE=0
if [ $# -lt 1 ] ;then
  HELP
  exit 0
fi

while [ $# -gt 0 ]; do
  case $1 in
    -h|-help|--help)
       HELP
       exit 0
    ;;
    -v|-verbose)
       OPT_VERBOSE=1
       shift
    ;;
    -c|-clean|--clean)
       CLEAN
       exit 0
    ;;
    -mb1|mohid_base_1)
       MOHID_BASE_1
       shift
    ;;
    -mb2|mohid_base_2)
       MOHID_BASE_2
       shift
    ;;
    -ml|mohid_land)
       MOHID_LAND
       shift
    ;;
    -mw|mohid_water)
       MOHID_WATER
       shift
    ;;
    -mr|mohid_river)
       MOHID_RIVER
       shift
    ;;
    -mt|mohid_tools)
       MOHID_TOOLS
       shift
    ;;
    *)
       echo -e " ${ERROR}: unrecognized command line option $1. Use \"$0 -h\" for help."
       echo
       shift
  esac
done

END_OF_COMPILE=`date`
echo "=========================================================================="
echo "build started:    $START_OF_COMPILE"
echo "build completed:  $END_OF_COMPILE"
echo
echo "--->                  Executables ready                               <---"
echo
ls -l $SRCREP_BIN --color=auto
echo
echo "=========================================================================="

exit 0
