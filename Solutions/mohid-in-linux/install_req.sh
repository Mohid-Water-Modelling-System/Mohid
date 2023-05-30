
#!/bin/bash
#==============================================================================
#title        : install_req.sh
#description  : This script is an attempt to compile all the necessary libraries
#               to compile MOHID in a machine with Ubuntu or CentOS linux distro
#               and Intel compiler. For more information consult
#               http://wiki.mohid.com and http://forum.mohid.com
#author       : Jorge Palma (jorgempalma@tecnico.ulisboa.pt)
#date         : 20180712
#usage        : bash install_req.sh
#notes        :
#==============================================================================

### Make the changes to fit your setup ###

source /opt/intel/oneapi/setvars.sh

## default path to compiler
export CC=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/icc
export CXX=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/icpc
export FC=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/ifort
export F77=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/ifort


## intel flags compiler
export CFLAGS='-O3 -xHost -ip -DIFORT'
export CXXFLAGS='-O3 -xHost -ip -DIFORT'
export FCFLAGS='-O3 -xHost -ip -DIFORT'
export FFLAGS='-O3 -xHost -ip -DIFORT'

## install libraries path
DIRINSTALL=$HOME/apps_intel

## have you root permissions? yes or no (uncomment what is best for you)
#sudo=sudo     ## yes
sudo=          ## no

## libraries to intall
zlib='zlib-1.2.13'
hdf5='hdf5-1.8.17'
#netcdf='netcdf-4.4.1.1'
netcdf='netcdf-c-4.8.1'
#netcdff='netcdf-fortran-4.4.4'
netcdff='netcdf-fortran-4.5.4'
proj='proj-4.9.3'
proj4fortran='proj4-fortran'
iphreeqc='iphreeqc-3.3.11-12535'
phreeqcrm='phreeqcrm-3.3.11-12535'
mpi=mpich
mpi_version='4.0a2'
#openmpi='openmpi-2.0.1'


####################################################
###### -Don't touch anything below this line- ######
####################################################
set -e

# echo colors ##
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color
OK=${GREEN}OK${NC}
NOK=${RED}NOK${NC}
ERROR=${RED}ERROR${NC}
WARNING=${RED}warning${NC}

#### FUNCTIONS ####
PAUSE(){
   read -rp "$*"
}

## program help ##
HELP(){
           echo
           echo "auto install mohid library"
           echo "Usage: $0 [-option]"
           echo "    -h|-help                : Show this help"
           echo "    -req|requirements       : install compile requirements"
           echo "    -zlib                   : install zlib library"
           echo "    -hdf5                   : install hdf5 library"
           echo "    -nc|netcdf              : install netcdf and netcdff library"
           echo "    -ncf|netcdff            : install netcdf fortran library"
	   echo "    -proj4                  : install proj4 library (optional)"
	   echo "    -proj4f|proj4fortran    : install proj4-fortran wrapper library (optional)"
           echo "    -phqc|iphreeqc          : install iphreeqc library               (optional)"
           echo "    -phqcrm|phreeqcrm       : install phreeqcrm library              (optional)"
           echo "    -mpich                  : install mpich                          (optional)"
           #echo "    -openmpi                : install open-mpi"
           echo "    -rm                     : remove tmp install dir (~/.mohid)"
           echo
}

WARNING(){
    echo
    echo -e " ${GREEN}You must close and reopen your session in another shell to load new env variables added in .bashrc${NC}"
    echo
}

## install requirements ##
OPT_REQ(){
    echo
    echo " #### Install requirements ####"
    #PAUSE 'Press [Enter] key to continue...'
    if ! exist="$(type -p "wget")" || [ -z "$exist" ]; then
      $sudo $PACKMANAG install wget
    fi
    ## to install netcdf library
    if ! exist="$(type -p "m4")" || [ -z "$exist" ]; then
      $sudo $PACKMANAG install m4
    fi
    ## to install proj4-fortran library
    if ! exist="$(type -p "git")" || [ -z "$exist" ]; then
      $sudo $PACKMANAG install git
    fi
    if ! exist="$(type -p "autoconf")" || [ -z "$exist" ]; then
      $sudo $PACKMANAG install autoconf
    fi
    if ! exist="$(type -p "automake")" || [ -z "$exist" ]; then
      $sudo $PACKMANAG install automake
    fi

    if ! exist="$(type -p "gcc")" || [ -z "$exist" ]; then
      $sudo $PACKMANAG install gcc gcc-c++ gcc-gfortran
    fi
    echo -e "${GREEN} All basic requirements are installed${NC}"
    echo
}

## install mpich ##
OPT_MPICH(){
       echo
    echo " #### Install $mpich ####"
    #PAUSE 'Press [Enter] key to continue...'
    mpiv="$mpi-$mpi_version"
    RM_DIR $DIRINSTALL/$mpiv
    if [ ! -e $TMPDIRINSTALL/$mpiv.tar.gz ]; then
      wget "http://www.mpich.org/static/downloads/$mpi_version/$mpiv.tar.gz"
    fi
    tar -xf "$mpiv.tar.gz"
    cd $mpiv || exit 1
    CC=$CC CXX=$CXX F77=$F77 FC=$FC  \
    ./configure --prefix=$DIRINSTALL/$mpiv \
                --enable-fast=O3 --disable-error-checking --without-timing --without-mpit-pvars --with-device=ch4:ofi || exit 1
                ##--with-pm=smpd --with-pmi=smpd
    make || exit 1
    $sudo make install || exit 1
    echo "############ $mpi-$mpi_version ############" >> ~/.bashrc
    echo "MPICH=$DIRINSTALL/$mpiv" >> ~/.bashrc
    echo 'export MPIRUN=$MPICH/bin/mpirun' >> ~/.bashrc
    echo 'export PATH=$PATH:$MPICH/bin:$MPICH/lib:$MPICH/include' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MPICH/lib' >> ~/.bashrc
    echo >> ~/.bashrc
}

## install openmpi ##
OPT_OPENMPI(){
    echo
    echo " #### Install $openmpi ####"
    #PAUSE 'Press [Enter] key to continue...'
    RM_DIR $DIRINSTALL/$openmpi
    wget "https://www.open-mpi.org/software/ompi/v2.0/downloads/$openmpi.tar.gz"
    tar -xf "$openmpi.tar.gz"
    cd $openmpi || exit
    ./configure --prefix=$DIRINSTALL/$openmpi  CC=$CC CXX=$CXX F77=$F77 FC=$FC
    make
    $sudo make install
    echo "############ $openmpi ############" >> ~/.bashrc
    echo "OPENMPI=$DIRINSTALL/$openmpi" >> ~/.bashrc
    echo 'export MPIRUN=$OPENMPI/bin/mpirun' >> ~/.bashrc
    echo 'export PATH=$PATH:$OPENMPI/bin:$OPENMPI/lib:$OPENMPI/include' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OPENMPI/lib' >> ~/.bashrc
    echo >> ~/.bashrc
}

## install zlib library ##
OPT_ZLIB(){
    echo
    echo " #### Install $zlib ####"
    #PAUSE 'Press [Enter] key to continue...'
    RM_DIR $DIRINSTALL/$zlib
    if [ ! -e $TMPDIRINSTALL/$zlib.tar.gz ]; then
      wget "http://zlib.net/$zlib.tar.gz"
    fi
    tar -xf "$zlib.tar.gz"
    cd $zlib
    ./configure --prefix=$DIRINSTALL/$zlib || exit 1
    make || exit 1
    $sudo make install || exit 1
    echo "###### $zlib #########" >> ~/.bashrc
    echo "ZLIB=$DIRINSTALL/$zlib" >> ~/.bashrc
    echo 'export PATH=$PATH:$ZLIB/lib:$ZLIB/include' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ZLIB/lib' >> ~/.bashrc
    echo >> ~/.bashrc
}

## install hdf5 library ##
OPT_HDF5(){
    echo
    echo " #### Install $hdf5 ####"
    echo -e " see ${GREEN}https://software.intel.com/en-us/articles/performance-tools-for-software-developers-building-hdf5-with-intel-compilers${NC}"
    PAUSE "Press [Enter] key to continue..."
    RM_DIR $DIRINSTALL/$hdf5
    if [ ! -e $TMPDIRINSTALL/$hdf5.tar.gz ]; then
      wget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/$hdf5/src/$hdf5.tar.gz"
    fi
    tar -xf "$hdf5.tar.gz"
    cd $hdf5 || exit 1
    CC=$CC FC=$FC  ./configure --with-zlib=$DIRINSTALL/$zlib --prefix=$DIRINSTALL/$hdf5 --enable-fortran --enable-fortran2003 || exit 1
    make || exit 1
    #make check || exit 1
    $sudo make install || exit 1
    echo "####### $hdf5 #######" >> ~/.bashrc
    echo "HDF5=$DIRINSTALL/$hdf5" >> ~/.bashrc
    echo 'export PATH=$PATH:$HDF5/bin:$HDF5/lib:$HDF5/include' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5/lib' >> ~/.bashrc
    echo >> ~/.bashrc
}

## install netcdf C library ##
OPT_NC(){
    RM_DIR $DIRINSTALL/$netcdf
    export CPPFLAGS="-I$DIRINSTALL/$hdf5/include -I$DIRINSTALL/$zlib/include"
    export LDFLAGS="-L$DIRINSTALL/$hdf5/lib -L$DIRINSTALL/$zlib/lib"
    export CFLAGS="$CFLAGS -fPIC"
    if [ ! -e $TMPDIRINSTALL/$netcdf.tar.gz ]; then
     # wget "ftp://ftp.unidata.ucar.edu/pub/netcdf/$netcdf.tar.gz"
     wget "https://downloads.unidata.ucar.edu/netcdf-c/4.8.1/$netcdf.tar.gz"
    fi
    tar -xf "$netcdf.tar.gz"
    cd $netcdf || exit 1
    CC=$CC FC=$FC ./configure --prefix=$DIRINSTALL/$netcdf || exit 1
    make || exit 1
    $sudo make install || exit 1
    echo
    echo "########## $netcdf #######" >> ~/.bashrc
    echo "export NETCDF=$DIRINSTALL/${netcdf}" >> ~/.bashrc
    echo 'export PATH=$PATH:$NETCDF/bin:$NETCDF/lib:$NETCDF/include' >> ~/.bashrc
    echo >> ~/.bashrc
    echo 'export NETCDF_ROOT=$NETCDF' >> ~/.bashrc
    echo 'export NETCDF4_ROOT=$NETCDF' >> ~/.bashrc
    echo 'export NETCDF_LIB=$NETCDF/lib' >> ~/.bashrc
    echo 'export NETCDF_INC=$NETCDF/include' >> ~/.bashrc
    echo >> ~/.bashrc
    echo 'export NETCDF_GF_ROOT=$NETCDF' >> ~/.bashrc
    echo 'export NETCDF4_GF_ROOT=$NETCDF' >> ~/.bashrc
    echo 'export NETCDF_GF_LIB=$NETCDF/lib' >> ~/.bashrc
    echo 'export NETCDF_GF_INC=$NETCDF/include' >> ~/.bashrc
    echo >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_LIB' >> ~/.bashrc
    echo >> ~/.bashrc
    echo 'export CPPFLAGS="$CPPFLAGS -I$NETCDF_INC"' >> ~/.bashrc
    echo 'export LDFLAGS="$LDFLAGS -L$NETCDF_LIB"' >> ~/.bashrc
    echo >> ~/.bashrc
}

## install netcdf fortran library ##
OPT_NCF(){
    export CPPFLAGS="-I$DIRINSTALL/$hdf5/include -I$DIRINSTALL/$zlib/include"
    export LDFLAGS="-L$DIRINSTALL/$hdf5/lib -L$DIRINSTALL/$zlib/lib"
    export NETCDF=$DIRINSTALL/${netcdf}
    export PATH=$PATH:$NETCDF/bin:$NETCDF/lib:$NETCDF/include
    export NETCDF_ROOT=$NETCDF
    export NETCDF4_ROOT=$NETCDF
    export NETCDF_LIB=$NETCDF/lib
    export NETCDF_INC=$NETCDF/include
    export NETCDF_GF_ROOT=$NETCDF
    export NETCDF4_GF_ROOT=$NETCDF
    export NETCDF_GF_LIB=$NETCDF/lib
    export NETCDF_GF_INC=$NETCDF/include
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_LIB
    export CPPFLAGS="$CPPFLAGS -I$NETCDF_INC"
    export LDFLAGS="$LDFLAGS -L$NETCDF_LIB"
    if [ ! -e $TMPDIRINSTALL/$netcdff.tar.gz ]; then
     # wget "ftp://ftp.unidata.ucar.edu/pub/netcdf/$netcdff.tar.gz"
     wget "https://downloads.unidata.ucar.edu/netcdf-fortran/4.5.4/$netcdff.tar.gz"
    fi
    tar -xf "$netcdff.tar.gz"
    cd $netcdff || exit
    CC=$CC FC=$FC ./configure --prefix=$DIRINSTALL/$netcdf || exit 1
    make  || exit 1
    $sudo make install  || exit 1
}

## install proj4 library ##
OPT_PROJ4(){
    cd $TMPDIRINSTALL || exit 1
    echo
    echo " #### Install $proj ####"
    #PAUSE 'Press [Enter] key to continue...'
    RM_DIR $DIRINSTALL/$proj
    if [ ! -e $TMPDIRINSTALL/$proj.tar.gz ]; then
      wget "http://download.osgeo.org/proj/$proj.tar.gz"
    fi
    tar -xf "$proj.tar.gz"
    cd $proj || exit 1
    ./configure --prefix=$DIRINSTALL/$proj CC=$CC FC=$FC || exit 1
    make || exit 1
    $sudo make install || exit 1
     ## must create libproj4.so for proj4-fortran install process
    $sudo ln -sf $DIRINSTALL/$proj/lib/libproj.so $DIRINSTALL/$proj/lib/libproj4.so
    echo "############ $proj ############" >> ~/.bashrc
    echo "PROJ4=$DIRINSTALL/$proj" >> ~/.bashrc
    echo 'export PATH=$PATH:$PROJ4/bin:$PROJ4/lib:$PROJ4/include' >> ~/.bashrc
    echo 'export PROJ_PREFIX=$PROJ4/lib' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PROJ4/lib' >> ~/.bashrc
    echo >> ~/.bashrc
}

## install proj4-fortran library ##
OPT_PROJ4F(){
    cd $TMPDIRINSTALL || exit 1
    echo
    echo " #### Install $proj4fortran ####"
    #PAUSE 'Press [Enter] key to continue...'
    if [ ! -d $TMPDIRINSTALL/$proj4fortran ]; then
      git clone https://github.com/mhagdorn/proj4-fortran.git
    fi
    cd $proj4fortran || exit 1
    git checkout 2865227446959983dbda81d52f999921d8b84ad5
    ./bootstrap
    mv configure.in configure.ac
    ./configure --with-proj4=$DIRINSTALL/$proj --prefix=$DIRINSTALL/$proj4fortran  CC=$CC FC=$FC || exit 1
    ## IF error, add -DIFORT or -D<COMPILER_TAG> in the compile expression
    echo
    echo -e " ${RED}If error occur, you need to specify which F90 compiler you use e.g. for INTEL compiler add -DIFORT; For GNU compiler add -Df2cFortran. See cfortran.h ${NC}"
    echo -e " ${RED}cd $TMPDIRINSTALL/$proj4fortran; Copy compile command, add correct compiler flag and run it again with \"make; make install\"${NC}"
    echo -e " ${RED}After that, you must close and reopen your session in another shell to load new env variables added in .bashrc ${NC}"
    echo
    PAUSE 'Press [Enter] key to continue...'
    echo "############ $proj4fortran ############" >> ~/.bashrc
    echo "PROJ4FORTRAN=$DIRINSTALL/$proj4fortran" >> ~/.bashrc
    echo 'export PATH=$PATH:$PROJ4FORTRAN/lib:$PROJ4FORTRAN/include' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PROJ4FORTRAN/lib' >> ~/.bashrc
    echo >> ~/.bashrc
    make || exit 1
    $sudo make install || exit 1
}

## install iphreeqc library ##
## https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
OPT_IPHC(){
    cd $TMPDIRINSTALL || exit
    echo
    echo " #### Install $iphreeqc ####"
    #PAUSE 'Press [Enter] key to continue...'
    RM_DIR $DIRINSTALL/$iphreeqc
    wget "ftp://brrftp.cr.usgs.gov/pub/charlton/iphreeqc/$iphreeqc.tar.gz"
    tar -xf "$iphreeqc.tar.gz"
    cd $iphreeqc || exit
    ## The module IPhreeqc has been revised when using Fortran.
    ## IPhreeqc is now a Fortran “module”. The old IPhreeqc.f.inc and IPhreeqc.f90.inc files are no longer used to define the interfaces for the subroutines of IPhreeqc.
    ## The include files (.inc) are now replaced with “USE” statements (USE IPhreeqc). In addition, an interface file (IPhreeqc_interface.F90) must be included in the user’s Fortran project.
    ## IPhreeqc now uses ISO_C_BINDING, which is available in the Fortran 2003 standard. Use of this standard removes some ambiguities in argument types when calling C.
    ##
    ## In Fortran, you will need to include the source file IPhreeqc_interface.F90 in your project files. This file defines the IPhreeqc Fortran module.
    ## This is the preferred method to use IPhreeqc from a Fortran program.
    ./configure --prefix=$DIRINSTALL/$iphreeqc --enable-fortran-module=yes --enable-fortran-test=yes  CC=$CC FC=$FC
    make
    make check
    $sudo make install
    echo "############ $iphreeqc ############" >> ~/.bashrc
    echo "IPHREEQC=$DIRINSTALL/$iphreeqc" >> ~/.bashrc
    echo 'export PATH=$PATH:$IPHREEQC/lib:$IPHREEQC/include' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$IPHREEQC/lib' >> ~/.bashrc
    echo >> ~/.bashrc
}

## install phreeqcrm library ##
## https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
OPT_PHCRM(){
    echo
    echo " #### Install $phreeqcrm ####"
    #PAUSE 'Press [Enter] key to continue...'
    RM_DIR $DIRINSTALL/$phreeqcrm
    wget "ftp://brrftp.cr.usgs.gov/pub/charlton/phreeqcrm/$phreeqcrm.tar.gz"
    tar -xf "$phreeqcrm.tar.gz"
    cd $phreeqcrm || exit
    ./configure --prefix=$DIRINSTALL/$phreeqcrm  CC=$CC FC=$FC
    make
    $sudo make install
    echo "############ $phreeqcrm ############" >> ~/.bashrc
    echo "PHREEQCRM=$DIRINSTALL/$phreeqcrm" >> ~/.bashrc
    echo 'export PATH=$PATH:$PHREEQCRM/lib:$PHREEQCRM/include' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PHREEQCRM/lib' >> ~/.bashrc
    echo >> ~/.bashrc
}

## remove install directory ##
OPT_RM(){
    $sudo rm -rf $TMPDIRINSTALL
    echo " $TMPDIRINSTALL removed with success"
    echo
}

## remove install directory if exist ##
RM_DIR(){
    LINK_OR_DIR=$1
    if [ -d "$LINK_OR_DIR" ]; then
        read -r -p " Install dir $LINK_OR_DIR already exist. Are you sure you want to continue? [y/N] " response
        case "$response" in
          [yY])
            $sudo rm -rf "$LINK_OR_DIR"
            ;;
          *)
            exit 0
            ;;
        esac
    fi
}

# ----------- Main -------------- ##

## create install and tmp directory ##
TMPDIRINSTALL=~/.mohid
if [ ! -d "$DIRINSTALL" ]; then
    $sudo mkdir -p $DIRINSTALL
fi
if [ ! -d "$TMPDIRINSTALL" ]; then
    mkdir -p $TMPDIRINSTALL
fi

##
OS=$(cat /etc/*-release)
case $OS in
    *Ubuntu*)
        PACKMANAG=apt
        ;;
    *CentOS*)
        PACKMANAG=yum
        ;;
    *)
        echo
        echo " This script only work for Ubuntu and CentOS Linux distribuition"
        echo " Maybe you can manually adapt for your linux distribution..."
        echo
        exit 1
        ;;
esac

## install options ##
if [ $# -lt 1 ] ;then
  HELP
  exit 0
fi

cd $TMPDIRINSTALL || exit 1

case $1 in
    -h|-help|--help)
       HELP
    ;;
    -req|requirements)
       OPT_REQ
       exit 0
    ;;
    -zlib)
       OPT_ZLIB
       WARNING
    ;;
    -hdf5)
       OPT_HDF5
       WARNING
    ;;
    -nc|netcdf)
       OPT_NC
       WARNING
    ;;
    -ncf|netcdff)
       OPT_NCF
       WARNING
    ;;
    -proj4)
       OPT_PROJ4
       WARNING
    ;;
    -proj4f|proj4fortran)
       OPT_PROJ4F
       WARNING
    ;;
    -phqc|iphreeqc)
       OPT_IPHC
       WARNING
    ;;
    -phqcrm|phreeqcrm)
       OPT_PHCRM
       WARNING
    ;;
    -mpich)
       OPT_MPICH
       WARNING
    ;;
    -openmpi)
       OPT_OPENMPI
       WARNING
    ;;
    -rm)
       OPT_RM
       exit 0
    ;;
    *)
       echo
       echo -e " ${ERROR}: unrecognized command line option $1. Use  $0 -h  for help."
       exit 0
esac

exit 0

