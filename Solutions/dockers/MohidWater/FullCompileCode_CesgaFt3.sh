# Versao FT3 - Meteogalicia

# Load ifort
module load intel
module load impi
echo "Loaded ifort..."

# instalar zlib
chmod +x ../externalSoftware/install_zlib.sh && cd ../externalSoftware/ && ./install_zlib.sh

# install hdf5
chmod +x install_hdf5.sh && ./install_hdf5.sh


## install libraries path
cd ../MohidWater/

DIRINSTALL=$HOME/apps_intel
zlib='zlib-1.2.11'
hdf5='hdf5-1.8.17'

RW=${PWD%}
OUTPUT=$RW/bin
mkdir $OUTPUT
echo "Working Folder ... $RW"

echo 'export PATH=$PATH:$ZLIB/lib:$ZLIB/include' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ZLIB/lib' >> ~/.bashrc
echo "HDF5=$DIRINSTALL/$hdf5" >> ~/.bashrc
echo 'export PATH=$PATH:$HDF5/bin:$HDF5/lib:$HDF5/include' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5/lib' >> ~/.bashrc

## -g -traceback only usefull for code debuging
SRCPATH=../../../Software
MPI=/opt/cesga/2020/software/Compiler/intel/2021.3.0/impi/2021.3.0/mpi/latest
MPIINC=$MPI/include
MPILIB=$MPI/lib
MPI_LIB="-lmpi -lmpifort"

ifort -O3 -w -cpp -real-size 64 -convert little_endian -fPIC -heap-arrays 64 -fp-model source -xHost -ip -fpe0 -fpp -D_NO_NETCDF -D_INCREASE_MAXINSTANCES_EXTRA -D_LAGRANGIAN_GLOBAL_ -D_USE_NIX -D_STACK_LIMITS_ -D_BIG_LINE_LENGTH -D_INCREASE_MAXINSTANCES -D_USE_MPI $SRCPATH/MOHIDBase1/ModuleGlobalData.F90  $SRCPATH/MOHIDBase1/ModuleLUD.F90  $SRCPATH/MOHIDBase1/ModuleTriangulation.F90  $SRCPATH/MOHIDBase1/ModuleTime.F90  $SRCPATH/MOHIDBase1/ModuleEnterData.F90  $SRCPATH/MOHIDBase1/ModuleWWTPQ.F90  $SRCPATH/MOHIDBase1/ModuleStopWatch.F90  $SRCPATH/MOHIDBase1/ModuleFunctions.F90  $SRCPATH/MOHIDBase1/ModuleMacroAlgae.F90  $SRCPATH/MOHIDBase1/ModuleWaterQuality.F90  $SRCPATH/MOHIDBase1/ModuleSedimentQuality.F90  $SRCPATH/MOHIDBase1/ModuleHydroIntegration.F90  $SRCPATH/MOHIDBase1/ModuleSeagrassWaterInteraction.F90  $SRCPATH/MOHIDBase1/ModuleHDF5.F90  $SRCPATH/MOHIDBase1/ModuleHDF5_OO.F90  $SRCPATH/MOHIDBase1/ModuleSeagrassSedimInteraction.F90  $SRCPATH/MOHIDBase1/ModuleLife.F90  $SRCPATH/MOHIDBase1/ModuleCEQUALW2.F90  $SRCPATH/MOHIDBase1/ModuleBenthos.F90  $SRCPATH/MOHIDBase1/ModuleDrawing.F90  $SRCPATH/MOHIDBase1/ModuleProfile.F90  $SRCPATH/MOHIDBase1/ModuleBivalve.F90  $SRCPATH/MOHIDBase1/ModuleBenthicEcology.F90 $SRCPATH/MOHIDBase1/ModuleInterface.F90  $SRCPATH/MOHIDBase1/ModuleTimeSerie.F90  $SRCPATH/MOHIDBase1/ModuleDischarges.F90  $SRCPATH/MOHIDBase1/ModuleLightExtinction.F90  $SRCPATH/MOHIDBase1/ModuleDrainageNetwork.F90  $SRCPATH/MOHIDBase1/ModuleMPImanagement.F90  $SRCPATH/MOHIDBase1/ModuleTask2000.F $SRCPATH/MOHIDBase2/ModuleHorizontalGrid.F90  $SRCPATH/MOHIDBase2/ModuleStatistic.F90  $SRCPATH/MOHIDBase2/ModuleGridData.F90  $SRCPATH/MOHIDBase2/ModuleBasinGeometry.F90  $SRCPATH/MOHIDBase2/ModuleHorizontalMap.F90  $SRCPATH/MOHIDBase2/ModuleBoxDif.F90  $SRCPATH/MOHIDBase2/ModuleGeometry.F90  $SRCPATH/MOHIDBase2/ModuleMap.F90  $SRCPATH/MOHIDBase2/ModuleAdvectionDiffusion.F90  $SRCPATH/MOHIDBase2/ModuleInterpolation.F90    $SRCPATH/MOHIDBase2/ModuleTwoWay.F90  $SRCPATH/MOHIDBase2/ModuleField4D.F90  $SRCPATH/MOHIDBase2/ModuleFillMatrix.F90  $SRCPATH/MOHIDBase2/ModuleChainReactions.F90  $SRCPATH/MOHIDBase2/ModuleAtmosphere.F90 $SRCPATH/MOHIDWater/ModuleTurbine.F90  $SRCPATH/MOHIDWater/ModuleGOTM.F90  $SRCPATH/MOHIDWater/ModuleTurbGOTM.F90  $SRCPATH/MOHIDWater/ModuleFreeVerticalMovement.F90  $SRCPATH/MOHIDWater/ModuleToga.F90  $SRCPATH/MOHIDWater/ModuleGauge.F90  $SRCPATH/MOHIDWater/ModuleOil.F90  $SRCPATH/MOHIDWater/ModuleOil_0D.F90  $SRCPATH/MOHIDWater/ModuleHNS.F90  $SRCPATH/MOHIDWater/ModuleOpenBoundary.F90  $SRCPATH/MOHIDWater/ModuleTurbulence.F90  $SRCPATH/MOHIDWater/ModuleHydrodynamicFile.F90  $SRCPATH/MOHIDWater/ModuleAssimilation.F90  $SRCPATH/MOHIDWater/ModuleWaves.F90  $SRCPATH/MOHIDWater/ModuleJet.F90  $SRCPATH/MOHIDWater/ModuleSand.F90  $SRCPATH/MOHIDWater/ModuleConsolidation.F90  $SRCPATH/MOHIDWater/ModuleHydrodynamic.F90  $SRCPATH/MOHIDWater/ModuleWaterProperties.F90  $SRCPATH/MOHIDWater/ModuleLagrangian.F90  $SRCPATH/MOHIDWater/ModuleLagrangianGlobal.F90  $SRCPATH/MOHIDWater/ModuleSedimentProperties.F90  $SRCPATH/MOHIDWater/ModuleSediment.F90  $SRCPATH/MOHIDWater/ModuleInterfaceSedimentWater.F90  $SRCPATH/MOHIDWater/ModuleInterfaceWaterAir.F90  $SRCPATH/MOHIDWater/ModuleModel.F90 $SRCPATH/MOHIDWater/Main.F90 -I$HDF5/include -L$HDF5/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -ldl -I$ZLIB/include -L$ZLIB/lib -lz -lm -I${MPIINC} -L${MPILIB} $MPI_LIB -lmpi -lmpifort -o MohidWater.exe
