# Setup a MOHID simulation using MPI

## Windows

### Requirements

In order to run a simulation using MPI in MOHID you must fulfill the following requirements:

1. These instructions are specific of a MOHID application running under the MS Windows.
2. Remember this can only be applied if you have at least one nested domain.
3. Make sure you have an executable file which was compiled together with MPI libraries and is able to be run in parallel. Learn more on [compiling Mohid with MPI](build/compiling-mohid-with-mpi.md).
4. In order to take advantage of the [parallel processing](programming/parallel-processing.md) features make sure you have one processor available for each domain and sub-domains coupled in your simulation.
5. Make sure that you have installed MPICH 2.0.x (download [here](http://www.mcs.anl.gov/research/projects/mpich2/))in each PC available for the simulation. To make sure MPICH is installed correctly and working properly, go to MS Windows Task Manager and look for a service called **mdm.exe**. In addition you should also follow all the instructions provided during the MPICH installation. Make sure that you've added MPICH installation folder to the PATH environmental variable. Just go to Control Panel -> System -> Advanced -> Environmental Variables. In the "System Variables" select the PATH variable and press Edit. Then add to the end of the string "C:\\Program Files\\MPICH2\\bin", which is the path to where the mpirun.exe is installed. You need to install MPICH in every computer that you which to use. You may need to restart the computers after the installation, in order to get everything running.
6. Since October 2018, MOHID is prepared to run with Microsoft MPI (download Microsoft MPI v10 [here](https://www.microsoft.com/en-us/download/details.aspx?id=57467)).

#### Setup the simulation

1. Setup the simulation with nested domains exactly the same way as in a normal simulation. Use [MOHID GUI](../tools/gui/mohid-gui.md) for an easier setup. No other specific configuration is necessary at this point
2. Test the simulation by executing the model in a single machine
3. Go to the **\\exe** folder in the *father* domain. There you must find at least 2 files.
    1.  [nomfich.dat](../io/nomfich-dat-file.md)
    2.  [Tree.dat](../io/tree-dat.md)

#### MPI command line

The command line shown bellow is the simplest way to run Mohid using MPI:

```
Location of mpiexec -np "number of processors" Location of Mohid Water Exe with MPI

An example with 2 processors:

"C:\Program Files\Microsoft MPI\Bin\mpiexec" -np 2 "D:\....\MOHIDWater_x64_release_double_MPI.exe"
```

In the above example MohidWater is running with 2 processes, using only local machine's resources. In case that you are using Microsoft MPI you should modify the path to your executable accordingly (i.e. "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec".

To run MohidWater using multiple computers use something like this:

```
>mpiexec  -hosts 2 RICARDOM-VM 1 FRANK-FIXO 2 -map M:\\RICARDOM-VM\projects M:\aplica\MPI\FATHER\exe\MohidWater.MPI.exe
```

Both machines are required to have the same executable. Here a explicit map is made so that both machines know where the executable is. M: must be mounted and mpiexec mus be invoqued from the mount point. -map will create the mount in slave machines and unmount it at the end of the run.

#### Tree.dat configuration

Define the number of processors for each domain by adding : and number of processors for that level after the level's exe folder location.

Paths to the domains exe folders must be writen with relative paths (i.e. ..\\..\\)

In the example below we are using 6 processors, 1 for the upper level domain and 5 for the nested domain:

```
+..\..\Level1\exe
++..\..\Level1\Level2\exe : 5
```

### Requirements

1. Remember this can only be applied if you have at least one nested domain.
2. Make sure you have an executable file which was compiled together with MPI libraries and is able to be run in parallel. Learn more on [compiling Mohid with MPI](build/compiling-mohid-with-mpi.md).
3. In order to take advantage of the [parallel processing](programming/parallel-processing.md) features make sure you have one processor available for each domain and sub-domains coupled in your simulation.
4. Make sure that you have installed the smpd daemon. Refer to the linux instructions in [compiling Mohid with MPI](build/compiling-mohid-with-mpi.md).
5. Add in the Bathymetry data file the keyword HALOPOINTS : 3 in the domains whose bathymetries would be decomposed in MPI domains
