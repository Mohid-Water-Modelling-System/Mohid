# MPIConfig.dat

The MPIConfig.dat file is no longer used in MPICH 1.2.x. You must write a command line with every option you intend to use.

---

The **MPIConfig.dat** is the setup file that distributes a [MOHID Water](mohid-water.md) simulation in several processors associating each nested domain to a processor. Below an example explanation will be given, considering that:

1. There is a simulation with a nested domain called **Rias** nested into a father domain called **Galicia**
2. There are 2 computers connected through a network called PC1 and PC2, each with one processor available
3. The simulation data files are all located in PC1
4. The full path to the simulation father domain working folder is **\\\\PC1\\Applications\\Galicia\\exe**
5. PC1 has a virtual Mapped Network Drive with the letter **T:\\** redirecting to the following full path **\\\\PC1\\Applications\\**.
6. The whole simulation path is relative to this **T:\\** drive.
7. The MOHID executable compiled with MPI options is located in PC1
8. The operative system being used is MS Windows XP Professional and the MPICH version installed is 1.1

The **MPIConfig.dat** must be configured like shown below. The construction of this file can be also done via [MOHID GUI](../tools/gui/mohid-gui.md).

**Line 1** - the first line contains a keyword **exe** containing the full path to the MOHID executable compiled with MPI options

```
exe \\PC1\EXE\Mohid_MPI.exe
```

**Line 2** - this line contains a keyword **dir** containing path to the father domain working folder (**\\exe**).

```
dir T:\Galicia\exe
```

**Line 3** - this line contains a keyword **map** containing the letter of the mapped drive **T:** followed by the full path to the father domain working folder (**\\exe**).

```
map T:\\PC1\Applications\Galicia\exe
```

**Line 4** - this line contains a keyword **hosts**. Below this line the definition of each computer specifying its name (e.g. PC1) and the number of processes to be launched in this computer.

```
hosts
PC1 1
PC2 1
```

For example, if PC1 is a computer carrying 2 processors the following definition could be written:

```
hosts
PC1 2
```
