# SeaDAS

[SeaDAS](http://oceancolor.gsfc.nasa.gov/seadas/) (SeaWiFS Data Analysis System) is a sattelite image processing package developed by [GSFC/NASA](http://www.gsfc.nasa.gov/).

It is developed for a UNIX/LINUX platform, so, for Windows users the best choice is to create a [VMWare](../archive/linux/vmware.md) virtual machine and install SeaDAS, or to use an existing linux machine and install it.

It has a GUI and a command line interface. To process a large set of files, the best choice is to use the command line(or terminal/SSH) and use Bash or Seadas/IDL scripts.

SeaDAS command line input in a terminal/SSH:

```
 >Seadas -em -b input.cmd
```

input.cmd is text file with SeaDAS(IDL) commands, for example:

```
msl12,IFILE='../l1/S1998011121940.L1A_MLAC.x.hdf',
OFILE1='../l1/l2_nonir/S1998011121940.L2.hdf',par='ns_nonir.par',metsearch=1,ozonesearch=1,sstsearch=1
```

## Links
