# MOHID instalation in Linux machines #
---
This repository is an attempt to compile MOHID in a linux machine with Ubuntu or CentOS linux distro and intel compiler ifort.

## Minimum system requirements:
Install the following softwares:
* intel compiler C and Fortran
* zlib
* hdf5
* netcdf (optional)
* proj4 (optional)
* proj4fortran (optional)
* mpich (optional)
* mohid src

or
run ``./install_req.sh -opt1 -opt2 ...``
**Install intel compiler first!** `install_req.sh` will download and install all libraries required to compile MOHID. **You must respect the instalation order**. For more information or see all option available, run:
```
./install_req.sh -h
```
Update all paths in script to fit your setup.

## Compile
To compile MOHID Water and/or MOHID Land, you must compile first MOHID Base 1 and MOHID Base 2. This will create all .mod files used to compile MOHID.
```
./compile.sh -mb1 -mb2
```
To compile MOHID Water run:
```
./compile.sh -mw
```
To compile MOHID Land run:
```
./compile.sh -ml
```
For more information or see all option available, run:
```
./compile.sh -h
```
Update all paths in script to fit your setup.

## Notes: ##
* You must download and install intel compiler C and Fortran first, before run scripts. intel compiler free for 30 days. You can do it here: <https://software.intel.com/en-us/intel-parallel-studio-xe>
* If you are a Academic Researcher or a Student you can get intel compiler for free here: <https://software.intel.com/en-us/qualify-for-free-software>
* See <https://software.intel.com/en-us/get-started-with-parallel-studio-xe-for-linux> to learn more about what to install.
* After installed intel compiler C and fortran, add paths in .bashrc:
```
export INTEL=/opt/intel
export PATH=$PATH:$INTEL/bin:$INTEL/lib/intel64:$INTEL/include:$INTEL/include/intel64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL/lib/intel64
export CC=$INTEL/bin/icc
export CXX=$INTEL/bin/icpc
export FC=$INTEL/bin/ifort
export F9X=$INTEL/bin/ifort
export CFLAGS='-O3 -xHost -ip'
export CXXFLAGS='-O3 -xHost -ip'
export FCFLAGS='-O3 -xHost -ip'
```
* Update all paths in scripts
* MOHID must be compiled with the same compiler used in hdf5 and netcdf libraries
* If necessary use ext2uppercase.sh to uppercase all .f90 to .F90 extensions
* If necessary use dos2unix linux command:
` find . -name *.F90 -exec dos2unix {} \; `

## Help, Bugs, Feedback
If you need help with compile MOHID in linux, you can hang out by mail: <general@mohid.com> or consult our [wiki](http://wiki.mohid.com). You can also subscribe to our [forum](http://forum.mohid.com). To report bugs, please contact any developers. More information consult <http://www.mohid.com>

## License
GNU General Public License. See the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) web page for more information.