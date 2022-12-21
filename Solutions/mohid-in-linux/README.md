# MOHID instalation in Linux machines #
---
How to compile MOHID in a linux machine with Ubuntu linux distro and intel compiler ifort.

## Install Intel® oneAPI Base Toolkit for Linux

### Configure CPU System

Install CMake*, pkg-config, and GNU* (Ubuntu*):

`sudo apt update`

`sudo apt -y install cmake pkg-config build-essential`

Verify the installation by displaying the installation location with this command:

`which cmake pkg-config make gcc g++`

One or more of these locations will display:
/usr/bin/cmake
/usr/bin/pkg-config
/usr/bin/make
/usr/bin/gcc
/usr/bin/g++

### Pre-installation Steps

Download the key to system keyring:

`wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
| gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null`


Add signed entry to apt sources and configure the APT client to use Intel repository:

`echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list`

Update packages list and repository index:

`sudo apt update`

### Instalation

`sudo apt install intel-basekit`

`sudo apt install intel-hpckit`

Verify the installation:

`source /opt/intel/oneapi/setvars.sh`

`ifort --version`

`icc --version`

### After Installation

After you have installed the oneAPI kits, they can be activated using the following command written in the ~/.bashrc file.

`vim ~/.bashrc`

`. /opt/intel/oneapi/setvars.sh`

### References

https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-base-linux/top/before-you-begin.html

https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/apt.html


## Install Git and clone MOHID repository

`sudo apt install git-all`

`git clone https://github.com/Mohid-Water-Modelling-System/Mohid.git`


## Install required libraries

`cd /home/mohid/Mohid/Solutions/mohid-in-linux`

`./install_req.sh -zlib`

`./install_req.sh -hdf5`

`sudo apt install m4`

`./install_req.sh -nc`

`./install_req.sh -ncf`

**You must respect the instalation order**. 

*To see help menu run:*

./install_req.sh -h

*After installing each library you can refresh the current shell environment with the command:*

`source ~/.bashrc`


## Compile MOHID

To compile MOHID Water and/or MOHID Land, you must compile first MOHID Base 1 and MOHID Base 2. This will create all .mod files used to compile MOHID.

`./compile_mohid.sh -mb1`

`./compile_mohid.sh -mb2`

To compile MOHID Water run:

`./compile_mohid.sh -mw`

To compile MOHID Land run:

`./compile.sh -ml`

For more information or see all option available, run:

`./compile.sh -h`


## Test the installation ##
Inside test directory there are two test cases, one for mohid water and another one for mohid land. Just execute `mohid-in-linux/test/mohidwater/25m_deep/exe/MohidWater.exe` or `mohid-in-linux/test/mohidland/schematicWatershed/exe/MohidLand.exe` to see if there's no errors.

If you see `Program Mohid Water successfully terminated` or `Program Mohid Land successfully terminated`, your installation was sucessfull.

If you want to run your own project, make a link for `bin/MohidWater.exe` or `bin/MohidLand.exe`

```
ln -s <MOHID_INSTALLATION_PATH>/Mohid/Solutions/mohid-in-linux/bin/MohidWater.exe <YOUR_PROJECT_PATH>/exe/MohidWater.exe
```

## Notes ##
* Update all paths in scripts
* MOHID must be compiled with the same compiler used in hdf5 and netcdf libraries
* If necessary use ext2uppercase.sh to uppercase all .f90 to .F90 extensions
* If necessary use dos2unix linux command:
` find . -name *.F90 -exec dos2unix {} \; `
* Intel MPI commonly changes the value of MPI_TAG_UB which sometimes results in errors like this:

  *Fatal error in PMPI_Send: Invalid tag, error stack:
  PMPI_Send(159): MPI_Send(buf=0x7ffdad70954c, count=1, MPI_INTEGER, dest=0, tag=999003, MPI_COMM_WORLD) failed
  PMPI_Send(97).: Invalid tag, value is 999003*

  It's possible to overcome this error by exporting the following environment variables:

  export MPIR_CVAR_CH4_OFI_TAG_BITS=31

  export MPIR_CVAR_CH4_OFI_RANK_BITS=8

## Help, Bugs, Feedback ##
If you need help with compile MOHID in linux, you can hang out by mail: <general@mohid.com> or consult our [wiki](http://wiki.mohid.com). You can also subscribe to our [forum](http://forum.mohid.com). To report bugs, please contact any developers. More information consult <http://www.mohid.com>

## License ##
GNU General Public License. See the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) web page for more information.
