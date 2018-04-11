NETCDF: 
- UNIDATA supports Windows support via CMake. Download and install Cmake [4]. Cmake 
is used to control the software compilation process using simple platform and compiler independent 
configuration files, and generate native makefiles and workspaces that can be used in the compiler 
environment of your choice. It is used by UNIDATa and also HDF5 Group;
- The UNIDATA provides all the source code necessary to compile the NETCDF libraries and include files that 
work as na interface between the the NETCDF and fortran. You nedd to download the NetCDF-Fortran Releases
from the follow link [5]
- The UNIDATA provides all the the pre-built releases of NETCDF libraries and dll's for Windows (win32 or x64). 
  Need to download the version that corresponds to the NetCDF-Fortran Releases and install 
  the NETCDF pre-build for windows 

 The steps to compile netcdf are the follow (work for VS2015 + Intel 18 for win32 and x64):
  Step 1 - Extract all the files from the compress file provide by UNIDATA with the NetCDF-Fortran 
Releases maintaining the directories tree. 
  Step 2 - Create a new diretory inside the parente diretory where all NetCDF-Fortran 
Releases are present (e.g. test_paulo_hdf5.1.8.17_x64)
  Step 3 - Create inside of this new diretory the a bat file (e.g. Run.bat) with the follow 
 instructions in one code line
  "C:\Program Files\CMake\bin\cmake" 
   [-G"NMake Makefiles" -DCMAKE_BUILD_TYPE:STRING=RELEASE 
    -DCMAKE_Fortran_COMPILER:STRING="C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.2.18\windows\bin\intel64/ifort.exe" 
     -DBUILD_SHARED_LIBS=FALSE 
     -DBUILD_V2=FALSE -DENABLE_TESTS=FALSE 
     -DnetCDF_INCLUDE_DIR:STRING="D:\Projects\paulo\mohid\ExternalLibs\HDF5_1.8.17_x64\Include\fortran" 
      -DNETCDF_C_LIBRARY:STRING="C:\Program Files\netCDF 4.4.1.1\lib\netcdf.lib"] ..
  Step 5 - First path need to point to the cmake.exe; 
           Second path need to point to the ifort.exe of interest (fortran compiler);
           Third path need to point to the include path of the HDF5 version and platform (win32 or x64) of interest;
           Fourth path need to point to the netcdf library netcdf.lib (win32 or x64) placed in the pc by the installation 
           of the pre-built releases of NETCDF;
  Step 6 - In the parent directory there is file called CMakeLists.txt. Probably need to comment some instructions (I did). 
           In my case the cmake only generate the VS2015 solution when I comment the follow instructions 
                      #CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nc_def_opaque "" USE_NETCDF4)
                      #CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nccreate "" USE_NETCDF_V2)
                      #CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nc_set_log_level "" USE_LOGGING) 
                      #CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} oc_open "" BUILD_DAP)
                      #CHECK_LIBRARY_EXISTS(${NETCDF_C_LIBRARY} nc_use_parallel_enabled "" BUILD_PARALLEL)
  Step 7 - For details look at  https://github.com/Unidata/netcdf-fortran/issues/48  
  Step 8 - After executing the bat file a solution NC4F.sln is generate in the same place of the bat file. Open the solution
  Step 9 - In my case (VS2015 + Intel 18) it was missing in the solution to add in the c library ncfortran 
           "Properties\C/C++\General\Additional Include Directories" the path where the NETCDF installation of NETCDF pre-build
             for Windows place the include files (*.h) "C:\Program Files\netCDF 4.4.1.1\include"
  Step 10 - Need to change the "netcdf" name to "netcdf90" to be compatible with the MOHID instruction 
            for Windows (use ModuleNetcdf90)
  Step 11 - For the x64 platform the follow a configuration occure associate with the "target machine". Is define as X86 and  
            should be x64. Need to change in the "Properties\Library\Comand line"
  Step 12 - Save the solution and build the solution. Copy the netcdff.lib and all the *.mod to a diretory that you want to use
            after in the MOHID build (e.g. ..\..\mohid\ExternalLibs\Netcdf_4.4.1.1\VS2015\x64);
  
 
 The steps to compile MOHID with the Netcdf include/library files are:
  Step 1 - In the MOHID solution configuration is necessary to add in the 
           "Properties\PreProcessor\Additional Include Directories" 
           the path diretory where the include files of NetCDF for the platform in question (win32 or x64) are placed      
            (e.g. ..\..\mohid\ExternalLibs\Netcdf_4.4.1.1\VS2015\x64)
  Step 2 - Additionaly it is necessary to add in the 
           "Linker\General\Additional Library Directories" the path diretory where the library files of NetCDf for the 
           platform in question (win32 or x64) are placed (e.g. ..\..\mohid\ExternalLibs\Netcdf_4.4.1.1\VS2015\x64)
  Step 3  - Finally it is necessary to add the Netcdf libraries of interest (e.g. netcdf.lib and netcdff.lib)                 
           in "Linker\Input\Additional Dependencies"