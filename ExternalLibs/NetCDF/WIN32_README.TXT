Prebuilt NetCDF-3.6.1 Beta1 WIN32_README

Origional Date   : May 25, 2006
Revised   Date   : Jan 7,2007
Revised   Date   : Jul 19, 2007
Revised   Date   : Dec 12, 2007 'add F90 interface to NetCDF
Author           : Dr. YAN Haoming
Email            : wfllib@gmail.com
Download page    : http://yhmyhm.googlepages.com/fortran.htm 
                   http://www.whigg.ac.cn/yanhm/fortran.htm

OS               : Windows2000 / XP
Support software : CVF6.*        --Compaq Visual Fortran compiler version 6.*
                   IVF9.*, 10.*  --Intel  Visual Fortran compiler version 9.*, 10.*
NetCDF homepage  : http://www.unidata.ucar.edu/software/netcdf/

*************************************************************************

Contents:

*************************************************************************


0. Tree of the prebuilt NetCDF-3.6.1 Beta1 Win32 file "NetCDF3.6.1Beta1.Win32.for.CVF.IVF.rar"
1. Using pre-built binaries
2. Using NetCDF libraries in Compaq Visual Fortran (CVF) 6.* and Intel "Visual" Fortran(IVF) 9.*, 10.*
3. Using NetCDF Fortran90 interface in CVF6.6 and IVF 9.* or 10.*



*************************************************************************
0. Tree of the prebuilt NetCDF-3.6.1 Beta1 Win32 file "NetCDF3.6.1Beta1.Win32.for.CVF.IVF.rar"

    NETCDF
    ©¦  
    ©¦  netcdf.pdf                 !document for NetCDF
    ©¦  netcdf-f77.pdf             !f77 document for NetCDF
    ©¦  netcdf-f90.pdf             !f90 document for NetCDF
    ©¦  WIN32_README.TXT
    ©¦  
    ©À©¤bin
    ©¦      ncdump.exe
    ©¦      ncgen.exe
    ©¦      netcdf.dll
    ©¦      
    ©À©¤F90                         !Fortran90 interface directory for NetCDF
    ©¦  ©¦  netcdf_test.f90
    ©¦  ©¦  
    ©¦  ©À©¤CVF
    ©¦  ©¦  ©À©¤include
    ©¦  ©¦  ©¦      netcdf90.mod
    ©¦  ©¦  ©¦      typesizes.f90
    ©¦  ©¦  ©¦      typesizes.mod
    ©¦  ©¦  ©¦      
    ©¦  ©¦  ©¸©¤lib
    ©¦  ©¦        netcdf90.lib
    ©¦  ©¦          
    ©¦  ©¸©¤IVF
    ©¦      ©À©¤include
    ©¦      ©¦      netcdf90.mod
    ©¦      ©¦      typesizes.f90
    ©¦      ©¦      typesizes.mod
    ©¦      ©¦      
    ©¦      ©¸©¤lib
    ©¦             netcdf90.lib
    ©¦              
    ©À©¤include                    !Fortran77 include files for NetCDF
    ©¦       netcdf.h
    ©¦       netcdf.inc
    ©¦          
    ©¸©¤lib                        !Fortran77 lib files for NetCDF
            netcdf.lib
            netcdfs.lib




*************************************************************************
1. Using pre-built binaries

I have pre-built Win32 binary versions of the netcdf dll and static library. 
DLL versions of ncgen.exe and ncdump.exe are download from NetCDF home page
(http://www.unidata.ucar.edu/software/netcdf/).  
(Note: we do NOT have a C++ interface here).

To use:

 To use the Dynamic Library (shared) version of the NetCDF library:
 Put below files in a directory which is in your system search PATH:
	netcdf.dll		library dll
	ncgen.exe		uses the dll
	ncdump.exe		uses the dll




*************************************************************************
2. Using NetCDF libraries in Compaq Visual Fortran (CVF) 6.* and Intel "Visual" Fortran (IVF) 9.*, 10.*

I had revised the Fortran interface of NetCDF, 
so it can be used in both CVF and IVF with no additional modification. 

To use NetCDF in both CVF6.* and IVF9.*, 10.*:

a. copy 'netcdf.inc' to the include file directory of CVF or IVF. 
   e.g., c:\Program Files\Microsoft Visual Studio\df98\inclde\    
   or    c:\Program Files\Intel\Compiler\Fortran\9.*\IA32\include\
   or    c:\Program Files\Intel\Compiler\Fortran\10.*\IA32\include\

b. copy 'netcdf.lib' and 'netcdfs.lib' to the library file directory of CVF or IVF. 
   e.g., c:\Program Files\Microsoft Visual Studio\df98\lib\    
   or    c:\Program Files\Intel\Compiler\Fortran\9.*\IA32\lib\ 
   or    c:\Program Files\Intel\Compiler\Fortran\10.*\IA32\lib\

c. copy 'netcdf.dll', 'ncgen.exe', 'ncdump.exe' to 
   c:\Program Files\winnt\system32\  
   or  c:\Program Files\windows\system32\

d. add " include 'netcdf.inc' " in the first line of your Fortran main program, 
   then you can use NetCDF library freely.




*************************************************************************
3. Using NetCDF Fortran90 interface in CVF6.6 and IVF 9.* or 10.*
I had successfully build the NetCDF Fortran90 interface and library in CVF6.6 and IVF9.* or 10.*.
The original Fortran90 interface source codes of NetCDF was provided by the NetCDF site, you also can 
find it in the NetCDF source code under f90 directory. The author is Robert Pincus, thank you (see 
more details in typesizes.f90 comments). I revised some of these source codes and build
the prebuilt binaries.

Note, the NetCDF Fortran90 interface based on the NetCDF Fortran77 library, so you must install NetCDF77 library
correctly first in the contents 2.

To use NetCDF Fortran90 functions in both CVF6.6 and IVF 9.* or 10.*

a.   copy the files in ..\netcdf\f90\cvf\include\ to the suitable directory (reference contents 2) for CVF6.6
  or copy the files in ..\netcdf\f90\ivf\include\ to the suitable directory for IVF

b.   copy the file in ..\netcdf\f90\cvf\lib\ to the suitable directory for CVF6.6
  or copy the file in ..\netcdf\f90\ivf\lib\ to the suitable directory for IVF

C. add " use netcdf90 " in the first line of your Fortran main program
   add " include netcdf.inc " after the first line or at the suitable position (if you have many " use ..." lines)

d. run ..\netcdf\f90\netcdf_test.f90 to test if Fortran90 interface run smoothly. Enjoy the NetCDF now.





*************************************************************************

Good luck!

Dr. YAN Haoming

Any questions and suggestions, please send email to: wfllib@gmail.com


Dec. 12, 2007.



