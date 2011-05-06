This directory contains the binary (release) distribution of 
HDF4.2.5 that was compiled on an Intel PC running Windows XP, 
using Visual Studio 2008/Intel Fortran 10.1.  It was built with the 
following options: 
	-- C/Fortran libraries, both static and shared
	-- SZIP (encoder enabled), ZLIB, and JPEG
	-- Static HDF4 tools

The contents of this directory are:

    COPYING             - Copyright notice
    INSTALL_Windows.txt - Install instructions from the source code.
                          Section II discusses how to compile an
                          application
    README.txt          - This file
    RELEASE.txt         - Detailed information regarding this release
    bin\                - HDF4 static Utilities
    dll\                - HDF4 dlls 
    examples\           - HDF4 examples
    include\            - HDF4 include files
    lib\                - HDF4 libraries 
  

The binaries were built with ZLIB compression enabled (zlib 1.2.3). Therefore 
you MUST link with the GNU ZLIB library when compiling with these binaries. We 
provide the pre-compiled binary distribution for Windows for ZLIB 1.2.3 from 
our ftp server at:

    ftp://ftp.hdfgroup.org/lib-external/zlib/1.2/bin/windows

The binaries were built with JPEG Library (jpeg 6b). Therefore 
you MUST link with the JPEG library when compiling with these binaries. We 
provide the pre-compiled binary distribution for Windows for JPEG 6b from 
our ftp server at:

    ftp://ftp.hdfgroup.org/lib-external/jpeg/bin/windows

These binaries were also built with SZIP 2.1 compression (encoder ENABLED).  
When compiling with these binaries, you must link with the SZIP library.  
Information regarding SZIP can be found at:

   http://www.hdfgroup.org/doc_resource/SZIP/

You can obtain the SZIP source code and binaries from our ftp server at:
  
   ftp://ftp.hdfgroup.org/lib-external/szip/2.1/bin/windows/

Source code can be found on the THG ftp server in:

   ftp://ftp.hdfgroup.org/HDF/HDF_Current/src	

Please send questions, comments, and suggestions to:

    help@hdfgroup.org

