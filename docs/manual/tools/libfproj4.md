# Libfproj4

libfproj4 is a Fortran wrapper for [PROJ4](proj4.md) library. It was developed as of the [GLIMMER](http://www.cbcb.umd.edu/software/glimmer/) project and can be found in [libfproj4-1.0.tar.gz](http://forge.nesc.ac.uk/download.php/10/libfproj4-1.0.tar.gz). This wrapper is necessary if we want to call PROJ4 conversion subroutines (in C) from a Fortran program. It is used in MOHID by the ConvertToHDF5 tool to be able to deal with different projections from other models.

## Compilation (Static Library)

libfproj4 it's a library build on top of [PROJ4](proj4.md). It can se used dynamically or statically. Here, only the static compilation is explained, as the instructions for dynamic compilation are included in the library source code. The first steps in compiling are:

- Download libfproj4 source code from [

$$
1
$$

](http://forge.nesc.ac.uk/download.php/10/libfproj4-1.0.tar.gz) or \\\\kepler\\Software\\PROJ4\\libfproj4-1.0.tar.
- Download proj4 source code from [

$$
2
$$

](http://download.osgeo.org/proj/proj-4.7.0.tar.gz) or \\\\kepler\\Software\\PROJ4\\proj-4.7.0.tar.
- Extract all contents to a place in your computer, from now referred to as ROOT\_DIR.

### Windows

#### Compaq Visual Studio

(Already compiled in \\\\kepler\\Software\\PROJ4\\libfproj4.lib)

1. New project W32 static lib
2. Add fort-proj.c (from libfproj4-1.0.tar.gz) and all \*.c files (from proj4 source), except all that have main function : css2css, geod, nad2bin, nad2nad, proj.c.
3. Project Settings / C/C++ / General / Preprocessor definitions : VISUAL\_CPLUSPLUS ou IFORT
4. Project Settings / C/C++ / Preprocessor / include : whatever\\proj-4.6.0\\src (the folder that has projects.h, because of directives #include <projects.h>)

#### Visual Studio 2008 with Intel Fortran

(Already compiled in [

$$
3
$$

](http://mohid.codeplex.com) ExternalLibs/Proj4)

1. Open Visual Studio and select File/New/Project/Visual C++/W32/Win32 Console Application. Name it libfproj4, and select location as $ROOT\_DIR\\libfproj4-1.0. Click ok.
2. In the wizard, click Next and then choose "Console application" and "Empty project" in the following window. Then click Finish.
3. Go to Project / Properties / Configuration Properties / Configuration Type and select "Static Library (.lib)". Click ok. (we don't choose this type in the project creation because if so VS2008 automatic header files that we don't want).
4. Go to Project / Add Existing Item. Add $ROOT\_DIR\\libfproj4-1.0\\fort-proj.c.
5. Go to Project / Add Existing Item. Add all \*.c files in $ROOT\_DIR\\proj-4.7.0\\src except the ones that have the main function: cs2cs.c, geod.c, nad2bin.c, nad2nad.c, proj.c.
6. Rename proj\_config.h.in to proj\_config.h in ROOT\_DIR\\proj-4.7.0\\src
7. Go to Project / libfproj4 Properties:
    -   C/C++ / General / Additional Include Directories: ..\\..\\..\\proj-4.7.0\\src (where the \*.h files are)
    -   C/C++ / Preprocessor / Preprocessor Definitions: IFORT
    -   C/C++ / Code Generation / Runtime Library : Multi-threaded Debug (/MTd) (para versao Debug)
    -   C/C++ / Code Generation / Runtime Library : Multi-threaded (/MT) (para versao Release)
8. Build Solution

### Linux

Prepare files and folders:

```
cd $ROOT_DIR
tar -zxvf proj-4.7.0.tar.gz
tar -zxvf libfproj4-1.0.tar.gz
mkdir src
cp libfproj4-1.0/fort-proj.c src/
cp libfproj4-1.0/cfortran.h src/
cp -r proj-4.7.0/src/*.c src/
cp -r proj-4.7.0/src/*.h src/
cp proj-4.7.0/src/proj_config.h.in src/proj_config.h
```

Compile:

```
cd $ROOT_DIR/src
rm -rf cs2cs.c geod.c nad2nad.c nad2bin.c proj.c  # with main functions
icc -DIFORT -c *.c
ar rc libfproj4.a *.o
```

Install:

```
cd $ROOT_DIR
mkdir lib
mkdir include
cp src/libfproj4.a lib/
cp libfproj4-1.0/proj4.inc include
```

## Usage

### Windows

1. Create $ROOT\_DIR\\teste
2. Copy test-proj.f90, proj4.f90 and proj4.inc from $ROOT\_DIR\\libfproj4-1.0 to $ROOT\_DIR\\teste.
3. Copy libfproj4.lib from $ROOT\_DIR\\libfproj4-1.0\\libfproj4\\Debug to $ROOT\_DIR\\teste.
4. Open Visual Studio and select File/New/Project/Intel(R) Visual Fortran / Console Application / Empty Project. Name it teste and select location as $ROOT\_DIR\\teste. Click ok.
5. Project / Add Existing Item: Add test-proj.f90 and proj4.f90
6. Project / Properties:
    -   Linker / General / Additional Library Directories: ..\\..\\..\\teste (or where libfproj4.lib file is).
    -   Linker / Input / Additional Dependencies: libfproj4.lib

### Linux

```
cd $ROOT_DIR/
mkdir teste
cd teste
cp ../libfproj4-1.0/test-proj.f90 .
cp ../libfproj4-1.0/proj4.f90 .
cp ../libfproj4-1.0/proj4.inc .
cp ../src/libfproj4.a .
ifort -c proj4.f90      #-I../libfproj4-1.0       # needs proj4.inc
ifort -c test-proj.f90
ifort proj4.o test-proj.o libfproj4.a -o test-proj.exe
./test-proj.exe
```

## InternalLinks

- [PROJ4](proj4.md)
- [ConvertToHDF5](converttohdf5.md)

## External Links

- [libfproj4-1.0.tar.gz](http://forge.nesc.ac.uk/download.php/10/libfproj4-1.0.tar.gz) (deprecated)
- [proj4-fortran](https://github.com/mhagdorn/proj4-fortran) (new link)
