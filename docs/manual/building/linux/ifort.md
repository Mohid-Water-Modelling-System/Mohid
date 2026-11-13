# Ifort

**Ifort** is the intel fortran compiler. Intel bought and continued the former compaq Fortran, which in turn was originally known as DEC Fortran and used in VAX machines.

## Installing

### Windows

### \*nix

In order to install Ifort für \*nix, you need to do the following steps:

1. get the correct Ifort sources from 

$$
1
$$

 (registration required)
2. make sure, the required libraries are installed: *g++ libstdc++5 zlib1g-bin zlib1g-dev libc6* (Note: the names may differ depending on your distribution. The names above apply to Debian Sid)
3. compile the Ifort sources: un-tar the sources and execute the installation script. With Debian, you need to manually convert the rpms (see below).
4. integrate the compiler into your working environment:
    1.  Symlink the installation version, to avoid having to edit all the config files again, once you install a new compiler version.
    2.  Add the *source* lines to your bashrc (see below).

**Debian specific installation notes**: With Debian (Sid), the shipped *install.sh* script from the Ifort sources does not work 

$$
2
$$

,

$$
3
$$

. You need to execute the following commands as root (VERSION is the version number):

```
aptitude install alien
tar xzvf l_fc_p_VERSION.tar.gz
cd l_fc_VERSION/data/
alien intel-ifortVERSION-VERSION-1.i386.rpm
alien intel-iidbVERSION-VERSION-1.i386.rpm
dpkg -i --force-overwrite *.deb
```

Note: The resources below mention the additional command *./install.sh -finalize-only*, executed in the source folder. This parameter belongs to versions < 10 and is not supported by the install.sh script any more.

Now, replace every occurance of *<INSTALLDIR>* in the newly created files with the path to the Ifort installation. The files to edit are:

```
/opt/intel/fc/VERSION/bin/ifort
/opt/intel/fc/VERSION/bin/ifort.cfg
/opt/intel/fc/VERSION/bin/ifc
/opt/intel/fc/VERSION/bin/ifc.cfg
/opt/intel/fc/VERSION/bin/ifortvars.csh
/opt/intel/fc/VERSION/bin/ifortvars.sh
/opt/intel/idb/VERSION/bin/idbvars.csh
/opt/intel/idb/VERSION/bin/idbvars.sh
```

You can check if there are other files with the following command:

```
grep -r INSTALLDIR /opt/intel/
```

**Links and *source* lines**: The path to the compiler executables includes the version. Hard coding it has the disadvantage that you need to edit every file, once you install a new version. Therefore, link the current compiler version to *current* and add the virtual link in all of your path lines:

```
cd /opt/intel/fc/
ln -s ./VERSION current
cd ../idb/
ln -s ./VERSION current
```

Add the following lines to your *.bashrc*, or use the global */etc/bashrc* (Debian: */etc/bash.bashrc*):

```
source /opt/intel/fc/current/bin/ifortvars.sh
source /opt/intel/idb/current/bin/idbvars.sh
```

**Resources**: 

$$
[| 1](http://www.intel.com/cd/software/products/asmo-na/eng/219771.htm)
$$

 

$$
[| 2](http://groups.google.de/group/comp.lang.fortran/browse_thread/thread/d5efc992211e517b/25331fb8ea8ec01a?lnk=st&q=ifortvars%20%20%20%20.sh+debian+installdir&rnum=1&hl=en#25331fb8ea8ec01a)
$$

 

$$
[| 3](http://groups.google.de/group/linux.debian.user/browse_thread/thread/9076eb8ca951f0ce/6408013d92230ca1?lnk=st&q=Debian.te%20%20%20%20sting+Intel.Fortran&rnum=1&hl=en#6408013d92230ca1)
$$



## Tools for linking

### [dumpbin](../windows/dumpbin.md)

For example, here's an extract of the external references of ftest.obj after compiling with **ifort**:

```
ifort /c /nologo /Zi /Od ftest.for
dumpbin /SYMBOLS ftest.obj
```

```
094 00000000 UNDEF  notype ()    External     | _NCDREN@20
095 00000000 UNDEF  notype ()    External     | _NCVID@16
096 00000000 UNDEF  notype ()    External     | _NCVREN@20
097 00000000 UNDEF  notype ()    External     | _NCAREN@28
098 00000000 UNDEF  notype ()    External     | _NCVDEF@28
099 00000000 UNDEF  notype ()    External     | _NCVGT@24
09A 00000000 UNDEF  notype ()    External     | _NCVGTC@32
09B 00000000 UNDEF  notype ()    External     | _NCVGT1@20
09C 00000000 UNDEF  notype ()    External     | _NCVG1C@24
09D 00000000 UNDEF  notype ()    External     | _NCVPT@24
09E 00000000 UNDEF  notype ()    External     | _NCVPTC@32
09F 00000000 UNDEF  notype ()    External     | _NCVPT1@20
0A0 00000000 UNDEF  notype ()    External     | _NCVP1C@24
```

## Troubleshoots

### Windows

- Q: I get this error. What should I do?

```
Fatal error cannot open "ifconsol"
ifort: error: problem during multi-file optimization compilation (code 1)
```

A: This error happens because some environment variables need to be declared. To solve this simply open from the intel fortran folder in your *Start->Programs* the *Build Environment for applications* and run make from there.

### \*nix

## Help screen

Here's the help screen for windows:

```
ifort /?
```

```
                     Intel(R) Fortran Compiler Help
                     ==============================

usage: ifort [options] file1 [file2 ...] [/link linker_options]

  where options represents zero or more compiler options

        fileN is a Fortran source (.f .for .ftn .f90 .fpp .i .i90),
            assembly (.asm), object (.obj), static library (.lib), or
            other linkable file

        linker_options represents zero or more linker options

Notes
-----
1. Many FL32 options are supported; a warning is printed for unsupported
   options.

2. Intel Fortran compiler options may be placed in your ifort.cfg file.

                       Compiler Option List
                       --------------------
Performance
-----------
 /O1    optimize for maximum speed, but disable some optimizations which
       increase code size for a small speed benefit: /Ob1gysi-
/O2    optimize for maximum speed (same as /Ox)
/O3    enable /O2 plus more aggressive optimizations that may not improve
       performance for all programs
/Ox    enable maximum optimizations (DEFAULT):       /Ob1gyti
/Od    disable optimizations; useful for selective optimizations (i.e. /Od /Og)
/optimize:<val>
       0 - disable optimization (same as /Od)
       1 - optimize for maximum speed, but disable some optimizations
           which increase code size for a small speed benefit (same as /O1)
       2 - same as 1, or /O1
       3 - optimize for maximum speed (default) (same as /O2)
       4 - same as 3, or /O2
       5 - enable /optimize:3 plus more aggressive optimizations that
           may not improve performance for all programs (same as /O3)

/fast  enable /QxP /O3 /Qipo /Qprec-div-
/Ob<n> control inline expansion:
       n=0  disable inlining (same as /inline:none)
       n=1  no inlining (unless /Qip specified)
       n=2  inline any function, at the compiler's discretion
            (same as /Qip or /inline:all)
/Og    enable global optimizations
/Op[-] enable/disable better floating-point precision
/[no]fltconsistency
       specify that improved floating-point consistency should be used
/Os    enable speed optimizations, but disable some optimizations which
       increase code size for small speed benefit (overrides /Ot)
/Ot    enable speed optimizations (overrides /Os)
/Oy[-] enable/disable using EBP as general purpose register (no frame pointer)

/Qax<codes> generate code specialized for processors specified by <codes>
            while also generating generic IA-32 code.  <codes> includes
            one or more of the following characters:
    K  Intel Pentium III and compatible Intel processors
    W  Intel Pentium 4 and compatible Intel processors
    N  Intel Pentium 4 and compatible Intel processors.  Enables new
       optimizations in addition to Intel processor-specific optimizations
    P  Intel Pentium 4 processors with SSE3 extensions
    B  Intel Pentium M and compatible Intel processors
/Qx<codes>  generate specialized code to run exclusively on processors
            indicated by <codes> as described above.
/tune:<keyword>
     pn1 - optimize for Pentium(R) processor
     pn2 - optimize for Pentium(R) Pro, Pentium(R) II, and
           Pentium(R) III processors
     pn3 - same as pn2
     pn4 - optimize for Pentium(R) 4 processor (DEFAULT)
/architecture:<keyword>
     pn1 - optimize for Pentium(R) processor
     pn2 - optimize for Pentium(R) Pro, Pentium(R) II, and
           Pentium(R) III processors
     pn3 - same as pn2
     pn4 - optimize for Pentium(R) 4 processor (DEFAULT)
     SSE - same as QxK
     SSE2  same as QxW

Advanced Performance
--------------------
Enable and specify the scope of Interprocedural (IP) Optimizations:
/Qip     enable single-file IP optimizations (within files, same as /Ob2)
/Qipo[n] enable multi-file IP optimizations (between files)
/Qipo-c  generate a multi-file object file (ipo_out.obj)
/Qipo-S  generate a multi-file assembly file (ipo_out.asm)

Modify the behavior of IP:
/Qip-no-inlining    disable full and partial inlining (requires /Qip or /Qipo)
/Qip-no-pinlining   disable partial inlining (requires /Qip or /Qipo)
/Qipo-separate      create one object file for every source file
                    (overrides /Qipo[n])

Other Advanced Performance Options:
/reentrancy:<keyword>
                 specify whether the threaded, reentrant run-time support
                 should be used
                 Keywords:  none (same as /noreentrancy), threaded
/noreentrancy    don't use reentrant run-time support
/Qpad[-]         enable/disable(DEFAULT) changing variable and array memory
                 layout
/Qunroll[n]      set maximum number of times to unroll loops.  Omit n to use
                 default heuristics.  Use n=0 to disable loop unroller.
/unroll[:n]      set maximum number of times to unroll loops.  Omit n to use
                 default heuristics.  Use n=0 to disable loop unroller.

/Qprof-dir <d>   specify directory for profiling output files (*.dyn and *.dpi)
/Qprof-file <f>  specify file name for profiling summary file
/Qprof-gen       instrument program for profiling
/Qprof-use       enable use of profiling information during optimization
/Qprof-gen-sampling
                 prepare code for use with profrun sample gathering tool
/Qfnsplit[-]     enable/disable function splitting (enabled with /Qprof-use)

/Qpc32           set internal FPU precision to 24 bit significand
/Qpc64           set internal FPU precision to 53 bit significand (DEFAULT)
/Qpc80           set internal FPU precision to 64 bit significand
/QIfist[-]       enable/disable(DEFAULT) fast float-to-int conversions
/Qrcd            same as /QIfist
/Qrct            set internal FPU rounding control to truncate
/rounding-mode:chopped
                 set internal FPU rounding control to truncate
/Qprec           improve floating-point precision (speed impact less than /Op)
/Qfp-port        round fp results at assignments & casts (some speed impact)
/Qfpstkchk       enable fp stack checking after every function/procedure call
/Qprec-div[-]    improve precision of FP divides (some speed impact)
/Qprec-sqrt[-]
                 determine if certain square root optimizations are enabled

/Qscalar-rep[-]  enable(DEFAULT)/disable scalar replacement (requires /O3)
/Qvec-report[n]  control amount of vectorizer diagnostic information:
            n=0    no diagnostic information
            n=1    indicate vectorized loops (DEFAULT)
            n=2    indicate vectorized/non-vectorized loops
            n=3    indicate vectorized/non-vectorized loops and prohibiting
                   data dependence information
            n=4    indicate non-vectorized loops
            n=5    indicate non-vectorized loops and prohibiting data
                   dependence information

/Qprefetch[-]              enable(DEFAULT)/disable prefetch insertion
/Qsafe-cray-ptr            Cray pointers do not alias with other variables
/Qopt-report               generate an optimization report to stderr
/Qopt-report-file<file>    specify the filename for the generated report
/Qopt-report-level[level]  specify the level of report verbosity (min|med|max)
/Qopt-report-phase<name>   specify the phase that reports are generated against
/Qopt-report-routine<name> reports on routines containing the given name
/Qopt-report-help          display the optimization phases available for
                           reporting

/Qtcheck               generate instrumentation to detect multi-threading bugs
                       (requires Intel(R) Threading Tools; cannot be used with
                       compiler alone)
/Qopenmp               enable the compiler to generate multi-threaded code
                       based on the OpenMP directives
/Qopenmp-profile       link with instrumented OpenMP runtime library to
                       generate OpenMP profiling information for use with the
                       OpenMP component of the VTune(TM) Performance Analyzer
/Qopenmp-stubs         enables the user to compile OpenMP programs in
                       sequential mode.  The openmp directives are ignored and
                       a stub OpenMP library is linked (sequential)
/Qopenmp-report{0|1|2} control the OpenMP parallelizer diagnostic level
/fpe:{0|1|3}           specifies behavior on floating point exceptions

/Qparallel             enable the auto-parallelizer to generate multi-threaded
                       code for loops that can be safely executed in parallel
/Qpar-report{0|1|2|3}  control the auto-parallelizer diagnostic level
/Qpar-threshold[n]     set threshold for the auto-parallelization of loops
                       where n is an integer from 0 to 100
/Qansi-alias[-]        enable(DEFAULT)/disable use of ANSI aliasing rules in
                       optimizations; user asserts that the program adheres to
                       these rules
/Qcomplex-limited-range[-]
                       enable/disable(DEFAULT) the use of the basic
                       algebraic expansions of some complex arithmetic
                       operations.  This can allow for some performance
                       improvement in programs which use a lot of complex
                       arithmetic at the loss of some exponent range.
/Qftz[-]               enable/disable flush denormal results to zero
/[no]recursive         compile all procedures for possible recursive execution
/Qssp                  enable software-based speculative pre-computation

Output, Debug
-------------
/c, /nolink, /compile-only
                  compile to object (.obj) only, do not link
/S                compile to assembly (.asm) only, do not link
/FA[cs]           produce assembly file with optional code or source
                  annotations
/asmattr:<keyword>
                  specify the information to generate in the assembleable
                  listing file
                  keywords:  [no]source, [no]machine
/noasmattr        generate default information in the assembleable listing
/Fa[file]         name assembly file (or directory for multiple files; i.e.
                  /FaMYDIR\)
/asmfile[:filename]
                  specifies that an assembleable listing file should be
                  generated
/Fo[file]         name object file (or directory for multiple files; i.e.
                  /FoMYDIR2\)
/object:<filename>
                  specify the name of the object file
/Fe<file>         name executable file or directory
/exe:<filename>   specifies the name to be used for the built program (.EXE)
                  or dynamic-link (.DLL) library
/map:<filename>   specify that a link map file should be generated
/pdbfile[:filename]
                  specify that debug related information should be generated
                  to a program database file
/nopdbfile        do not generate debug related information to a program
                  database file
/RTCu             report use of variable that was not initialized
/Zi, /Z7          produce symbolic debug information in object file
/debug[:keyword]
                  enable debug information and control output of enhanced
                  debug information.
                  keywords:  all, full, minimal, none (same as /nodebug)
                             inline-debug-info, variable-locations,
                             semantic-stepping, extended
/nodebug          do not enable debug information
/debug-parameters[:keyword]
                  control output of debug information for PARAMETERS
                  keywords: all, used, none (same as /nodebug-parameters)
/nodebug-parameters
                  do not output debug information for PARAMETERS
/Quse-vcdebug     emit debug information compatible with Visual C++ debugger;
                  this disables emission of extended information used by
                  Fortran debuggers
/Qinline-debug-info
                  preserve the source position of inlined code instead of
                  assigning the call-site source position to inlined code
/Qtrapuv          trap uninitialized variables
/Qmap-opts        enable option mapping tool

Fortran Preprocessor
--------------------
/D<name>[{=|#}<text>], /define:symbol[=value]
           Defines the symbol specified for use with conditional
           compilation metacommands.  Optionally, the symbol can
           be defined to be an integer or string value
/nodefines specifies that any /D macros go to the preprocessor only, and not
           to the compiler too
/allow:[no]fpp-comments
           If a Fortran end-of-line comment is seen within a #define,
           treat it as part of the definition.  Default is nofpp-comments
/E         preprocess to stdout
/EP        preprocess to stdout, omitting #line directives
/EP /P     preprocess to file, omitting #line directives
/P, /preprocess-only
           preprocess to file
/Qfpp[n], /[no]fpp
           run Fortran preprocessor on source files prior to compilation
     n=0   disable running the preprocessor, equivalent to nofpp
     n=1,2,3 run preprocessor
/Qcpp[n]   same as /Qfpp[n]
/I<dir>, /[no]include:<dir>
           define the directories in which the compiler searches for INCLUDE
           and mod files
/module[:path]
           specify path where mod files should be placed and
           first location to look for mod files
/U<name>, /undefine:<name>
           remove predefined macro
/u         remove all predefined macros
/X         remove standard directories from include file search path

Language
--------
/4I{2|4|8}       set default KIND of integer and logical variables to 2, 4,
                 or 8
/integer-size:<size>
                 specifies the default size of integer and logical variables
                 size:  16, 32, 64
/4R{8|16}        set default size of REAL to 8 or 16 bytes
/real-size:<size>
                 specify the size of REAL and COMPLEX declarations, constants,
                 functions, and intrinsics
                size: 32, 64, 128
/Qautodouble     same as /4R8
/double-size:<size>
                 defines the size of DOUBLE PRECISION and DOUBLE COMPLEX
                 declarations, constants, functions, and intrinsics
                 size:  64, 128
/[no]fpconstant  extends the precision of single precision constants assigned
                 to double precision variables to double precision
/[no]intconstant use Compaq Fortran 77 semantics, rather than Compaq Fortran
                 90, to determine kind of integer constants

/Qd-lines, /[no]d-lines
                 compile debug statements (indicated by D in column 1)

/4{Y|N}f         enable/disable free format source
/[no]fixed, /FI  specifies source files are in fixed format
/[no]free , /FR  specifies source files are in free format
/4L{72|80|132}   specify 72, 80, or 132 column lines for fixed form sources
/Qextend-source  same as /4L132
/extend-source[:size]
                 determines the column that ends the statement field of each
                 source line in a fixed-format file.  /extend-source default
                 size is 132
                 size:  72 (same as /noextend-source), 80, 132
/noextend-source source lines are expected to be 72 characters
/ccdefault:<keyword>
                 specify default carriage control for units 6 and *
                 keywords:  default, fortran, list, or none
/stand:<keyword> specifies level of conformance with ANSI standard to check
                 for
                 keywords:  f90, f95, none
/Qpad-source, /[no]pad-source
                 make compiler acknowledge blanks at the end of a line
/Qdps[-], /[no]altparam
                 specify if alternate form of parameter constant declarations
                 (without parenthesis) is recognized

/1, /Qonetrip execute any DO loop at least once
/f66          allow extensions that enhance FORTRAN-66 compatibility
/[no]f77rtl   specify that the Fortran 77 specific run-time support should
              be used
/4{Y|N}a      enable/disable putting local variables on the run-time stack
/Qauto        same as /4Ya or /automatic
/Qauto-scalar make scalar local variables AUTOMATIC
/Qsave        save all variables (static allocation)
               (same as /noautomatic or /4Na, opposite of /Qauto)
/4{Y|N}d      enable/disable default IMPLICIT NONE
/Qcommon-args assume "by reference" subprogram arguments may alias one another
              (same as /assume:dummy_aliases)
/[Q]vms       enable VMS I/O statement extensions
/assume:<keyword>
              specify assumptions made by the optimizer and code generator
              keywords: none, [no]accuracy_sensitive, [no]byterecl,
                        [no]buffered_io, [no]dummy_aliases, [no]minus0,
                        [no]protect_constants, [no]source_include,
                        [no]underscore
/Qzero[-]     enable/disable(DEFAULT) implicitly initialize all data to zero

/Qdyncom"blk1,blk2,..."  make given common blocks dynamically-allocated

/4{Y|N}b enable/disable extensive runtime error checking
/CB      runtime checks for out-of-bounds array subscript/substring refs
         Same as /check:bounds
/CU      runtime check for uninitialized variables
         same as /check:uninit
/check:<keyword>
         check run-time conditions
         keywords: all (same as /4Yb), none (same as /nocheck, /4Nb),
                   [no]arg_temp_created, [no]bounds,[no]format, [no]overflow,
                   [no]output_conversion, [no]power, [no]uninit, [no]args
/nocheck perform no extra run-time checks

/nbs          treat backslash as a normal character, not an escape character
/us           append an underscore to external subroutine names
/Qlowercase   change routine names to all lowercase (same as /names:lowercase)
/Quppercase   change routine names to all uppercase (DEFAULT)
              (same as /names:uppercase)
/names:<keyword>
              specify how source code identifiers and external names are
              interpreted
              keywords:  as_is, lowercase, uppercase
/Gm           enable CVF and Powerstation calling convention compatibility
              (same as /iface:cvf)
/Gz           make STDCALL the default calling convention
              (same as /iface:stdcall)
/iface:<keyword>
              specify the calling conventions to be used
              keywords: [no]mixed_str_len_arg, cref, cvf,
                        default, stdcall, stdref
/align[-]     analyze and reorder memory layout for variables and arrays
/Zp[n]        specify alignment constraint for structures (n=1,2,4,8,16)
/align:<keyword>
              specify how data items are aligned
              keywords: [no]commons, dcommons, [no]records, rec1byte,
                        rec2byte, rec4byte, rec8byte, rec16byte, sequence
/Zs, /syntax-only
              perform syntax check only

Compiler Diagnostics
--------------------
/cm        same as /warn:nousage
/w         disable all warnings
/Wn        disable warnings (n=0) or show warnings
           (n=1, DEFAULT, same as /warn:general)
/w90, /w95 suppress messages about use of non-standard Fortran-90/95
/4{Y|N}s   enable/disable issuing of errors for non-standard Fortran-90
/warn:<keyword>
            specifies the level of warning messages issued.
           keywords: all, none (same as /nowarn),
                     [no]alignments, [no]argument_checking, [no]declarations,
                     [no]errors, [no]general, [no]ignore_loc,
                     [no]interfaces, [no]stderrors, [no]truncated_source,
                     [no]uncalled, [no]uninitialized, [no]unused, [no]usage
/nowarn    suppress all warning messages
/WB        turn a compile-time bounds check into a warning
/[no]gen-interfaces
           create interface blocks for all routines in the file.  Can then be
           checked with warn:interfaces
/[no]traceback
           specify whether the compiler generates PC correlation data used
           to display a symbolic traceback rather than a hexidecimal traceback
           at runtime failure
/error-limit:<size>
           specify the maximum number of error-level or fatal-level compiler
           errors allowed
/noerror-limit
            set no maximum number on error-level or fatal-level error messages

Miscellaneous
-------------
/?, /help     print this help message
/[no]logo     [do not] display compiler version information
/watch:<keyword>
              show driver tool commands but do not execute tools
/Tf<file>     compile file as Fortran source
/source:<filename>
              specify that the file indicated should be processed as a Fortran
              source file
/extfor:<ext> specify extension of file to be recognized as a Fortran file
/extfpp:<ext> specify extension of file to be recognized as a preprocessor
              file
/convert:<keyword>
              specify the format of unformatted files containing numeric data
              keywords: big_endian, cray, ibm, little_endian, native, vaxd,
                        vaxg
/fpscomp[:keyword]
              specify the level of compatibility to adhere to with Fortran
              PowerStation
              keywords: all, none (same as /nofpscomp),
                        [no]filesfromcmd, [no]general, [no]ioformat,
                        [no]ldio_spacing, [no]libs, [no]logicals
/nofpscomp    no specific level of compatibility with Fortran Powerstation
/libdir[:keyword]
              control the library names that should be emitted into the object
              file
              Keywords:  all, none (same as /nolibdir), [no]automatic, [no]user
/nolibdir     no library names should be emitted into the object file
/inline[:keyword]
              Specifies the level of inline function expansion.
              keywords: all (same as /Ob2 /Ot), size (same as /Ob2 /Os)
                        speed (same as /Ob2 /Ot), none or manual (same as /Ob0)
/Qoption,<str>,<opts>
              pass options <opts> to tool specified by <str>
/Qlocation,<str>,<dir>
              set <dir> as the location of tool specified by <str>
/bintext:<string>
              place the string specified into the object file and executable
/Qsox[-]      enable/disable(DEFAULT) saving of compiler options and version
              in the executable
/Qnobss-init  disable placement of zero-initialized vars in BSS (use DATA)

Linking/Linker
--------------
/link           specify that all following options are for the linker
/extlnk:<ext>   specify extension of file to be passed directly to linker
/F<n>           set the stack reserve amount specified to the linker
/Qvc6           use Visual C++ 6.0 linker and debug format (DEFAULT)
/Qvc7           enable Visual Studio .NET 2002 linker and debug format
/Qvc7.1         enable Visual Studio .NET 2003 linker and debug format
/LD[d], /dll    produce a DLL instead of an EXE ('LDd' = debug version)
/dbglibs        use the debug version of runtime libraries, when appropriate
/MD             use dynamically-loaded, multithread runtime
/MDs            use dynamically-loaded, single thread runtime
/MG, /winapp    use Windows API runtime libraries
/ML             use statically-loaded, single thread runtime (/ML is DEFAULT)
/MT             use statically-loaded, multithread runtime
/MW             use statically-loaded, multidocument QuickWin runtime
/MWs            use statically-loaded, single document QuickWin runtime
/libs:<keyword> determine the default libraries to be linked against
                dll    - use dynamically-loaded, multithread runtime
                qwin   - use statically-loaded, multidocument QuickWin runtime
                qwins  - use statically-loaded, single document QuickWin runtime
                static - use statically-loaded, single thread runtime (DEFAULT)
/static         synonym to /libs:static
/4{Y|N}portlib  enable/disable linking with portability library
/Zl             omit library names from object file (same as /nolibdir)
/[no]threads    specify whether or not multi-threaded libraries should be
                 linked against

Copyright (C) 1985-2005, Intel Corporation.  All rights reserved.

- Other brands and names are the property of their respective owners.
```
