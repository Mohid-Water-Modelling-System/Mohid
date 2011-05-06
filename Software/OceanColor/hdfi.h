/****************************************************************************
 * NCSA HDF                                                                 *
 * Software Development Group                                               *
 * National Center for Supercomputing Applications                          *
 * University of Illinois at Urbana-Champaign                               *
 * 605 E. Springfield, Champaign IL 61820                                   *
 *                                                                          *
 * For conditions of distribution and use, see the accompanying             *
 * hdf/COPYING file.                                                        *
 *                                                                          *
 ****************************************************************************/

/* $Id: hdfi.h,v 1.148 2001/09/20 21:37:31 ymuqun Exp $ */

#ifndef HDFI_H
#define HDFI_H

#ifdef GOT_MACHINE
#undef GOT_MACHINE
#endif

/*--------------------------------------------------------------------------*/
/*                              MT/NT constants                             */
/*  Four MT nibbles represent double, float, int, uchar (from most          */
/*      significant to least significant).                                  */
/*  Each "column" in the "table" below is essentially independant of the    */
/*      other "columns", for example the CONVEXNATIVE entry means that the  */
/*      floating point formats are in Convex native format but the integers */
/*      are big-endian and standard sizes                                   */
/*  If you add another value to this "table", you need to add another       */
/*      DFNTF_xxx entry in hntdefs.h                                        */
/*  The values for each nibble are:                                         */
/*      1 - Big Endian                                                      */
/*          (i.e. Big-Endian, 32-bit architecture w/IEEE Floats)            */
/*      2 - VAX                                                             */
/*          (i.e. Middle-Endian, 32-bit architecture w/VAX Floats)          */
/*      3 - Cray                                                            */
/*          (i.e. Big-Endian, all 64-bit architecture w/Cray Floats)        */
/*      4 - Little Endian                                                   */
/*          (i.e. Little-Endian, 32-bit architecture w/IEEE Floats)         */
/*      5 - Convex                                                          */
/*          (i.e. Big-Endian, 32-bit architecture w/Convex Native Floats)   */
/*      6 - Fujitsu VP                                                      */
/*          (i.e. Big-Endian, 32-bit architecture w/Fujitsu Native Floats)  */
/*      7 - Cray MPP                                                        */
/*          (i.e. Big-Endian, 32-bit architecture w/IEEE Floats, but no 16-bit type)            */
/*      8 - Cray IEEE                                                       */
/*          (i.e. Big-Endian, all 64-bit architecture w/IEEE Floats)        */
/*--------------------------------------------------------------------------*/
#define     DFMT_SUN            0x1111
#define     DFMT_ALLIANT        0x1111
#define     DFMT_IRIX           0x1111
#define     DFMT_APOLLO         0x1111
#define     DFMT_IBM6000        0x1111
#define     DFMT_HP9000         0x1111
#define     DFMT_CONVEXNATIVE   0x5511
#define     DFMT_CONVEX         0x1111
#define     DFMT_UNICOS         0x3331
#define     DFMT_UNICOSIEEE     0x1831
#define     DFMT_CTSS           0x3331
#define     DFMT_VAX            0x2221
#define     DFMT_MIPSEL         0x4441
#define     DFMT_PC             0x4441
#define     DFMT_MAC            0x1111
#define     DFMT_SUN386         0x4441
#define     DFMT_NEXT           0x1111
#define     DFMT_MOTOROLA       0x1111
#define     DFMT_ALPHA          0x4441
#define     DFMT_VP             0x6611
#define     DFMT_I860           0x4441
#define     DFMT_CRAYMPP        0x1171
#define     DFMT_IA64           0x4441

/* I/O library constants */
#define UNIXUNBUFIO 1
#define UNIXBUFIO   2
#define MACIO       3
#define PCIO        4    /* 16-bit MS-DOS File I/O (deprecated) */
#define WINIO       5    /* 16-bit Windows File I/O (deprecated) */
#define PAGEBUFIO   6    /* page buffering - fmpool */
#define WINNTIO     7    /* 32-bit Windows File I/O (deprecated, WinNT now uses UNIXBUFIO) */ 

/* IBM RS6000 AIX hack */
#if defined(IBM6000) || defined(_AIX)
#define _POSIX_SOURCE
#endif 

/* Standard header files needed all the time */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

/* PABLO support files */
#ifdef HAVE_PABLO
#define HDFIOTRACE
#include "HDFIOTrace.h"
#endif  /* HAVE_PABLO */

/**
 * Provide the macros to adapt the HDF public functions to
 * dll entry points.
 * In addition it provides error lines if the configuration is incorrect.
 **/

#ifdef WIN32
/**
 * Under WIN32 we have single threaded static libraries, or
 * mutli-threaded DLLs using the multithreaded runtime DLLs.
 **/
#	if defined(_MT) &&	defined(_DLL) &&!defined(_HDFDLL_)
/*		If the user really ment to use _HDFDLL_, but he forgot, just define it. */
#		define _HDFDLL_
#	endif

#	if defined(_MT) && !defined(_HDFDLL_)
#		error To use the HDF libraries from a multithreaded project, you must use the HDF DLLs
#		error Use the Mutlithreaded DLL runtime libraries (prefered), or define the macro "_HDFDLL_"
#	endif

#	if !defined(_MT) && defined(_HDFDLL_)
#		error To use the HDF libraries from a single-threaded project, you must use static HDF libraries
#		error Undefine the macro "_HDFDLL_"
#	endif

#	if defined(_HDFDLL_)
#		pragma warning( disable: 4273 )	/* Disable the stupid dll linkage warnings */

#		if !defined(_HDFLIB_)
#			define HDFPUBLIC __declspec(dllimport)
#		else
#			define HDFPUBLIC __declspec(dllexport)
#		endif

#		if !defined(_MFHDFLIB_) && !defined(_HDFLIB_)
#			define HDFLIBAPI __declspec(dllimport) extern
#		else
#			define HDFLIBAPI __declspec(dllexport) extern
#		endif 
#	else
#		define HDFPUBLIC
#		define HDFLIBAPI extern
#	endif
#else	/* !defined( WIN32 ) */
#	define HDFPUBLIC
#	define HDFLIBAPI extern
#endif


/*-------------------------------------------------------------------------
 * Define options for each platform
 *-------------------------------------------------------------------------*/

/*
 * Meaning of each defined macros (not completed yet)
 *
 * BIG_LONGS--Define when long is not "equal" to int32.  True in cases
 *      where (int32 *) is not compatible with (long *).  Should
 *      be renamed as LONGNEINT32.
 */

#if (defined(SUN) || defined(sun) || defined(__sun__) || defined(__SUNPRO_C)) & !defined(__i386)
#ifdef __STDC__
#define ANSISUN
#else /* __STDC__ */
#define KNRSUN
#endif /* __STDC__ */
#endif /* SUN || sun */

#if defined(ANSISUN)

#if !defined(SUN)
#define SUN
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#include <unistd.h>                 /* for some file I/O stuff */
#include <sys/time.h>
#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT             DFMT_SUN
typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
#ifdef _LP64 /* 64-bit environment */
typedef int               int32;
typedef unsigned int      uint32;
#else /* 32-bit environment */
typedef long int          int32;
typedef unsigned long int uint32;
#endif
typedef int               intn;
typedef unsigned int      uintn;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef float             float32;
typedef double            float64;
#ifdef _LP64 /* 64-bit environment */
typedef long              hdf_pint_t;   /* an integer the same size as a pointer */
#else /* 32-bit environment */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#endif
#define FNAME_POST_UNDERSCORE
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI
#define HAVE_STDC
#define INCLUDES_ARE_ANSI

#endif /* ANSISUN */

#if defined(KNRSUN)

#if !defined(SUN)
#define SUN
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#   define BSD
#define DUMBCC 	/* because it is.  for later use in macros */
#ifndef __GNUC__
#include <memory.h>
#endif /* __GNUC__ */
#include <unistd.h>
#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT             DFMT_SUN
typedef void              VOID;
typedef char              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
typedef long int          int32;
typedef unsigned long int uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef float             float32;
typedef double            float64;
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define FNAME_POST_UNDERSCORE
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#ifdef __GNUC__
#define HAVE_STDC
#define INCLUDES_ARE_ANSI
#endif

#endif /* SUN */


#if defined(IBM6000) || defined(_AIX)

#ifndef IBM6000
#define IBM6000
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#   define BSD

#ifndef __GNUC__
#include <memory.h>
#endif /* __GNUC__ */
#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT             DFMT_IBM6000
typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
#ifndef _ALL_SOURCE       
typedef char              int8;
typedef short int         int16; 
typedef long int          int32; 
#endif  
typedef char              char8;
typedef unsigned char     uchar8;
typedef unsigned char     uint8;
typedef unsigned short int uint16;
typedef unsigned long int uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef float             float32;
typedef double            float64;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#define HAVE_STDC
#define INCLUDES_ARE_ANSI

#endif /* IBM6000 */

#if defined(HP9000) || (!defined(__convexc__) && (defined(hpux) || defined(__hpux)))

#ifndef HP9000
#define HP9000
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#ifndef HAVE_UNISTD_H
#define HAVE_UNISTD_H  /* unistd.h - close, fork,..etc */
#endif

#   define BSD
#ifndef __GNUC__
#include <memory.h>
#endif /* __GNUC__ */
#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT             DFMT_HP9000
typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
typedef long int          int32;
typedef unsigned long int uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef float             float32;
typedef double            float64;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#endif /* HP9000 */


#if defined(IRIX) || defined(IRIS4) || defined(sgi) || defined(__sgi__) || defined(__sgi)

#ifndef IRIX
#define IRIX
#endif

#if (_MIPS_SZLONG == 64)
/* IRIX 64 bits objects.  It is nearly the same as the conventional
 * 32 bits objects.  Let them share IRIX definitions for now.
 */
#define IRIX64
#endif


#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1

/*
 * BSD was originally defined with no value.  But some newer SGI system
 * header files (e.g., resolv.h) assume it has a value and evaluate it
 * in expressions, thus causing compiling errors.  This has been reported
 * to SGI as bug #781568.  SGI could not provide a list of the semantics
 * of BSD values and suggested a work around of setting BSD to 1.
 */
#   define BSD 1
#ifndef __GNUC__
#include <memory.h>
#endif /* __GNUC__ */
#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT              DFMT_IRIX
typedef void               VOID;
typedef void               *VOIDP;
typedef char               *_fcd;
typedef signed char        char8;
typedef unsigned char      uchar8;
typedef signed char        int8;
typedef unsigned char      uint8;
typedef short int          int16;
typedef unsigned short int uint16;
typedef int                int32;
typedef unsigned int       uint32;
typedef int                intn;
typedef unsigned int       uintn;
typedef float              float32;
typedef double             float64;
typedef int                intf;     /* size of INTEGERs in Fortran compiler */
typedef long               hdf_pint_t;   /* an integer the same size as a pointer */
#define FNAME_POST_UNDERSCORE
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif
#ifdef IRIX64
#define BIG_LONGS
#endif


/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#define HAVE_STDC
#define INCLUDES_ARE_ANSI

#endif /* IRIX */

#if (defined(UNICOS) || defined(_UNICOS)) && !defined(_CRAYMPP)

#ifndef UNICOS
#define UNICOS
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1

#include <memory.h>
#include <fortran.h>
#ifndef O_RDONLY
#include <fcntl.h>              /* for unbuffered i/o stuff */
#define L_INCR  1
#include <sys/stat.h>
#endif /*O_RDONLY*/

#ifdef _CRAYIEEE
#define DF_MT   DFMT_UNICOSIEEE
#else
#define DF_MT   DFMT_UNICOS
#endif
typedef void            VOID;
typedef void            *VOIDP;
#ifdef OLD_WAY /* May need to be included on other machines than the C-90 */
typedef char            *_fcd;
#endif /* OLD_WAY */
typedef signed char     char8;
typedef unsigned char   uchar8;
typedef signed char     int8;
typedef unsigned char   uint8;
typedef int             int16;
typedef unsigned int    uint16;
typedef int             int32;
typedef unsigned int    uint32;
typedef int             intn;
typedef unsigned int    uintn;
typedef float           float32;
typedef double          float64;
typedef int             intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */

#define DF_CAPFNAMES            /* fortran names are in all caps */
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI
#define RIGHT_SHIFT_IS_UNSIGNED
#define CHAR_IS_UNSIGNED

#endif /* UNICOS */

#if defined(_CRAYMPP)

#ifndef CRAYMPP
#define CRAYMPP
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1

#include <string.h>
#include <limits.h>
#include <memory.h>
#include <fortran.h>
#ifndef O_RDONLY
#include <fcntl.h>              /* for unbuffered i/o stuff */
#define L_INCR  1
#include <sys/stat.h>
#endif /*O_RDONLY*/

#define DF_MT   DFMT_CRAYMPP
typedef void            VOID;
typedef void            *VOIDP;
#ifdef OLD_WAY /* May need to be included on other machines than the C-90 */
typedef char            *_fcd;
#endif /* OLD_WAY */
typedef signed char     char8;
typedef unsigned char   uchar8;
typedef signed char     int8;
typedef unsigned char   uint8;
typedef short           int16;
typedef unsigned short  uint16;
typedef short           int32;
typedef unsigned short  uint32;
typedef int             intn;
typedef unsigned int    uintn;
typedef float           float32;
typedef double          float64;
typedef int             intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */

#define _HUGE              /* This should only be defined to a value on the PC */
#define DF_CAPFNAMES            /* fortran names are in all caps */
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI
#define RIGHT_SHIFT_IS_UNSIGNED
#define CHAR_IS_UNSIGNED

#endif /* CRAYMPP */

#if defined(VMS) || defined(vms)

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1
#include <file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT              DFMT_VAX
typedef void               VOID;
typedef void               *VOIDP;
typedef char               *_fcd;
typedef char               char8;
typedef unsigned char      uchar8;
typedef char               int8;
typedef unsigned char      uint8;
typedef short int          int16;
typedef unsigned short int uint16;
#ifdef __alpha
typedef int                int32;
typedef unsigned int       uint32;
typedef long              hdf_pint_t;   /* an integer the same size as a pointer */
#else
typedef long int           int32;
typedef unsigned long int  uint32;
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#endif
typedef int                intn;
typedef unsigned int       uintn;
typedef float              float32;
typedef double             float64;
typedef int                intf;     /* size of INTEGERs in Fortran compiler */
#define _fcdtocp(desc)  ((char *) *((char **) &desc[4]))

/* 
  Redef a couple of C routine names to avoid conflicts
  since the VMS link command is case-insensitive
*/
#define FILELIB UNIXBUFIO
#define DF_CAPFNAMES            /* fortran names are in all caps */
#include "dfivms.h"


/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#endif /* VMS */

#if defined(CONVEX) || defined(CONVEXNATIVE) || defined(__convexc__)

#ifndef CONVEX
#define CONVEX
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#include <sys/types.h>
#include <sys/stat.h>
/* For Convex machines with native format floats */
#ifdef CONVEXNATIVE
#define DF_MT             DFMT_CONVEXNATIVE
#else
#define DF_MT             DFMT_CONVEX
#endif
typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
typedef long int          int32;
typedef unsigned long int uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef float             float32;
typedef double            float64;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI
#define RIGHT_SHIFT_IS_UNSIGNED
#define INCLUDES_ARE_ANSI
#define HAVE_STDC

#endif /* CONVEX */


#if defined(MIPSEL) || ((defined(mips) || defined(__mips)) && (defined(ultrix) || defined(__ultrix)))

#ifndef MIPSEL
#define MIPSEL
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1

#ifndef __GNUC__
#define DUMBCC 	/* because it is.  for later use in macros */
#endif /* __GNUC__ */

#include <sys/types.h>
#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT   DFMT_MIPSEL
typedef void            VOID;
typedef void            *VOIDP;
typedef char            *_fcd;
typedef char            char8;
typedef unsigned char   uchar8;
typedef char            int8;
typedef unsigned char   uint8;
typedef short           int16;
typedef unsigned short  uint16;
typedef int             int32;
typedef unsigned int    uint32;
typedef int             intn;
typedef unsigned int    uintn;
typedef float           float32;
typedef double          float64;
typedef int             intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define _fcdtocp(desc) (desc)
#define FNAME_POST_UNDERSCORE
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#endif /* MIPSEL */

#if defined(MAC) || defined(macintosh) || defined (SYMANTEC_C)

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1

#include <memory.h>             /* malloc stuff for MPW */
#include <fcntl.h>              /* unbuffered I/O stuff for MPW */
#ifdef __MWERKS__				/* Metrowerks */
#include <sioux.h>
#include <console.h>
#endif
#ifdef SYMANTEC_C				/* for SYMANTEC C */
#include <unix.h>
#define isascii(c)  (isprint(c) || iscntrl(c))
#else  /* MPW, possibly others */
#include <Files.h>              /* for unbuffered I/O stuff */
#endif /* SYMANTEC_C*/
#ifndef ABSOFT
#define DF_CAPFNAMES            /* fortran names are in all caps */
#endif /* ABSOFT */
#define DF_DYNAMIC              /* use dynamic allocation */
#define DF_MT   DFMT_MAC

typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
typedef long int          int32;
typedef unsigned long int uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef float             float32;
typedef double            float64;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define _fcdtocp(desc) (desc)
void exit(int status);

#define FILELIB MACIO

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#endif /*MAC*/

/* Metrowerks Mac compiler defines some PC stuff so need to exclude this on the Mac */
#if !(defined(macintosh) || defined(MAC))

#if defined _M_ALPHA || defined _M_IX86 || defined INTEL86 || defined M_I86 || defined M_I386 || defined DOS386 || defined __i386 || defined UNIX386 || defined i386
#ifndef INTEL86
#define INTEL86
#endif /* INTEL86 */

#if !defined UNIX386 && (defined unix || defined __unix)
#define UNIX386
#endif /* UNIX386 */

#if !defined DOS386 && defined M_I386
#define DOS386
#endif /* M_I386 && !DOS386 */

#if defined _WINDOWS || defined WIN32
#define WIN386
#endif  /* _WINDOWS | WIN32 */

#if defined WIN386 || defined DOS386 || defined UNIX386
#define INTEL386
#endif /* WIN386 | DOS386 | UNIX386 */

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1

#if defined _WINDOWS || defined WIN32
#pragma comment( lib, "oldnames" )
#endif

#include <fcntl.h>
#ifdef UNIX386
#include <sys/types.h>      /* for unbuffered file I/O */
#include <sys/stat.h>
#include <unistd.h>
#else /* !UNIX386 */
#include <sys\types.h>      /* for unbuffered file I/O */
#include <sys\stat.h>
#include <io.h>
#include <conio.h>          /* for debugging getch() calls */
#include <malloc.h>
#endif /* UNIX386 */
#include <ctype.h>          /* for character macros */
#ifdef __WATCOMC__
#include <stddef.h>         /* for the 'fortran' pragma */
#endif

#if defined(_MSC_VER) && !defined(_MFHDFLIB_) && !defined(_HDFLIB_)	/* Auto-link when possible */
#	define HDF_LIB_VER	"425"
#	if !defined(_DEBUG)
#		if !defined(_HDFDLL_)
#			define HDF_LIB_NAME	"HD" HDF_LIB_VER ".lib"
#			pragma message( "Automatic linking with the static single-threaded HDF library - " HDF_LIB_NAME )
#		else
#			define HDF_LIB_NAME	"HD" HDF_LIB_VER "m.lib"
#			pragma message( "Automatic linking with the multithreaded HDF DLL - " HDF_LIB_NAME )
#		endif
#	else
#		if !defined(_HDFDLL_)
#			define HDF_LIB_NAME	"HD" HDF_LIB_VER "d.lib"
#			pragma message( "Automatic linking with the debug static single-threaded HDF library - " HDF_LIB_NAME  )
#		else
#			define HDF_LIB_NAME	"HD" HDF_LIB_VER "md.lib"
#			pragma message( "Automatic linking with the debug multithreaded HDF DLL - " HDF_LIB_NAME  )
#		endif
#	endif
#	pragma comment(lib, HDF_LIB_NAME )
#endif /* defined(_MSC_VER) && !defined(_MFHDFLIB_) && !defined(_HDFLIB_) */

#if defined WIN386
#ifndef GMEM_MOVEABLE       /* check if windows header is already included */
#include <windows.h>        /* include the windows headers */
#include <winnt.h>
#define HAVE_BOOLEAN
#endif /* GMEM_MOVEABLE */
#endif /* WIN386 */

#define DF_MT             DFMT_PC

#ifndef VOID    /* The stupid windows.h header file uses a #define instead of a typedef */
typedef void              VOID;
#endif  /* end VOID */
typedef void *            VOIDP;
typedef char *            _fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
typedef long int          int32;
typedef unsigned long int uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef float             float32;
typedef double            float64;
typedef long              intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */

#if defined _M_ALPHA
#define FNAME_PRE_UNDERSCORE
#endif

#if defined UNIX386
#define FNAME_POST_UNDERSCORE
#elif defined INTEL386
#define DF_CAPFNAMES
#endif
#define _fcdtocp(desc) (desc)

#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI
#define HAVE_STDC
#define INCLUDES_ARE_ANSI

#endif /* INTEL86 */
#endif /* !(defined(macintosh) || defined(MAC)) */

#if defined(NEXT) || defined(NeXT)

#ifndef NEXT
#define NEXT
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#define isascii(c)  (isprint(c) || iscntrl(c))
#ifndef __GNUC__
#include <memory.h>
#endif /* __GNUC__ */
#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT             DFMT_NEXT
typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
typedef long int          int32;
typedef unsigned long int uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef float             float32;
typedef double            float64;
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define FNAME_POST_UNDERSCORE
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#define HAVE_STDC
#define INCLUDES_ARE_ANSI

#endif /* NEXT */

#if defined(MOTOROLA) || defined(m88k)

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#ifndef __GNUC__
#include <memory.h>
#endif /* __GNUC__ */
#include <unistd.h>
#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#ifndef O_RDONLY
#include <fcntl.h>              /* for unbuffered i/o stuff */
#endif /*O_RDONLY*/
#define DF_MT             DFMT_MOTOROLA
typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
typedef long int          int32;
typedef unsigned long int uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef float             float32;
typedef double            float64;
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define FNAME_POST_UNDERSCORE
#define _fcdtocp(desc) (desc)
#define FILELIB UNIXBUFIO

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#endif /* MOTOROLA */

#if defined DEC_ALPHA || (defined __alpha && defined __unix__)

#ifndef DEC_ALPHA
#define DEC_ALPHA
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT             DFMT_ALPHA
typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
#ifndef __rpc_types_h
typedef int               int32;
typedef unsigned int      uint32;
#endif /* __rpc_types_h */
typedef int               intn;
typedef unsigned int      uintn;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef float             float32;
typedef double            float64;
typedef long              hdf_pint_t;   /* an integer the same size as a pointer */
#define FNAME_POST_UNDERSCORE
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#ifdef __GNUC__
#define HAVE_STDC
#define INCLUDES_ARE_ANSI
#endif

#endif /* DEC_ALPHA */

#if defined VP | defined __uxpm__

#ifndef VP
#define VP
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1

#include <memory.h>
#include <sys/types.h>
#include <sys/stat.h>
#define DF_MT              DFMT_VP
typedef void                VOID;
typedef void               *VOIDP;
typedef char               *_fcd;
typedef char               char8;
typedef unsigned char      uchar8;
typedef char               int8;
typedef unsigned char      uint8;
typedef short int          int16;
typedef unsigned short int uint16;
typedef long int           int32;
typedef unsigned long int  uint32;
typedef int                intn;
typedef unsigned int       uintn;
typedef int                intf;     /* size of INTEGERs in Fortran compiler */
typedef float              float32;
typedef double             float64;
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define FNAME_POST_UNDERSCORE
#define _fcdtocp(desc) (desc)
#define FILELIB UNIXBUFIO

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#endif /* VP */

#if defined I860 | defined i860

#ifndef I860
#define I860
#endif

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE 1

#include <sys/types.h>
#include <sys/file.h>           /* for unbuffered i/o stuff */
#include <sys/stat.h>
#include <unistd.h>             /* mis-using def. for SEEK_SET, but oh well */
#define DF_MT   DFMT_I860
typedef void            VOID;
typedef void            *VOIDP;
typedef char            *_fcd;
typedef char            char8;
typedef unsigned char   uchar8;
typedef char            int8;
typedef unsigned char   uint8;
typedef short           int16;
typedef unsigned short  uint16;
typedef int             int32;
typedef unsigned int    uint32;
typedef int             intn;
typedef unsigned int    uintn;
typedef float           float32;
typedef double          float64;
typedef int             intf;     /* size of INTEGERs in Fortran compiler */
typedef int               hdf_pint_t;   /* an integer the same size as a pointer */
#define _fcdtocp(desc) (desc)
#define FNAME_POST_UNDERSCORE
#define FILELIB UNIXBUFIO

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#endif /* I860 */


/* IA64 running Linux */
#if defined IA64

#ifdef GOT_MACHINE
If you get an error on this line more than one machine type has been defined.
Please check your Makefile.
#endif
#define GOT_MACHINE

#include <sys/file.h>               /* for unbuffered i/o stuff */
#include <sys/stat.h>
#define DF_MT             DFMT_IA64
typedef void              VOID;
typedef void              *VOIDP;
typedef char              *_fcd;
typedef char              char8;
typedef unsigned char     uchar8;
typedef char              int8;
typedef unsigned char     uint8;
typedef short int         int16;
typedef unsigned short int uint16;
typedef int               int32;
typedef unsigned int      uint32;
typedef int               intn;
typedef unsigned int      uintn;
typedef int               intf;     /* size of INTEGERs in Fortran compiler */
typedef float             float32;
typedef double            float64;
typedef long              hdf_pint_t;   /* an integer the same size as a pointer */
#define FNAME_POST_UNDERSCORE
#define _fcdtocp(desc) (desc)
#ifdef  HAVE_FMPOOL
#define FILELIB PAGEBUFIO  /* enable page buffering */
#else
#define FILELIB UNIXBUFIO
#endif

/* JPEG #define's - Look in the JPEG docs before changing - (Q) */

/* Determine the memory manager we are going to use. Valid values are: */
/*  MEM_DOS, MEM_ANSI, MEM_NAME, MEM_NOBS.  See the JPEG docs for details on */
/*  what each does */
#define JMEMSYS         MEM_ANSI

#ifdef __GNUC__
#define HAVE_STDC
#define INCLUDES_ARE_ANSI
#endif

#endif /* IA64 */

#ifndef GOT_MACHINE
No machine type has been defined.  Your Makefile needs to have someing like
-DSUN or -DUNICOS in order for the HDF internal structures to be defined
correctly.
#endif

/*-----------------------------------------------------*/
/*              encode and decode macros               */
/*-----------------------------------------------------*/

#   define INT16ENCODE(p, i) \
{ *(p) = (uint8)(((uintn)(i) >> 8) & 0xff); (p)++; \
        *(p) = (uint8)((uintn)(i) & 0xff); (p)++; }

#   define UINT16ENCODE(p, i) \
{ *(p) = (uint8)(((uintn)(i) >> 8) & 0xff); (p)++; *(p) = (uint8)((i) & 0xff); (p)++; }

#   define INT32ENCODE(p, i) \
{ *(p) = (uint8)(((uint32)(i) >> 24) & 0xff); (p)++; \
        *(p) = (uint8)(((uint32)(i) >> 16) & 0xff); (p)++; \
        *(p) = (uint8)(((uint32)(i) >> 8) & 0xff); (p)++; \
        *(p) = (uint8)((uint32)(i) & 0xff); (p)++; }

#   define UINT32ENCODE(p, i) \
{ *(p) = (uint8)(((i) >> 24) & 0xff); (p)++; \
        *(p) = (uint8)(((i) >> 16) & 0xff); (p)++; \
        *(p) = (uint8)(((i) >> 8) & 0xff); (p)++; \
        *(p) = (uint8)((i) & 0xff); (p)++; }

#   define NBYTEENCODE(d, s, n) \
{   HDmemcpy(d,s,n); p+=n }

/* DECODE converts big endian bytes pointed by p to integer values and store
 * it in i.  For signed values, need to do sign-extension when converting
 * the 1st byte which carries the sign bit.
 * The macros does not require i be of a certain byte sizes.  It just requires
 * i be big enough to hold the intended value range.  E.g. INT16DECODE works
 * correctly even if i is actually a 64bit int like in a Cray.
 */

#   define INT16DECODE(p, i) \
{ (i) = ((*(p) & 0x80) ? ~0xffff : 0x00) | ((int16)(*(p) & 0xff) << 8); (p)++; \
        (i) |= (int16)((*(p) & 0xff)); (p)++; }

#   define UINT16DECODE(p, i) \
{ (i) = (uint16)((*(p) & 0xff) << 8); (p)++; \
        (i) |= (uint16)(*(p) & 0xff); (p)++; }

#   define INT32DECODE(p, i) \
{ (i) = ((*(p) & 0x80) ? ~0xffffffff : 0x00) | ((int32)(*(p) & 0xff) << 24); (p)++; \
        (i) |= ((int32)(*(p) & 0xff) << 16); (p)++; \
        (i) |= ((int32)(*(p) & 0xff) << 8); (p)++; \
        (i) |= (*(p) & 0xff); (p)++; }

#   define UINT32DECODE(p, i) \
{ (i) = ((uint32)(*(p) & 0xff) << 24); (p)++; \
        (i) |= ((uint32)(*(p) & 0xff) << 16); (p)++; \
        (i) |= ((uint32)(*(p) & 0xff) << 8); (p)++; \
        (i) |= (uint32)(*(p) & 0xff); (p)++; }

/* Note! the NBYTEDECODE macro is backwards from the memcpy() routine, */
/*      in the spirit of the other DECODE macros */
#   define NBYTEDECODE(s, d, n) \
{   HDmemcpy(d,s,n); p+=n }

/*----------------------------------------------------------------
** MACRO FCALLKEYW for any special fortran-C stub keyword
**
** MacIntosh MPW LS-fortran needs pascal since it can interface
**  best with pascal functions.
** Microsoft C and Fortran need __fortran for Fortran callable C
**  routines.
**
** MACRO FRETVAL for any return value from a fortran-C stub function
**  Replaces the older FCALLKEYW macro.
**---------------------------------------------------------------*/
#ifdef FRETVAL
#undef FRETVAL
#endif

#if defined(MAC)                /* with LS FORTRAN */
#ifndef ABSOFT
#   define FCALLKEYW    pascal
#   define FRETVAL(x)   pascal x
#endif /* ABSOFT */
#endif

#ifndef FRETVAL /* !MAC */
#   define FCALLKEYW    /*NONE*/
#   define FRETVAL(x)   x
#endif


/*----------------------------------------------------------------
** MACRO FNAME for any fortran callable routine name.
**
**  This macro prepends, appends, or does not modify a name
**  passed as a macro parameter to it based on the FNAME_PRE_UNDERSCORE,
**  FNAME_POST_UNDERSCORE macros set for a specific system.
**
**---------------------------------------------------------------*/
#if defined(FNAME_PRE_UNDERSCORE) && defined(FNAME_POST_UNDERSCORE)
#   define FNAME(x)     _##x##_
#endif
#if defined(FNAME_PRE_UNDERSCORE) && !defined(FNAME_POST_UNDERSCORE)
#   define FNAME(x)     _##x
#endif
#if !defined(FNAME_PRE_UNDERSCORE) && defined(FNAME_POST_UNDERSCORE)
#   define FNAME(x)     x##_
#endif
#if !defined(FNAME_PRE_UNDERSCORE) && !defined(FNAME_POST_UNDERSCORE)
#   define FNAME(x)     x
#endif

/**************************************************************************
*  Generally useful macro definitions
**************************************************************************/
#ifndef MIN
#define MIN(a,b)    (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b)    (((a)>(b)) ? (a) : (b))
#endif

/**************************************************************************
*  Debugging Allocation functions
**************************************************************************/
#ifdef MALDEBUG
#include "maldebug.h"
#endif

/**************************************************************************
*  Macros to work around ANSI C portability problems.
**************************************************************************/
#ifdef DUMBCC
#define CONSTR(v,s) char *v=s
#else
#define CONSTR(v,s) static const char v[]=s
#endif

/* Old-style memory allocation function aliases -QAK */
#define HDgetspace HDmalloc
#define HDclearspace HDcalloc
#define HDregetspace HDrealloc
#define HDfreespace HDfree

/**************************************************************************
*  Allocation functions defined differently 
**************************************************************************/
#if !defined(MALLOC_CHECK) && !defined(_HDFDLL_)
#  define HDmalloc(s)      (malloc((size_t)s))
#  define HDcalloc(a,b)    (calloc((size_t)a,(size_t)b))
#  define HDfree(p)        (free((void*)p))
#  define HDrealloc(p,s)   (realloc((void*)p,(size_t)s))
#endif /* !defined MALLOC_CHECK */
/* Macro to free space and clear pointer to NULL */
#define HDfreenclear(p) { if((p)!=NULL) HDfree(p); p=NULL; }

/**************************************************************************
*  String functions defined differently 
**************************************************************************/

#  define HDstrcat(s1,s2)   (strcat((s1),(s2)))
#  define HDstrcmp(s,t)     (strcmp((s),(t)))
#  define HDstrcpy(s,d)     (strcpy((s),(d)))
#  define HDstrlen(s)       (strlen((const char *)(s)))
#  define HDstrncmp(s1,s2,n)    (strncmp((s1),(s2),(n)))
#  define HDstrncpy(s1,s2,n)    (strncpy((s1),(s2),(n)))
#  define HDstrchr(s,c)         (strchr((s),(c)))
#  define HDstrrchr(s,c)        (strrchr((s),(c)))
#  define HDstrtol(s,e,b)       (strtol((s),(e),(b)))
/* non-standard function, not defined on the following machines - */
#if !(defined VMS || defined macintosh || defined MAC || defined SYMANTEC_C || defined MIPSEL || defined NEXT || defined CONVEX || defined IBM6000 || defined ANSISUN || defined IRIX || defined _HDFDLL_ )
#  define HDstrdup(s)      ((char *)strdup((const char *)(s)))
#endif /* !(VMS | etc..) */


/**************************************************************************
*  Memory functions defined differently
**************************************************************************/

# define HDmemcpy(dst,src,n)   (memcpy((void *)(dst),(const void *)(src),(size_t)(n)))
# define HDmemset(dst,c,n)     (memset((void *)(dst),(intn)(c),(size_t)(n)))
# define HDmemcmp(dst,src,n)   (memcmp((const void *)(dst),(const void *)(src),(size_t)(n)))


/**************************************************************************
*  Misc. functions
**************************************************************************/
#if defined (MAC) || defined (macintosh) || defined (SYMANTEC_C)
#define HDstat(path, result)	(mstat(path))
#else /* !macintosh */
#define HDstat(path, result)	(stat(path, result))
#endif /* !macintosh */
#define HDgetenv(s1)            (getenv(s1))
#define HDputenv(s1)            (putenv(s1))
#define HDltoa(v)               (ltoa(v))
#if defined (SUN) && defined(__GNUC__)
#define HDatexit(f)             (0) /* we punt on the SUN using gcc */
#else /* !SUN & GCC */
#define HDatexit(f)             (atexit(f))
#endif /* !SUN & GCC */

/* Compatibility #define for V3.3, should be taken out by v4.0 - QAK */
#define DFSDnumber DFSDndatasets

#endif /* HDFI_H */

