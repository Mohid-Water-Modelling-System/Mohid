
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

/* $Id: hdf.h,v 1.61 2001/07/06 02:55:48 bmribler Exp $ */

#ifndef HDF_H
#define HDF_H

#include "hdfi.h"
#include "hlimits.h"

/* Internal DF structure */
typedef struct
  {
      uint16      tag;          /* tag of element */
      uint16      ref;          /* ref of element */
  }
DFdi;

/* For annotations */
/* enumerated types of the varous annotation types */
typedef enum 
{ 
    AN_UNDEF = -1,
    AN_DATA_LABEL = 0, /* Data label */
    AN_DATA_DESC,      /* Data description */
    AN_FILE_LABEL,     /* File label */
    AN_FILE_DESC       /* File description */
} ann_type;

/* internal file access codes */

#define DFACC_READ 1
#define DFACC_WRITE 2
#define DFACC_CREATE 4
#define DFACC_ALL 7

#define DFACC_RDONLY 1
#define DFACC_RDWR 3
#define DFACC_CLOBBER 4

/* New file access codes (for Hstartaccess only, currently) */
#define DFACC_BUFFER 8  /* buffer the access to this AID */
#define DFACC_APPENDABLE 0x10 /* make this AID appendable */
#define DFACC_CURRENT 0x20 /* start looking for a tag/ref from the current */
			   /* location in the DD list (useful for continued */
			   /* searching ala findfirst/findnext) */

/* External Element File access mode */
/* #define DFACC_CREATE 4	is for creating new external element file */
#define DFACC_OLD	1	/* for accessing existing ext. element file */

/* The magic cookie for Hcache to cache all files */
#define CACHE_ALL_FILES (-2)

/* File access modes */
/* 001--007 for different serial modes */
/* 011--017 for different parallel modes */

#define DFACC_DEFAULT   000
#define DFACC_SERIAL    001
#define DFACC_PARALLEL  011

/* used by Hnextread to determine where to start searching for the
   next tag/ref to read */

#define DF_START 0
#define DF_CURRENT 1
#define DF_END 2

/* Used by Hfind to determine the direction to search for tag/ref's in the */
/* file. */

#define DF_FORWARD  1
#define DF_BACKWARD 2

/* return code - since some unix/c routines use 0 and -1 as their return
   code, and some assumption had been made in the code about that, it is
   important to keep these constants the same values.  For explicitly
   boolean functions, use TRUE and FALSE */

#define SUCCEED 0
#define FAIL (-1)

/* boolean values,  reminder: NEVER compare with numeric values */

#ifndef FALSE
#   define FALSE 0
#endif
#ifndef TRUE
#   define TRUE (!FALSE)
#endif

/* macros */
#define STREQ(s, t) (HDstrcmp((s), (t)) == 0)
#define NSTREQ(s, t, n) (HDstrncmp((s), (t), (n)) == 0)

/*
 * Macros used for variable and function scoping in code.....
 */
#ifndef EXPORT
#define EXPORT
#endif

#ifndef PRIVATE
#define PRIVATE static
#endif

/* Include the Number-type definitions */
#include "hntdefs.h"

/* Include the Tag definitions */
#include "htags.h"

/*
   * interlacing supported by the vset.
 */

#define FULL_INTERLACE  0
#define NO_INTERLACE    1

/* type for File ID to send to Hlevel from Vxx interface */
typedef int32 HFILEID;

typedef intn (*hdf_termfunc_t)(void);   /* termination function typedef */

/* .................................................................. */

/* Publically accessible functions declarations.  This includes all the
   functions that are used by application programs.  */

#include "hbitio.h"
#include "hcomp.h"
#include "herr.h"
#include "hproto.h"
#include "vg.h"         /* Add the Vgroup/Vdata header so the users don't have to */
#include "mfgr.h"       /* Add the GR header so the users don't have to */

/* For Pablo Instrumentation */
#ifdef HAVE_PABLO
#include "ProcIDs.h"
#include "trace.h"
#endif /* PABLO */

/* these may eventaully evolve into real-life functions but not yet */
#define HDFopen(f,a,d)      Hopen((f), (a), (d))
#define HDFclose(f,a,d)     Hclose((f), (a), (d))
#define Vstart(f)           Vinitialize((f))
#define Vend(f)             Vfinish((f))

/* Misc. macros for backward compability */
#define HDgettagname(tag)   HDgettagdesc(tag)

/* This is also defined in fmpio.h */
#define MP_PAGEALL    0x01  /* page the whole file i.e. no limit on 'maxcache'*/

#endif /* HDF_H */

