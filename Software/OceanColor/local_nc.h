/*
 *	Copyright 1993, University Corporation for Atmospheric Research
 *      See netcdf/COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: local_nc.h,v 1.50 1999/02/25 19:02:43 koziol Exp $ */
#ifndef _LOCAL_NC_
#define _LOCAL_NC_

/*
 *	netcdf library 'private' data structures, objects and interfaces
 */

#if (defined macintosh) || (defined MPW) || (defined __MWERKS__)
#ifndef HDF
#define HDF  /* For Mac we need to define this, to avoid putting on compile line */
#endif
#define NO_SYS_XDR_INC /* use stuff in "::xdr" */
#define NO_ACCESS
#define NO_GETPID
#endif /* non command line compilers */

#include	<stddef.h> /* size_t */
#include	<stdio.h> /* FILENAME_MAX */
#if (defined MPW)
#include   <memory.h>
#endif /* MPW */

#ifndef FILENAME_MAX
#define FILENAME_MAX  255
#endif

/* Do we have systeme XDR files */
#ifndef  NO_SYS_XDR_INC 
#ifdef VMS
#    define  STDC_INCLUDES
#endif   /* VMS */
#ifdef __ultrix
#define GCC_FIX
#endif /* __ultrix */
#include	<rpc/types.h>
#ifdef __ultrix
#undef GCC_FIX
#endif /* __ultrix */
#include	<rpc/xdr.h>
#else    /* NO_SYS_XDR_INC */
#if defined(macintosh) || defined (SYMANTEC_C)
     /* For the mac reference types.h specifically
        because we don't want it to pick up the system one */
#include      "::xdr:types.h"  /* "../xdr/types.h" */
#include      "::xdr:xdr.h"    /* "../xdr/xdr.h" */
#else /* !macintosh */
#include      <types.h>  /* <types.h */
#include      <xdr.h>    /* <xdr.h> */
#endif /* !macintosh */
#endif /* NO_SYSTEM_XDR_INCLUDES */

#include	"netcdf.h" /* needed for defs of nc_type, ncvoid, ... */

/* ptr argument type in internal functions */
#define Void    char

/*
** Include HDF stuff
*/
#ifdef HDF

#include "hdf.h"
#include "vg.h"
#include "hfile.h"

#define ATTR_TAG  DFTAG_VH
#define DIM_TAG   DFTAG_VG
#define VAR_TAG   DFTAG_VG
#define DATA_TAG  DFTAG_SD
#define BOGUS_TAG ((uint16) 721)

#if 0
#define ATTRIBUTE         "Attr0.0"
#define VARIABLE          "Var0.0"
#define DIMENSION         "Dim0.0"
#define UDIMENSION        "UDim0.0"
#define DIM_VALS          "DimVal0.0" 
#define DIM_VALS01        "DimVal0.1"
#define CDF               "CDF0.0"
/* DATA is defined in DTM. Change DATA to DATA0 *
#define DATA              "Data0.0"
*/
#define DATA0             "Data0.0"
#define ATTR_FIELD_NAME   "VALUES"
#endif

#define DIMVAL_VERSION00  0  /* <dimsize> fake values */
#define DIMVAL_VERSION01  1  /* 1 elt with value of <dimsize>  */
#define BLOCK_MULT  64    /* multiplier for bytes in linked blocks */
#define MAX_BLOCK_SIZE  65536    /* maximum size of block in linked blocks */
#define BLOCK_COUNT 128   /* size of linked block pointer objects  */

#endif /* HDF */

/* from cdflib.h CDF 2.3 */
#ifndef MAX_VXR_ENTRIES
#define MAX_VXR_ENTRIES                 10
#endif /* MAX_VXR_ENTRIES */

#ifdef HDF
/* VIX record for CDF variable data storage */
typedef struct vix_t_def {
    int32              nEntries;                    /* number of entries in this vix */
    int32              nUsed;                       /* number of entries containing valid data */
    int32              firstRec[MAX_VXR_ENTRIES];   /* number of first records */
    int32              lastRec[MAX_VXR_ENTRIES];    /* number of last records */
    int32              offset[MAX_VXR_ENTRIES];     /* file offset of records */
    struct vix_t_def * next;                        /* next one in line */
} vix_t;
#endif /* HDF */

/* like, a discriminated union in the sense of xdr */
typedef struct {
	nc_type type ;		/* the discriminant */
	size_t len ;		/* the total length originally allocated */
	size_t szof ;		/* sizeof each value */
	unsigned count ;	/* length of the array */
	Void *values ;		/* the actual data */
} NC_array ;

/* Counted string for names and such */
/* 

  count is the actual size of the buffer for the string
  len is the length of the string in the buffer
  
  count != len when a string is resized to something smaller

*/
#ifdef HDF
#define NC_compare_string(s1,s2) ((s1)->hash!=(s2)->hash ? 1 : HDstrcmp((s1)->values,(s2)->values))
#endif /* HDF */

typedef struct {
	unsigned count ;
    unsigned len ; 
#ifdef HDF
    uint32 hash;        /* [non-perfect] hash value for faster comparisons */
#endif /* HDF */
	char *values ;
} NC_string ;

/* Counted array of ints for assoc list */
typedef struct {
	unsigned count ;
	int *values ;
} NC_iarray ;

/* NC dimension stucture */
typedef struct {
	NC_string *name ;
    long size ;
#ifdef HDF
    int32 dim00_compat;   /* compatible with Dim0.0 */
	int32 vgid;   /* id of the Vgroup representing this dimension */
    int32 count;  /* Number of pointers to this dimension */
#endif
} NC_dim ;

/* NC attribute */
typedef struct {
	NC_string	*name ;
	NC_array	*data ;
#ifdef HDF
	int32           HDFtype; /* it should be in NC_array *data. However, */
                             /* NC.dims and NC.vars are NC_array too. */
#endif
} NC_attr ;

typedef struct {
	char path[FILENAME_MAX + 1] ;
	unsigned flags ;
	XDR *xdrs ;
	long begin_rec ; /* (off_t) postion of the first 'record' */
	unsigned long recsize ; /* length of 'record' */
	int redefid ;
	/* below gets xdr'd */
	unsigned long numrecs ; /* number of 'records' allocated */
	NC_array *dims ;
	NC_array *attrs ;
	NC_array *vars ;
#ifdef HDF
	int32 hdf_file;
    int file_type;
    int32 vgid;
    int hdf_mode; /* mode we are attached for */
    hdf_file_t cdf_fp; /* file pointer used for CDF files */
#endif
} NC ;

/* NC variable: description and data */
typedef struct {
	NC_string *name ;
	NC_iarray *assoc ; /* user definition */
	unsigned long *shape ; /* compiled info */
	unsigned long *dsizes ; /* compiled info */
	NC_array *attrs;
	nc_type type ;		/* the discriminant */
	unsigned long len ;		/* the total length originally allocated */
	size_t szof ;		/* sizeof each value */
	long begin ;  /* seek index, often an off_t */
#ifdef HDF
	NC *cdf;    /* handle of the file where this var belongs to  */
	int32 vgid;     /* id of the variable's Vgroup */
    uint16 data_ref;  /* ref of the variable's data storage (if exists) */
    uint16 data_tag;  /* tag of the variable's data storage (if exists) */
    uint16 ndg_ref;   /* ref of ndg for this dataset */
    intn   data_offset; /* non-traditional data may not begin at 0 */
    int32  block_size;  /* size of the blocks for unlimited dim. datasets */
    int numrecs;  /* number of records this has been filled to */
    int32 aid;    /* aid for DFTAG_SD data */
    int32 HDFtype; /* type of this variable as HDF thinks */
    int32 HDFsize; /* size of this variable as HDF thinks */
    /* These next two flags control when space in the file is allocated
        for a new dataset.  They are used (currently) in SDwritedata() and
        hdf_get_vp_aid() to allocate the full length of a new fixed-size dataset
        which is not writing fill values, instead of letting them get created
        as an "appendable" dataset and probably get converted into a linked-
        block special element when they don't need to be one */
    int32   created;    /* BOOLEAN == is newly created */
    int32   set_length; /* BOOLEAN == needs length set */
    int32   is_ragged; /* BOOLEAN == is a ragged array */
    int32 * rag_list;  /* size of ragged array lines */
    int32   rag_fill;  /* last line in rag_list to be set */
    vix_t * vixHead;   /* list of VXR records for CDF data storage */
#endif
} NC_var ;

#define IS_RECVAR(vp) \
	((vp)->shape != NULL ? (*(vp)->shape == NC_UNLIMITED) : 0 )

#define netCDF_FILE  0
#define HDF_FILE     1
#define CDF_FILE     2

extern const char *cdf_routine_name ; /* defined in lerror.c */

               /*  C D F 1 */
#define	NCMAGIC	0x43444601
                       /*  C D L 1 */
#define	NCLINKMAGIC	0x43444c01

/* #ifndef HDF *//* HDF has already worked out if we have prototypes */
#ifdef HDF
#define PROTOTYPE
#endif
#undef PROTO
#ifndef NO_HAVE_PROTOTYPES 
#   define	PROTO(x)	x
#else
#   define	PROTO(x)	()
#endif
/* #endif */ /* HDF */

#ifdef __cplusplus
extern "C" {
#endif

/* If using the real netCDF library and API (use -DHAVE_NETCDF)
   need to mangle the HDF versions of netCDF API function names 
   to not conflict w/ oriinal netCDF ones */
#ifdef HAVE_NETCDF
#define nc_serror        HNAME(nc_serror)
#define NCadvise         HNAME(NCadvise)
#define NC_computeshapes HNAME(NC_computeshapes)
#define NC_xtypelen      HNAME(NC_xtypelen)
#define NC_xlen_array    HNAME(NC_xlen_array)
#define NC_xlen_attr     HNAME(NC_xlen_attr)
#define NC_xlen_cdf      HNAME(NC_xlen_cdf)
#define NC_xlen_dim      HNAME(NC_xlen_dim)
#define NC_xlen_iarray   HNAME(NC_xlen_iarray)
#define NC_xlen_string   HNAME(NC_xlen_string)
#define NC_xlen_var      HNAME(NC_xlen_var)
#define NCmemset         HNAME(NCmemset)
#define NC_arrayfill     HNAME(NC_arrayfill)
#define NC_copy_arrayvals HNAME(NC_copy_arrayvals)
#define NC_free_array    HNAME(NC_free_array)
#define NC_free_attr     HNAME(NC_free_attr)
#define NC_free_cdf      HNAME(NC_free_cdf)
#define NC_free_dim      HNAME(NC_free_dim)
#define NC_free_iarray   HNAME(NC_free_iarray)
#define NC_free_string   HNAME(NC_free_string)
#define NC_free_var      HNAME(NC_free_var)
#define NC_incr_array    HNAME(NC_incr_array)
#define NC_dimid         HNAME(NC_dimid)
#define NCcktype         HNAME(NCcktype)
#define NC_indefine      HNAME(NC_indefine)
#define xdr_cdf          HNAME(xdr_cdf)
#define xdr_numrecs      HNAME(xdr_numrecs)
#define xdr_shorts       HNAME(xdr_shorts)
#define xdr_NC_array     HNAME(xdr_NC_array)
#define xdr_NC_attr      HNAME(xdr_NC_attr)
#define xdr_NC_dim       HNAME(xdr_NC_dim)
#define xdr_NC_fill      HNAME(xdr_NC_fill)
#define xdr_NC_iarray    HNAME(xdr_NC_iarray)
#define xdr_NC_string    HNAME(xdr_NC_string)
#define xdr_NC_var       HNAME(xdr_NC_var)
#define NC_typelen       HNAME(NC_typelen)
#define NC_check_id      HNAME(NC_check_id)
#define NC_dup_cdf       HNAME(NC_dup_cdf)
#define NC_new_cdf       HNAME(NC_new_cdf)
#define NC_new_array     HNAME(NC_new_array)
#define NC_re_array      HNAME(NC_re_array)
#define NC_new_attr      HNAME(NC_new_attr)
#define NC_findattr      HNAME(NC_findattr)
#define NC_new_dim       HNAME(NC_new_dim)
#define NC_new_iarray    HNAME(NC_new_iarray)
#define NC_new_string    HNAME(NC_new_string)
#define NC_re_string     HNAME(NC_re_string)
#define NC_hlookupvar    HNAME(NC_hlookupvar)
#define NC_new_var       HNAME(NC_new_var)
#define NCvario          HNAME(NCvario)
#define NCcoordck        HNAME(NCcoordck)
#define xdr_NCvshort     HNAME(xdr_NCvshort)
#define NC_dcpy          HNAME(NC_dcpy)
#define NCxdrfile_sync   HNAME(NCxdrfile_sync)
#define NCxdrfile_create HNAME(NCxdrfile_create)
#ifdef HDF
#define NCgenio          HNAME(NCgenio)      /* from putgetg.c */
#define NC_var_shape     HNAME(NC_var_shape) /* from var.c */
#endif
#endif /* HAVE_NETCDF ie. NOT USING HDF version of netCDF ncxxx API */

extern void		nc_serror			PROTO((
	const char *fmt,
	...
)) ;
extern void		NCadvise			PROTO((
	int err,
	const char *fmt,
	...
)) ;

extern int        NC_computeshapes	PROTO((
    NC		*handle
));
extern int        NC_xtypelen		PROTO((
    nc_type	type
));
extern int        NC_xlen_array		PROTO((
    NC_array	*array
));
extern int        NC_xlen_attr		PROTO((
    NC_attr	**app
));
extern int        NC_xlen_cdf		PROTO((
    NC		*cdf
));
extern int        NC_xlen_dim		PROTO((
    NC_dim	**dpp
));
extern int        NC_xlen_iarray	PROTO((
    NC_iarray	*iarray
));
extern int        NC_xlen_string	PROTO((
    NC_string	*cdfstr
));
extern int        NC_xlen_var		PROTO((
    NC_var	**vpp
));

extern char       *NCmemset		PROTO((
    char	*s,
    int		c,
    int		n
));

extern void       NC_arrayfill		PROTO((
    void	*lo,
    size_t	len,
    nc_type	type
));
extern void       NC_copy_arrayvals	PROTO((
    char	*target,
    NC_array	*array
));
extern int       NC_free_array		PROTO((
    NC_array	*array
));
extern int       NC_free_attr		PROTO((
    NC_attr	*attr
));
extern int       NC_free_cdf		PROTO((
    NC		*handle
));
extern int       NC_free_dim		PROTO((
    NC_dim	*dim
));
extern int       NC_free_iarray	PROTO((
    NC_iarray	*iarray
));
extern int       NC_free_string	PROTO((
    NC_string	*cdfstr
));
extern int       NC_free_var		PROTO((
    NC_var	*var
));

extern Void      *NC_incr_array		PROTO((
    NC_array	*array,
    Void	*tail
));

extern int       NC_dimid               PROTO((
    NC          *handle,
    char        *name
));
extern bool_t     NCcktype		PROTO((
    nc_type	datatype
));
extern bool_t     NC_indefine		PROTO((
    int		cdfid,
    bool_t	iserr
));
extern bool_t     xdr_cdf		PROTO((
    XDR		*xdrs,
    NC		**handlep
));
extern bool_t     xdr_numrecs		PROTO((
    XDR		*xdrs,
    NC		*handle
));
extern bool_t     xdr_shorts		PROTO((
    XDR		*xdrs,
    short	*sp,
    u_int	cnt
));
extern bool_t     xdr_NC_array		PROTO((
    XDR		*xdrs,
    NC_array	**app
));
extern bool_t     xdr_NC_attr		PROTO((
    XDR		*xdrs,
    NC_attr	**app
));
extern bool_t     xdr_NC_dim		PROTO((
    XDR		*xdrs,
    NC_dim	**dpp
));
extern bool_t     xdr_NC_fill		PROTO((
    XDR		*xdrs,
    NC_var	*vp
));
extern bool_t     xdr_NC_iarray		PROTO((
    XDR		*xdrs,
    NC_iarray	**ipp
));
extern bool_t     xdr_NC_string		PROTO((
    XDR		*xdrs,
    NC_string	**spp
));
extern bool_t     xdr_NC_var		PROTO((
    XDR		*xdrs,
    NC_var	**vpp
));

extern size_t     NC_typelen		PROTO((
    nc_type	type
));

extern NC        *NC_check_id		PROTO((
    int		cdfid
));
extern NC        *NC_dup_cdf		PROTO((
    const char *name,
	int     mode,
    NC		*old
));
extern NC        *NC_new_cdf		PROTO((
    const char *name,
    int		mode
));
extern NC_array  *NC_new_array		PROTO((
    nc_type	type,
    unsigned	count,
    const void	*values
));
extern NC_array  *NC_re_array		PROTO((
    NC_array	*old,
    nc_type	type,
    unsigned	count,
    const void	*values
));
extern NC_attr  *NC_new_attr        PROTO((
    const char *name,
    nc_type type,
    unsigned count ,
    const void *values
));
extern NC_attr  **NC_findattr		PROTO((
    NC_array	**ap,
    const char	*name
));
extern NC_dim    *NC_new_dim		PROTO((
    const char	*name,
    long	size
));
extern NC_iarray *NC_new_iarray		PROTO((
    unsigned	count,
    const int		values[]
));
extern NC_string *NC_new_string		PROTO((
    unsigned	count,
    const char	*str
));
extern NC_string *NC_re_string		PROTO((
    NC_string	*old,
    unsigned	count,
    const char	*str
));
extern NC_var    *NC_hlookupvar		PROTO((
    NC		*handle,
    int		varid
));
extern NC_var    *NC_new_var		PROTO((
    const char	*name,
    nc_type	type,
    int		ndims,
    const int		*dims
));
extern int	NCvario			PROTO((
	NC *handle,
	int varid,
	const long *start,
	const long *edges,
	void *values
));
extern bool_t	NCcoordck	PROTO((
	NC *handle,
	NC_var *vp, 
	const long *coords
));
extern bool_t xdr_NCvshort    PROTO((
        XDR *xdrs,
        unsigned which,
        short *values
));
extern bool_t	NC_dcpy			PROTO((
	XDR *target,
	XDR *source,
	long nbytes
));
extern int NCxdrfile_sync
    PROTO((XDR *xdrs));

extern int NCxdrfile_create
    PROTO((XDR *xdrs,const char *path,int ncmode));

#ifdef HDF
/* this routine is found in 'xdrposix.c' */
extern void hdf_xdrfile_create
    PROTO(( XDR *xdrs, int ncop));

extern intn hdf_fill_array
    PROTO((Void  * storage,int32 len,Void  * value,int32 type));

extern intn hdf_get_data
    PROTO((NC *handle,NC_var *vp));

extern int32 hdf_get_vp_aid
    PROTO((NC *handle, NC_var *vp));

extern int hdf_map_type
    PROTO((nc_type ));

extern nc_type hdf_unmap_type
    PROTO((int ));

extern intn hdf_get_ref
    PROTO((NC *,int ));

extern intn hdf_create_dim_vdata
    PROTO((XDR *,NC *,NC_dim *));

extern intn hdf_create_compat_dim_vdata
    PROTO((XDR *xdrs, NC *handle, NC_dim *dim, int32 dimval_ver));

extern intn hdf_write_attr
    PROTO((XDR *,NC *,NC_attr **));

extern int32 hdf_write_dim
    PROTO((XDR *,NC *,NC_dim **,int32));

extern int32 hdf_write_var
    PROTO((XDR *,NC *,NC_var **));

extern intn hdf_write_xdr_cdf
    PROTO((XDR *,NC **));

extern intn hdf_conv_scales
    PROTO((NC **));

extern intn hdf_read_dims
    PROTO((XDR *,NC *,int32 ));

extern NC_array *hdf_read_attrs
    PROTO((XDR *,NC *,int32 ));

extern intn hdf_read_vars
    PROTO((XDR *,NC *,int32 ));

extern intn hdf_read_xdr_cdf
    PROTO((XDR *,NC **));

extern intn hdf_xdr_cdf
    PROTO((XDR *,NC **));

extern intn hdf_vg_clobber
    PROTO((NC *,int ));

extern intn hdf_cdf_clobber
    PROTO((NC *));

extern intn hdf_close
    PROTO((NC *));

extern intn hdf_read_sds_dims
    PROTO((NC *));

extern intn hdf_read_sds_cdf
    PROTO((XDR *,NC **));

extern intn SDPfreebuf PROTO((void));

extern intn NCgenio
    PROTO((NC *handle, int varid, const long *start, const long *count,
        const long *stride, const long *imap, void *values));

extern intn NC_var_shape
    PROTO((NC_var *var,NC_array *dims));

/* CDF stuff. don't need anymore? -GV */
extern nc_type cdf_unmap_type
    PROTO((int type));

extern bool_t nssdc_read_cdf
    PROTO((XDR *xdrs, NC **handlep));

extern bool_t nssdc_write_cdf
   PROTO((XDR *xdrs, NC **handlep));

extern bool_t nssdc_xdr_cdf
    PROTO((XDR *xdrs, NC **handlep));

#endif /* HDF */

#ifdef __cplusplus
}
#endif

#endif /* _LOCAL_NC_ */
