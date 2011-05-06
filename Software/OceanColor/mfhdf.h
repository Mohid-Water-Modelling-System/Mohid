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

/* $Id: mfhdf.h,v 1.49 2001/09/26 20:50:52 bmribler Exp $ */

#ifndef _MFSD_H_
#define _MFSD_H_

#ifndef HDF
#define HDF 1
#endif

/* change this back if it causes problems on other machines than the Alhpa-QAK */
/* Reverse back to the previous way. AKC */
#include "hdf.h"
#include "netcdf.h"
#ifdef OLD_WAY
#include "local_nc.h"
#endif /* OLD_WAY */

#if defined(_MSC_VER) && !defined(_MFHDFLIB_) && !defined(_HDFLIB_)	/* Auto-link when possible */
#	define MFHDF_LIB_VER	"425"
#	if !defined(_DEBUG)
#		if !defined(_HDFDLL_)
#			define MFHDF_LIB_NAME	"HM" MFHDF_LIB_VER ".lib"
#			pragma message( "Automatic linking with the static single-threaded MFHDF library - " MFHDF_LIB_NAME )
#		else
#			define MFHDF_LIB_NAME	"HM" MFHDF_LIB_VER "m.lib"
#			pragma message( "Automatic linking with the multithreaded MFHDF DLL - " MFHDF_LIB_NAME )
#		endif
#	else
#		if !defined(_HDFDLL_)
#			define MFHDF_LIB_NAME	"HM" MFHDF_LIB_VER "d.lib"
#			pragma message( "Automatic linking with the debug static single-threaded MFHDF library - " MFHDF_LIB_NAME  )
#		else
#			define MFHDF_LIB_NAME	"HM" MFHDF_LIB_VER "md.lib"
#			pragma message( "Automatic linking with the debug multithreaded MFHDF DLL - " MFHDF_LIB_NAME  )
#		endif
#	endif
#	pragma comment(lib, MFHDF_LIB_NAME )
#endif /* defined(_MSC_VER) && !defined(_MFHDFLIB_) && !defined(_HDFLIB_) */

#define SD_UNLIMITED NC_UNLIMITED /* use this as marker for unlimited dimension */
#define SD_NOFILL    NC_NOFILL
#define SD_FILL      NC_FILL
#define SD_DIMVAL_BW_COMP   1
#define SD_DIMVAL_BW_INCOMP  0
#define SD_RAGGED    -1  /* marker for ragged dimension */

#ifdef __cplusplus
extern "C" {
#endif

HDFLIBAPI int32 SDstart
    (const char *name, int32 accs);

HDFLIBAPI intn SDend
    (int32 fid);

HDFLIBAPI intn SDfileinfo
    (int32 fid, int32 *datasets, int32 *attrs);

HDFLIBAPI int32 SDselect
    (int32 fid, int32 idx);

HDFLIBAPI intn SDgetinfo
    (int32 sdsid, char *name, int32 *rank, int32 *dimsizes, 
           int32 *nt, int32 *nattr);

#ifndef __CSTAR__
HDFLIBAPI intn SDreaddata
    (int32 sdsid, int32 *start, int32 *stride, int32 *end, void * data);
#endif

HDFLIBAPI uint16 SDgerefnumber
    (int32 sdsid);

HDFLIBAPI int32 SDnametoindex
    (int32 fid, const char *name);

HDFLIBAPI intn SDgetrange
    (int32 sdsid, void * pmax, void * pmin);

HDFLIBAPI int32 SDcreate
    (int32 fid, const char *name, int32 nt, int32 rank, int32 *dimsizes);

HDFLIBAPI int32 SDgetdimid
    (int32 sdsid, intn number);

HDFLIBAPI intn SDsetdimname
    (int32 id, const char *name);

HDFLIBAPI intn SDendaccess
    (int32 id);

HDFLIBAPI intn SDsetrange
    (int32 sdsid, void * pmax, void * pmin);

HDFLIBAPI intn SDsetattr
    (int32 id, const char *name, int32 nt, int32 count, const void * data);

HDFLIBAPI intn SDattrinfo
    (int32 id, int32 idx, char *name, int32 *nt, int32 *count);

HDFLIBAPI intn SDreadattr
    (int32 id, int32 idx, void * buf);

#ifndef __CSTAR__
HDFLIBAPI intn SDwritedata
    (int32 sdsid, int32 *start, int32 *stride, int32 *end, void * data);
#endif

HDFLIBAPI intn SDsetdatastrs
    (int32 sdsid, const char *l, const char *u, const char *f, const char *c);

HDFLIBAPI intn SDsetcal
    (int32 sdsid, float64 cal, float64 cale, float64 ioff,
               float64 ioffe, int32 nt);

HDFLIBAPI intn SDsetfillvalue
    (int32 sdsid, void * val);

HDFLIBAPI intn SDgetfillvalue
    (int32 sdsid, void * val);

HDFLIBAPI intn SDsetfillmode
    (int32 id, intn fillmode);

HDFLIBAPI intn SDgetdatastrs
    (int32 sdsid, char *l, char *u, char *f, char *c, intn len);

HDFLIBAPI intn SDgetcal
    (int32 sdsid, float64 *cal, float64 *cale, float64 *ioff, 
               float64 *ioffe, int32 *nt);

HDFLIBAPI intn SDsetdimstrs
    (int32 id, const char *l, const char *u, const char *f);

HDFLIBAPI intn SDsetdimscale
    (int32 id, int32 count, int32 nt, void * data);

HDFLIBAPI intn SDgetdimscale
    (int32 id, void * data);

HDFLIBAPI intn SDdiminfo
    (int32 id, char *name, int32 *size, int32 *nt, int32 *nattr);

HDFLIBAPI intn SDgetdimstrs
    (int32 id, char *l, char *u, char *f, intn len);

HDFLIBAPI intn SDsetexternalfile
    (int32 id, const char *filename, int32 offset);

HDFLIBAPI intn SDsetnbitdataset
    (int32 id, intn start_bit, intn bit_len, intn sign_ext, intn fill_one);

HDFLIBAPI intn SDsetcompress
    (int32 id, comp_coder_t type, comp_info *c_info);

HDFLIBAPI intn SDgetcompress
    (int32 id, comp_coder_t* type, comp_info *c_info);

HDFLIBAPI int32 SDfindattr
    (int32 id, const char *attrname);

HDFLIBAPI int32 SDidtoref
    (int32 id);

HDFLIBAPI int32 SDreftoindex
    (int32 fid, int32 ref);

HDFLIBAPI int32 SDisrecord
    (int32 id);

HDFLIBAPI intn SDiscoordvar
    (int32 id);

HDFLIBAPI intn SDsetaccesstype
    (int32 id, uintn accesstype);

HDFLIBAPI intn SDsetblocksize
    (int32 sdsid, int32 block_size);

HDFLIBAPI intn SDsetdimval_comp
    (int32 dimid, intn compt_mode);

HDFLIBAPI intn SDisdimval_bwcomp
    (int32 dimid);

HDFLIBAPI int32 SDcheckempty
    (int32 sdsid, intn *emptySDS);

/*====================== Chunking Routines ================================*/

/* For defintion of HDF_CHUNK_DEF union see hproto.h since 
   this defintion is also used by GRs. */

/******************************************************************************
 NAME
      SDsetchunk   -- make SDS a chunked SDS

 DESCRIPTION
      This routine makes the SDS a chunked SDS according to the chunk
      definition passed in.

      The dataset currently cannot be special already.  i.e. NBIT,
      COMPRESSED, or EXTERNAL. This is an Error.

      The defintion of the HDF_CHUNK_DEF union with relvant fields is:

      typedef union hdf_chunk_def_u
      {
         int32   chunk_lengths[MAX_VAR_DIMS];  Chunk lengths along each dimension

         struct 
          {   
            int32     chunk_lengths[MAX_VAR_DIMS]; Chunk lengths along each dimension
            int32     comp_type;                   Compression type 
            comp_info cinfo;                       Compression info struct 
          }comp;

      } HDF_CHUNK_DEF

      The simplist is the 'chunk_lengths' array specifiying chunk 
      lengths for each dimension where the 'flags' argument set to 
      'HDF_CHUNK';

      COMPRESSION is set by using the 'HDF_CHUNK_DEF' structure to set the
      appropriate compression information along with the required chunk lengths
      for each dimension. The compression information is the same as 
      that set in 'SDsetcompress()'. The bit-or'd'flags' argument' is set to 
      'HDF_CHUNK | HDF_COMP'.

      See the example in pseudo-C below for further usage.

      The maximum number of Chunks in an HDF file is 65,535.

      The dataset currently cannot have an UNLIMITED dimension.

      The performance of the SDxxx interface with chunking is greatly
      affected by the users access pattern over the dataset and by
      the maximum number of chunks set in the chunk cache. The cache contains 
      the Least Recently Used(LRU cache replacment policy) chunks. See the
      routine SDsetchunkcache() for further info on the chunk cache and how 
      to set the maximum number of chunks in the chunk cache. A default chunk 
      cache is always created.

      The following example shows the organization of chunks for a 2D array.
      e.g. 4x4 array with 2x2 chunks. The array shows the layout of
           chunks in the chunk array.

            4 ---------------------                                           
              |         |         |                                                 
        Y     |  (0,1)  |  (1,1)  |                                       
        ^     |         |         |                                      
        |   2 ---------------------                                       
        |     |         |         |                                               
        |     |  (0,0)  |  (1,0)  |                                      
        |     |         |         |                                           
        |     ---------------------                                         
        |     0         2         4                                       
        ---------------> X                                                       
                                                                                
        --Without compression--:
        {                                                                    
        HDF_CHUNK_DEF chunk_def;
                                                                            
        .......                                                                    
        -- Set chunk lengths --                                                    
        chunk_def.chunk_lengths[0]= 2;                                                     
        chunk_def.chunk_lengths[1]= 2; 

        -- Set Chunking -- 
        SDsetchunk(sdsid, chunk_def, HDF_CHUNK);                      
         ......                                                                  
        }                                                                           

        --With compression--:
        {                                                                    
        HDF_CHUNK_DEF chunk_def;
                                                                            
        .......                
        -- Set chunk lengths first --                                                    
        chunk_def.chunk_lengths[0]= 2;                                                     
        chunk_def.chunk_lengths[1]= 2;

        -- Set compression --
        chunk_def.comp.cinfo.deflate.level = 9;
        chunk_def.comp.comp_type = COMP_CODE_DEFLATE;

        -- Set Chunking with Compression --
        SDsetchunk(sdsid, chunk_def, HDF_CHUNK | HDF_COMP);                      
         ......                                                                  
        }                                                                           

 RETURNS
        SUCCEED/FAIL
******************************************************************************/
HDFLIBAPI intn SDsetchunk
    (int32 sdsid,             /* IN: sds access id */
     HDF_CHUNK_DEF chunk_def, /* IN: chunk definition */
     int32 flags              /* IN: flags */);

/******************************************************************************
 NAME
     SDgetchunkinfo -- get Info on SDS

 DESCRIPTION
     This routine gets any special information on the SDS. If its chunked,
     chunked and compressed or just a regular SDS. Currently it will only
     fill the array of chunk lengths for each dimension as specified in
     the 'HDF_CHUNK_DEF' union. It does not tell you the type of compression
     used or the compression parameters. You can pass in a NULL for 'chunk_def'
     if don't want the chunk lengths for each dimension.
     Additionaly if successfull it will return a bit-or'd value in 'flags' 
     indicating if the SDS is:

     Chunked                  -> flags = HDF_CHUNK
     Chunked and compressed   -> flags = HDF_CHUNK | HDF_COMP 
     Non-chunked              -> flags = HDF_NONE
  
     e.g. 4x4 array - Pseudo-C
     {
     int32   rcdims[3];
     HDF_CHUNK_DEF rchunk_def;
     int32   cflags;
     ...
     rchunk_def.chunk_lengths = rcdims;
     SDgetchunkinfo(sdsid, &rchunk_def, &cflags);
     ...
     }

 RETURNS
        SUCCEED/FAIL
******************************************************************************/
HDFLIBAPI intn SDgetchunkinfo
    (int32 sdsid,              /* IN: sds access id */
     HDF_CHUNK_DEF *chunk_def, /* IN/OUT: chunk definition */
     int32 *flags              /* IN/OUT: flags */);

/******************************************************************************
 NAME
     SDwritechunk  -- write the specified chunk to the SDS

 DESCRIPTION
     This routine writes a whole chunk of data to the chunked SDS 
     specified by chunk 'origin' for the given SDS and can be used
     instead of SDwritedata() when this information is known. This
     routine has less overhead and is much faster than using SDwritedata().

     Origin specifies the co-ordinates of the chunk according to the chunk
     position in the overall chunk array.

     'datap' must point to a whole chunk of data.

     See SDsetchunk() for a description of the organization of chunks in an SDS.

 RETURNS
        SUCCEED/FAIL
******************************************************************************/
HDFLIBAPI intn SDwritechunk
    (int32 sdsid,      /* IN: sds access id */
     int32 *origin,    /* IN: origin of chunk to write */
     const void *datap /* IN: buffer for data */);

/******************************************************************************
 NAME
     SDreadchunk   -- read the specified chunk to the SDS

 DESCRIPTION
     This routine reads a whole chunk of data from the chunked SDS
     specified by chunk 'origin' for the given SDS and can be used
     instead of SDreaddata() when this information is known. This
     routine has less overhead and is much faster than using SDreaddata().

     Origin specifies the co-ordinates of the chunk according to the chunk
     position in the overall chunk array.

     'datap' must point to a whole chunk of data.

     See SDsetchunk() for a description of the organization of chunks in an SDS.

 RETURNS
        SUCCEED/FAIL
******************************************************************************/
HDFLIBAPI intn SDreadchunk
    (int32 sdsid,      /* IN: sds access id */
     int32 *origin,    /* IN: origin of chunk to read */
     void  *datap      /* IN/OUT: buffer for data */);

/******************************************************************************
NAME
     SDsetchunkcache -- maximum number of chunks to cache 

DESCRIPTION
     Set the maximum number of chunks to cache.

     The cache contains the Least Recently Used(LRU cache replacment policy) 
     chunks. This routine allows the setting of maximum number of chunks that 
     can be cached, 'maxcache'.

     The performance of the SDxxx interface with chunking is greatly
     affected by the users access pattern over the dataset and by
     the maximum number of chunks set in the chunk cache. The number chunks 
     that can be set in the cache is process memory limited. It is a good 
     idea to always set the maximum number of chunks in the cache as the 
     default heuristic does not take into account the memory available for 
     the application.

     By default when the SDS is promoted to a chunked element the 
     maximum number of chunks in the cache 'maxcache' is set to the number of
     chunks along the last dimension.

     The values set here affects the current object's caching behaviour.

     If the chunk cache is full and 'maxcache' is greater then the
     current 'maxcache' value, then the chunk cache is reset to the new
     'maxcache' value, else the chunk cache remains at the current
     'maxcache' value.

     If the chunk cache is not full, then the chunk cache is set to the
     new 'maxcache' value only if the new 'maxcache' value is greater than the
     current number of chunks in the cache.

     Use flags argument of 'HDF_CACHEALL' if the whole object is to be cached 
     in memory, otherwise pass in zero(0). Currently you can only
     pass in zero.

    See SDsetchunk() for a description of the organization of chunks in an SDS.

RETURNS
     Returns the 'maxcache' value for the chunk cache if successful 
     and FAIL otherwise
******************************************************************************/
HDFLIBAPI intn SDsetchunkcache
    (int32 sdsid,     /* IN: sds access id */
     int32 maxcache,  /* IN: max number of chunks to cache */
     int32 flags      /* IN: flags = 0, HDF_CACHEALL */);

/* Define the FORTRAN names */

#ifndef MFSD_FNAMES
#   define  MFSD_FNAMES
#ifdef DF_CAPFNAMES
# if defined(UNIX386) || (!(defined INTEL86) && !(defined WIN32))
#   define nscstart    FNAME(SCSTART)
#   define nsfend      FNAME(SFEND)
#   define nsfendacc   FNAME(SFENDACC)
#   define nsffinfo    FNAME(SFFINFO)
#   define nsfselect   FNAME(SFSELECT)
#   define nscginfo    FNAME(SCGINFO)
#   define nscgainfo   FNAME(SCGAINFO)
#   define nscgdinfo   FNAME(SCGDINFO)
#   define nsfgcal     FNAME(SFGCAL)
#   define nsfscal     FNAME(SFSCAL)
#   define nsfgdscale  FNAME(SFGDSCALE)
#   define nsfsdscale  FNAME(SFSDSCALE)
#   define nsfgcfill   FNAME(SFGCFILL)
#   define nsfgfill    FNAME(SFGFILL)
#   define nsfscfill   FNAME(SFSCFILL)
#   define nsfsfill    FNAME(SFSFILL)
#   define nsfsflmd    FNAME(SFSFLMD)
#   define nsfgrange   FNAME(SFGRANGE)
#   define nsfsrange   FNAME(SFSRANGE)
#   define nscn2index  FNAME(SCN2INDEX)
#   define nsccreate   FNAME(SCCREATE)
#   define nscsdimstr  FNAME(SCSDIMSTR)
#   define nscsdimname FNAME(SCSDIMNAME)
#   define nscsdatstr  FNAME(SCSDATSTR)
#   define nsfdimid    FNAME(SFDIMID)
#   define nsfrcatt    FNAME(SFRCATT)
#   define nsfrnatt    FNAME(SFRNATT)
#   define nsfrattr    FNAME(SFRATTR)
#   define nsfrcdata   FNAME(SFRCDATA)
#   define nsfrdata    FNAME(SFRDATA)
#   define nsfwcdata   FNAME(SFWCDATA)
#   define nsfwdata    FNAME(SFWDATA)
#   define nscgdatstrs FNAME(SCGDATSTRS)
#   define nscgdimstrs FNAME(SCGDIMSTRS)
#   define nscscatt    FNAME(SCSCATT)
#   define nscsnatt    FNAME(SCSNATT)
#   define nscsattr    FNAME(SCSATTR)
#   define nscfattr    FNAME(SCFATTR)
#   define nsfid2ref  FNAME(SFID2REF)
#   define nsfref2index FNAME(SFREF2INDEX)
#   define nsfiscvar   FNAME(SFISCVAR)
#   define nscsextf    FNAME(SCSEXTF)
#   define nsfsacct    FNAME(SFSACCT)
#   define nsfsdmvc    FNAME(SFSDMVC)
#   define nsfisdmvc   FNAME(SFISDMVC)
#   define nsfisrcrd     FNAME(SFISRCRD)
#   define nscgichnk     FNAME(SCGICHNK)
#   define nscrcchnk     FNAME(SCRCCHNK)
#   define nscrchnk      FNAME(SCRCHNK)
#   define nscscchnk     FNAME(SCSCCHNK)
#   define nscschnk      FNAME(SCSCHNK)
#   define nscwcchnk     FNAME(SCWCCHNK) 
#   define nscwchnk      FNAME(SCWCHNK)
#   define nscscompress  FNAME(SCSCOMPRESS)
#   define nscgcompress  FNAME(SCGCOMPRESS)
#   define nsfsnbit      FNAME(SFSNBIT)
#   define nsfsblsz      FNAME(SFSBLSZ)
#   define nscchempty    FNAME(SCCHEMPTY)
# else /* Fortran PowerStation */
#   define nscstart    FNAME(SCSTART)
#   define nscend      FNAME(SCEND)
#   define nscendacc   FNAME(SCENDACC)
#   define nscfinfo    FNAME(SCFINFO)
#   define nscselect   FNAME(SCSELECT)
#   define nscginfo    FNAME(SCGINFO)
#   define nscgainfo   FNAME(SCGAINFO)
#   define nscgdinfo   FNAME(SCGDINFO)
#   define nscgcal     FNAME(SCGCAL)
#   define nscscal     FNAME(SCSCAL)
#   define nscgdscale  FNAME(SCGDSCALE)
#   define nscsdscale  FNAME(SCSDSCALE)
#   define nscgcfill   FNAME(SCGCFILL)
#   define nscgfill    FNAME(SCGFILL)
#   define nscscfill   FNAME(SCSCFILL)
#   define nscsfill    FNAME(SCSFILL)
#   define nscsflmd    FNAME(SCSFLMD)
#   define nscgrange   FNAME(SCGRANGE)
#   define nscsrange   FNAME(SCSRANGE)
#   define nscn2index  FNAME(SCN2INDEX)
#   define nsccreate   FNAME(SCCREATE)
#   define nscsdimstr  FNAME(SCSDIMSTR)
#   define nscsdimname FNAME(SCSDIMNAME)
#   define nscsdatstr  FNAME(SCSDATSTR)
#   define nscdimid    FNAME(SCDIMID)
#   define nscrcatt    FNAME(SCRCATT)
#   define nscrnatt    FNAME(SCRNATT)
#   define nscrattr    FNAME(SCRATTR)
#   define nscrcdata   FNAME(SCRCDATA)
#   define nscrdata    FNAME(SCRDATA)
#   define nscwcdata   FNAME(SCWCDATA)
#   define nscwdata    FNAME(SCWDATA)
#   define nscgdatstrs FNAME(SCGDATSTRS)
#   define nscgdimstrs FNAME(SCGDIMSTRS)
#   define nscscatt    FNAME(SCSCATT)
#   define nscsnatt    FNAME(SCSNATT)
#   define nscsattr    FNAME(SCSATTR)
#   define nscfattr    FNAME(SCFATTR)
#   define nscid2ref  FNAME(SCID2REF)
#   define nsciscvar   FNAME(SCISCVAR)
#   define nscsextf    FNAME(SCSEXTF)
#   define nscsacct    FNAME(SCSACCT)
#   define nscsdmvc    FNAME(SCSDMVC)
#   define nscisdmvc   FNAME(SCISDMVC)
#   define nscisrcrd     FNAME(SCISRCRD)
#   define nscgichnk     FNAME(SCGICHNK)
#   define nscrcchnk     FNAME(SCRCCHNK)
#   define nscrchnk      FNAME(SCRCHNK)
#   define nscscchnk     FNAME(SCSCCHNK)
#   define nscschnk      FNAME(SCSCHNK)
#   define nscwcchnk     FNAME(SCWCCHNK)
#   define nscwchnk      FNAME(SCWCHNK)
#   define nscscompress  FNAME(SCSCOMPRESS)
#   define nscgcompress  FNAME(SCGCOMPRESS)
#   define nscsnbit      FNAME(SCSNBIT)
#   define nscsblsz      FNAME(SCSBLSZ)
#   define nscselct      FNAME(SCSELCT)
#   define nscr2idx      FNAME(SCR2IDX)
#   define nscchempty    FNAME(SCCHEMPTY)
#  endif   /* Fortran PowerStation */
#else   /* DF_CAPFNAMES */
# if defined(UNIX386) || (!(defined INTEL86) && !(defined WIN32))
#   define nscstart    FNAME(scstart)
#   define nsfend      FNAME(sfend)
#   define nsfendacc   FNAME(sfendacc)
#   define nsffinfo    FNAME(sffinfo)
#   define nsfselect   FNAME(sfselect)
#   define nscginfo    FNAME(scginfo)
#   define nscgainfo   FNAME(scgainfo)
#   define nscgdinfo   FNAME(scgdinfo)
#   define nsfgcal     FNAME(sfgcal)
#   define nsfscal     FNAME(sfscal)
#   define nsfgdscale  FNAME(sfgdscale)
#   define nsfsdscale  FNAME(sfsdscale)
#   define nsfgcfill   FNAME(sfgcfill)
#   define nsfgfill    FNAME(sfgfill)
#   define nsfscfill   FNAME(sfscfill)
#   define nsfsfill    FNAME(sfsfill)
#   define nsfsflmd    FNAME(sfsflmd)
#   define nsfgrange   FNAME(sfgrange)
#   define nsfsrange   FNAME(sfsrange)
#   define nscn2index  FNAME(scn2index)
#   define nsccreate   FNAME(sccreate)
#   define nscsdimstr  FNAME(scsdimstr)
#   define nscsdimname FNAME(scsdimname)
#   define nscsdatstr  FNAME(scsdatstr)
#   define nsfdimid    FNAME(sfdimid)
#   define nsfrcatt    FNAME(sfrcatt)
#   define nsfrnatt    FNAME(sfrnatt)
#   define nsfrattr    FNAME(sfrattr)
#   define nscscatt    FNAME(scscatt)
#   define nscsnatt    FNAME(scsnatt)
#   define nscsattr    FNAME(scsattr)
#   define nscfattr    FNAME(scfattr)
#   define nsfrcdata   FNAME(sfrcdata)
#   define nsfrdata    FNAME(sfrdata)
#   define nsfwcdata   FNAME(sfwcdata)
#   define nsfwdata    FNAME(sfwdata)
#   define nscgdatstrs FNAME(scgdatstrs)
#   define nscgdimstrs FNAME(scgdimstrs)
#   define nsfid2ref   FNAME(sfid2ref)
#   define nsfref2index FNAME(sfref2index)
#   define nsfiscvar   FNAME(sfiscvar)
#   define nscsextf    FNAME(scsextf)
#   define nsfsacct    FNAME(sfsacct)
#   define nsfsdmvc    FNAME(sfsdmvc)
#   define nsfisdmvc   FNAME(sfisdmvc)
#   define nsfisrcrd     FNAME(sfisrcrd)
#   define nscgichnk     FNAME(scgichnk)
#   define nscrcchnk     FNAME(scrcchnk)
#   define nscrchnk      FNAME(scrchnk)
#   define nscscchnk     FNAME(scscchnk)
#   define nscschnk      FNAME(scschnk)
#   define nscwcchnk     FNAME(scwcchnk) 
#   define nscwchnk      FNAME(scwchnk)
#   define nscscompress  FNAME(scscompress)
#   define nscgcompress  FNAME(scgcompress)
#   define nsfsnbit      FNAME(sfsnbit)
#   define nsfsblsz      FNAME(sfsblsz)
#   define nscchempty    FNAME(scchempty)
#  else /* Powerstation */
#   define nscstart    FNAME(scstart)
#   define nscend      FNAME(scend)
#   define nscendacc   FNAME(scendacc)
#   define nscfinfo    FNAME(scfinfo)
#   define nscginfo    FNAME(scginfo)
#   define nscgainfo   FNAME(scgainfo)
#   define nscgdinfo   FNAME(scgdinfo)
#   define nscgcal     FNAME(scgcal)
#   define nscscal     FNAME(scscal)
#   define nscgdscale  FNAME(scgdscale)
#   define nscsdscale  FNAME(scsdscale)
#   define nscgcfill   FNAME(scgcfill)
#   define nscgfill    FNAME(scgfill)
#   define nscscfill   FNAME(scscfill)
#   define nscsfill    FNAME(scsfill)
#   define nscsflmd    FNAME(scsflmd)
#   define nscgrange   FNAME(scgrange)
#   define nscsrange   FNAME(scsrange)
#   define nscn2index  FNAME(scn2index)
#   define nsccreate   FNAME(sccreate)
#   define nscsdimstr  FNAME(scsdimstr)
#   define nscsdimname FNAME(scsdimname)
#   define nscsdatstr  FNAME(scsdatstr)
#   define nscdimid    FNAME(scdimid)
#   define nscrcatt    FNAME(scrcatt)
#   define nscrnatt    FNAME(scrnatt)
#   define nscrattr    FNAME(scrattr)
#   define nscscatt    FNAME(scscatt)
#   define nscsnatt    FNAME(scsnatt)
#   define nscsattr    FNAME(scsattr)
#   define nscfattr    FNAME(scfattr)
#   define nscrcdata   FNAME(scrcdata)
#   define nscrdata    FNAME(scrdata)
#   define nscwcdata   FNAME(scwcdata)
#   define nscwdata    FNAME(scwdata)
#   define nscgdatstrs FNAME(scgdatstrs)
#   define nscgdimstrs FNAME(scgdimstrs)
#   define nscid2ref   FNAME(scid2ref)
#   define nsciscvar   FNAME(sciscvar)
#   define nscsextf    FNAME(scsextf)
#   define nscsacct    FNAME(scsacct)
#   define nscsdmvc    FNAME(scsdmvc)
#   define nscisdmvc   FNAME(scisdmvc)
#   define nscisrcrd     FNAME(scisrcrd)
#   define nscgichnk     FNAME(scgichnk)
#   define nscrcchnk     FNAME(scrcchnk)
#   define nscrchnk      FNAME(scrchnk)
#   define nscscchnk     FNAME(scscchnk)
#   define nscschnk      FNAME(scschnk)
#   define nscwcchnk     FNAME(scwcchnk)
#   define nscwchnk      FNAME(scwchnk)
#   define nscgcompress  FNAME(scgcompress)
#   define nscsnbit      FNAME(scsnbit)
#   define nscsblsz      FNAME(scsblsz)
#   define nscselct      FNAME(scselct)
#   define nscr2idx      FNAME(scr2idx)
#   define nscchempty    FNAME(scchempty)
# endif /* powerstation  */
#endif /* capital */
#endif /* MFSD_FNAMES */

#ifdef __cplusplus
}
#endif

#endif /* _MFSD_H_ */
