/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF.  The full HDF copyright notice, including       *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at      *
 * http://hdfgroup.org/products/hdf4/doc/Copyright.html.  If you do not have *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* $Id: mfhdf.h 5167 2009-01-26 01:55:18Z epourmal $ */

#ifndef _MFSD_H_
#define _MFSD_H_

#ifndef HDF
#define HDF 1
#endif

/* change this back if it causes problems on other machines than the Alhpa-QAK */
/* Reverse back to the previous way. AKC */
#include "hdf.h"
#ifdef H4_HAVE_NETCDF
#include "netcdf.h"
#else
#include "hdf4_netcdf.h"
#endif
#ifdef OLD_WAY
#include "local_nc.h"
#endif /* OLD_WAY */

#include "mfhdfi.h"

#define SD_UNLIMITED NC_UNLIMITED /* use this as marker for unlimited dimension */
#define SD_NOFILL    NC_NOFILL
#define SD_FILL      NC_FILL
#define SD_DIMVAL_BW_COMP   1
#define SD_DIMVAL_BW_INCOMP  0
#define SD_RAGGED    -1  /* marker for ragged dimension */

/* used to indicate the type of the variable at an index */
typedef struct hdf_varlist
{
    int32 var_index;     /* index of the current variable */
    hdf_vartype_t var_type; /* type of a variable (IS_SDSVAR, IS_CRDVAR, or UNKNOWN */
} hdf_varlist_t;

/* enumerated types for various types of ids in SD interface */
typedef enum
{
    NOT_SDAPI_ID = -1,	/* not an SD API id */
    SD_ID = 0,		/* SD id */
    SDS_ID,		/* SDS id */
    DIM_ID		/* Dimension id */
} hdf_idtype_t;

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

HDFLIBAPI intn SDnametoindices
    (int32 fid, const char *name, hdf_varlist_t *var_list);

HDFLIBAPI intn SDgetnumvars_byname
    (int32 fid, const char *name, int32 *n_vars);

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

HDFLIBAPI intn SDgetcompinfo
    (int32 id, comp_coder_t* type, comp_info *c_info);

HDFLIBAPI intn SDgetcomptype
    (int32 id, comp_coder_t* type);

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

HDFLIBAPI hdf_idtype_t SDidtype
    (int32 an_id);

HDFLIBAPI intn SDreset_maxopenfiles
    (intn req_max);

HDFLIBAPI intn SDget_maxopenfiles
    (intn *curr_max, intn *sys_limit);

HDFLIBAPI intn SDget_numopenfiles
    ();

HDFLIBAPI intn SDgetdatasize
    (int32 sdsid, int32 *comp_size, int32 *uncomp_size);

HDFLIBAPI intn SDgetfilename
    (int32 fid, char *filename);

HDFLIBAPI intn SDgetnamelen
    (int32 sdsid, uint16 *name_len);

HDFLIBAPI intn HDiscdf
    (const char *filename);

HDFLIBAPI intn HDisnetcdf
    (const char *filename);

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
         int32   chunk_lengths[H4_MAX_VAR_DIMS];  Chunk lengths along each dimension

         struct 
          {   
            int32     chunk_lengths[H4_MAX_VAR_DIMS]; Chunk lengths along each dimension
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

#ifdef _WIN32 /* windows specific cstub functions*/

#define nscend          H4_F77_FUNC(scend,SCEND)
#define nscendacc       H4_F77_FUNC(scendacc,SCENDACC)
#define nscfinfo        H4_F77_FUNC(scfinfo,SCFINFO)
#define nscselct       H4_F77_FUNC(scselect,SCSELCT)
#define nscdimid        H4_F77_FUNC(scdimid,SCDIMID)
#define nscgcal         H4_F77_FUNC(scgcal,SCGCAL)
#define nscscal         H4_F77_FUNC(scscal, SCSCAL)
#define nscgdscale      H4_F77_FUNC(scgdscale, SCGDSCALE)
#define nscsdscale      H4_F77_FUNC(scsdscale, SCSDSCALE)
#define nscgcfill       H4_F77_FUNC(scgcfill, SCGCFILL)
#define nscgfill        H4_F77_FUNC(scgfill, SCGFILL)
#define nscscfill       H4_F77_FUNC(scscfill, SCSCFILL)
#define nscsfill        H4_F77_FUNC(scsfill, SCSFILL)
#define nscsrange       H4_F77_FUNC(scsrange, SCSRANGE)
#define nscgrange       H4_F77_FUNC(scgrange, SCGRANGE)
#define nscrcatt        H4_F77_FUNC(scrcatt, SCRCATT)
#define nscrnatt        H4_F77_FUNC(scrnatt, SCRNATT)
#define nscrattr        H4_F77_FUNC(scrattr, SCRATTR)
#define nscrcdata       H4_F77_FUNC(scrcdata, SCRCDATA)
#define nscrdata        H4_F77_FUNC(scrdata, SCRDATA)
#define nscwcdata       H4_F77_FUNC(scwcdata, SCWCDATA)
#define nscwdata        H4_F77_FUNC(scwdata, SCWDATA)
#define nscid2ref       H4_F77_FUNC(scid2ref, SCID2REF)
#define nscref2index    H4_F77_FUNC(scref2index,SCREF2INDEX)
#define nscr2idx        H4_F77_FUNC(scr2idx,SCR2IDX)
#define nsciscvar       H4_F77_FUNC(sciscvar,SCISCVAR)
#define nscsnbit        H4_F77_FUNC(scsnbit,SCSNBIT)
#define nscsacct        H4_F77_FUNC(scsacct,SCSACCT)
#define nscsdmvc        H4_F77_FUNC(scsdmvc,SCSDMVC)
#define nscisdmvc       H4_F77_FUNC(scisdmvc,SCISDMVC)
#define nscsflmd        H4_F77_FUNC(scsflmd,SCSFLMD)
#define nscisrcrd       H4_F77_FUNC(scisrcrd,SCISRCRD)
#define nscsblsz        H4_F77_FUNC(scsblsz,SCSBLSZ)
#endif

#   define nscstart      H4_F77_FUNC(scstart, SCSTART)
#   define nsfend        H4_F77_FUNC(sfend, SFEND)
#   define nsfendacc     H4_F77_FUNC(sfendacc, SFENDACC)
#   define nsffinfo      H4_F77_FUNC(sffinfo, SFFINFO)
#   define nsfselect     H4_F77_FUNC(sfselect, SFSELECT)
#   define nscginfo      H4_F77_FUNC(scginfo, SCGINFO)
#   define nscgainfo     H4_F77_FUNC(scgainfo, SCGAINFO)
#   define nscgdinfo     H4_F77_FUNC(scgdinfo, SCGDINFO)
#   define nsfgcal       H4_F77_FUNC(sfgcal, SFGCAL)
#   define nsfscal       H4_F77_FUNC(sfscal, SFSCAL)
#   define nsfgdscale    H4_F77_FUNC(sfgdscale, SFGDSCALE)
#   define nsfsdscale    H4_F77_FUNC(sfsdscale, SFSDSCALE)
#   define nsfgcfill     H4_F77_FUNC(sfgcfill, SFGCFILL)
#   define nsfgfill      H4_F77_FUNC(sfgfill, SFGFILL)
#   define nsfscfill     H4_F77_FUNC(sfscfill, SFSCFILL)
#   define nsfsfill      H4_F77_FUNC(sfsfill, SFSFILL)
#   define nsfsflmd      H4_F77_FUNC(sfsflmd, SFSFLMD)
#   define nsfgrange     H4_F77_FUNC(sfgrange, SFGRANGE)
#   define nsfsrange     H4_F77_FUNC(sfsrange, SFSRANGE)
#   define nscn2index    H4_F77_FUNC(scn2index, SCN2INDEX)
#   define nsccreate     H4_F77_FUNC(sccreate, SCCREATE)
#   define nscsdimstr    H4_F77_FUNC(scsdimstr, SCSDIMSTR)
#   define nscsdimname   H4_F77_FUNC(scsdimname, SCSDIMNAME)
#   define nscsdatstr    H4_F77_FUNC(scsdatstr, SCSDATSTR)
#   define nsfdimid      H4_F77_FUNC(sfdimid, SFDIMID)
#   define nsfrcatt      H4_F77_FUNC(sfrcatt, SFRCATT)
#   define nsfrnatt      H4_F77_FUNC(sfrnatt, SFRNATT)
#   define nsfrattr      H4_F77_FUNC(sfrattr, SFRATTR)
#   define nsfrcdata     H4_F77_FUNC(sfrcdata, SFRCDATA)
#   define nsfrdata      H4_F77_FUNC(sfrdata, SFRDATA)
#   define nsfwcdata     H4_F77_FUNC(sfwcdata, SFWCDATA)
#   define nsfwdata      H4_F77_FUNC(sfwdata, SFWDATA)
#   define nscgdatstrs   H4_F77_FUNC(scgdatstrs, SCGDATSTRS)
#   define nscgdimstrs   H4_F77_FUNC(scgdimstrs, SCGDIMSTRS)
#   define nscscatt      H4_F77_FUNC(scscatt, SCSCATT)
#   define nscsnatt      H4_F77_FUNC(scsnatt, SCSNATT)
#   define nscsattr      H4_F77_FUNC(scsattr, SCSATTR)
#   define nscfattr      H4_F77_FUNC(scfattr, SCFATTR)
#   define nsfid2ref     H4_F77_FUNC(sfid2ref, SFID2REF)
#   define nsfref2index  H4_F77_FUNC(sfref2index, SFREF2INDEX)
#   define nsfiscvar     H4_F77_FUNC(sfiscvar, SFISCVAR)
#   define nscsextf      H4_F77_FUNC(scsextf, SCSEXTF)
#   define nsfsacct      H4_F77_FUNC(sfsacct, SFSACCT)
#   define nsfsdmvc      H4_F77_FUNC(sfsdmvc, SFSDMVC)
#   define nsfisdmvc     H4_F77_FUNC(sfisdmvc, SFISDMVC)
#   define nsfisrcrd     H4_F77_FUNC(sfisrcrd, SFISRCRD)
#   define nscgichnk     H4_F77_FUNC(scgichnk, SCGICHNK)
#   define nscrcchnk     H4_F77_FUNC(scrcchnk, SCRCCHNK)
#   define nscrchnk      H4_F77_FUNC(scrchnk, SCRCHNK)
#   define nscscchnk     H4_F77_FUNC(scscchnk, SCSCCHNK)
#   define nscschnk      H4_F77_FUNC(scschnk, SCSCHNK)
#   define nscwcchnk     H4_F77_FUNC(scwcchnk, SCWCCHNK) 
#   define nscwchnk      H4_F77_FUNC(scwchnk, SCWCHNK)
#   define nscscompress  H4_F77_FUNC(scscompress, SCSCOMPRESS)
#   define nscgcompress  H4_F77_FUNC(scgcompress, SCGCOMPRESS)
#   define nsfsnbit      H4_F77_FUNC(sfsnbit, SFSNBIT)
#   define nsfsblsz      H4_F77_FUNC(sfsblsz, SFSBLSZ)
#   define nscchempty    H4_F77_FUNC(scchempty, SCCHEMPTY)
#   define nscgetfname   H4_F77_FUNC(scgetfname,SCGETFNAME)
#   define nscgetnamelen H4_F77_FUNC(scgetnamelen,SCGETNAMELEN)
#   define nscidtype     H4_F77_FUNC(scidtype,SCIDTYPE)
#   define nscgnvars_byname H4_F77_FUNC_(scgnvars_byname,SCGNVARS_BYNAME)
#   define nscn2indices  H4_F77_FUNC(scn2indices,SCN2INDICES)
#   define nscgmaxopenf  H4_F77_FUNC(scgmaxopenf,SCGMAXOPENF)
#   define nscgnumopenf  H4_F77_FUNC(scgnumopenf,SCGNUMOPENF)
#   define nscrmaxopenf  H4_F77_FUNC(scrmaxopenf,SCRMAXOPENF)


#ifdef __cplusplus
}
#endif

#endif /* _MFSD_H_ */
