!
! netcdf version 3.6.1 Beta1 fortran interface
!
! for using by Compaq Visual Fortran compiler 6.5, 6.6 and Intel Visual Fortran compiler 9.*, 10.* in Windows2000 / XP:
! revised date  : Jul 19, 2007
! revised date  : Jan 7, 2007
! revised author: Dr. YAN Haoming
! email         : wfllib@gmail.com
! homepage      : http://yhmyhm.googlepages.com/fortran.htm
!

!
!DEC$OBJCOMMENT LIB:'netcdf.LIB'
!DEC$OBJCOMMENT LIB:'netcdfs.LIB'
!

!
!DEC$ IF DEFINED(_DF_VERSION_)
!DEC$ ELSE
!DEC$ DEFINE _DF_VERSION_ =660
!DEC$ ENDIF



!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)

!
! default fill values:
!
      integer     nf_fill_byte
      integer     nf_fill_int1
      integer     nf_fill_char
      integer     nf_fill_short
      integer     nf_fill_int2
      integer     nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690e+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 512)
      parameter (nf_max_attrs = 4096)
      parameter (nf_max_vars = 4096)
      parameter (nf_max_name = 128)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
			INTERFACE
			
			FUNCTION NF_INQ_LIBVERS()
!DEC$ IF( _DF_VERSION_ > 650 )      
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_LIBVERS' :: NF_INQ_LIBVERS
!DEC$ ENDIF			
      character*80   nf_inq_libvers
      END FUNCTION   nf_inq_libvers
      
      FUNCTION NF_STRERROR(NCERR)
!DEC$ IF( _DF_VERSION_ > 650 )      
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_STRERROR' :: NF_STRERROR
!DEC$ ENDIF
      character*80   nf_strerror
   integer   ncerr
      END FUNCTION   nf_strerror
      
      FUNCTION NF_ISSYSERR(NCERR)
!DEC$ IF( _DF_VERSION_ > 650 )      
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_ISSYSERR' :: NF_ISSYSERR
!DEC$ ENDIF
      logical  nf_issyserr
			integer  ncerr
      END FUNCTION nf_issyserr
      
! control routines:
!
      FUNCTION NF_INQ_BASE_PE(NCID,PE)
!DEC$ IF( _DF_VERSION_ > 650 )      
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_BASE_PE' :: NF_INQ_BASE_PE
!DEC$ ENDIF
      integer  nf_inq_base_pe
			integer  ncid
      integer  pe
      END FUNCTION nf_inq_base_pe
      
      FUNCTION NF_SET_BASE_PE(NCID,PE)
!DEC$ IF( _DF_VERSION_ > 650 )      
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_SET_BASE_PE' :: NF_SET_BASE_PE
!DEC$ ENDIF
      integer  nf_set_base_pe
			integer  ncid,PE
      END FUNCTION nf_set_base_pe
      
      FUNCTION NF_CREATE(PATH,CMODE,NCID)
!DEC$ IF( _DF_VERSION_ > 650 )      
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_CREATE' :: NF_CREATE
!DEC$ ENDIF
      integer  nf_create
			character(*)  path
	    integer cmode, ncid
      END FUNCTION nf_create
      
      FUNCTION NF__CREATE(PATH, CMODE, INITIALSZ, CHUNKSIZEHINT, NCID)
!DEC$ IF( _DF_VERSION_ > 650 )      
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF__CREATE' :: NF__CREATE
!DEC$ ENDIF
      integer  nf__create
			character*(*)       path
			integer  cmode, initialsz, chunksizehint, ncid
			END FUNCTION  nf__create
			
			FUNCTION NF__CREATE_MP(PATH, CMODE, INITIALSZ, BASEPE, CHUNKSIZEHINT, NCID)
!DEC$ IF( _DF_VERSION_ > 650 )			
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF__CREATE_MP' :: NF__CREATE_MP
!DEC$ ENDIF
      integer nf__create_mp
      character*(*)       path
			integer  cmode, initialsz, basepe, chunksizehint, ncid
			END FUNCTION nf__create_mp
			
			FUNCTION NF_OPEN(PATH, MODE, NCID)
!DEC$ IF( _DF_VERSION_ > 650 )			
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_OPEN' :: NF_OPEN
!DEC$ ENDIF
			integer  nf_open
			character(*)  path
      integer  mode, ncid
      END FUNCTION nf_open
      
      FUNCTION NF__OPEN(PATH, MODE, CHUNKSIZEHINT, NCID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF__OPEN' :: NF__OPEN
!DEC$ ENDIF
      integer  nf__open
      character(*) path
      integer  mode,chunksizehint, ncid
      END FUNCTION  nf__open
      
      FUNCTION NF__OPEN_MP(PATH, MODE, BASEPE, CHUNKSIZEHINT, NCID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF__OPEN_MP' :: NF__OPEN_MP
!DEC$ ENDIF
      integer nf__open_mp
      character(*)   path
      integer   mode, basepe, chunksizehint, ncid
      END FUNCTION nf__open_mp
      
      FUNCTION NF_SET_FILL(NCID, FILLMODE, OLD_MODE)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_SET_FILL' :: NF_SET_FILL
!DEC$ ENDIF
      integer  nf_set_fill
      integer  ncid, fillmode, old_mode
      END FUNCTION  nf_set_fill
      
      FUNCTION NF_SET_DEFAULT_FORMAT(FORMAT, OLD_FORMAT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_SET_DEFAULT_FORMAT' :: NF_SET_DEFAULT_FORMAT
!DEC$ ENDIF
      integer  nf_set_default_format
			integer  format,  old_format
      END FUNCTION  nf_set_default_format
      
      FUNCTION NF_REDEF(NCID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_REDEF' :: NF_REDEF
!DEC$ ENDIF
      integer  nf_redef
      integer  ncid
      END FUNCTION   nf_redef
      
      FUNCTION NF_ENDDEF(NCID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_ENDDEF' :: NF_ENDDEF
!DEC$ ENDIF
      integer  nf_enddef
      integer  ncid
      END FUNCTION   nf_enddef
      
      FUNCTION NF__ENDDEF(NCID, H_MINFREE, V_ALIGN, V_MINFREE, R_ALIGN)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF__ENDDEF' :: NF__ENDDEF
!DEC$ ENDIF
      integer  nf__enddef
      integer  ncid, h_minfree, v_align, v_minfree, r_align
      END FUNCTION nf__enddef
      
      FUNCTION NF_SYNC(NCID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_SYNC' :: NF_SYNC
!DEC$ ENDIF
      integer  nf_sync
      integer  ncid
      END FUNCTION  nf_sync
      
      FUNCTION NF_ABORT(NCID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_ABORT' :: NF_ABORT
!DEC$ ENDIF
      integer   nf_abort
      integer   ncid
      END FUNCTION  nf_abort
      
      FUNCTION NF_CLOSE(NCID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_CLOSE' :: NF_CLOSE
!DEC$ ENDIF
      integer nf_close
      integer ncid
      END FUNCTION  nf_close
      
      FUNCTION NF_DELETE(NCID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_DELETE' :: NF_DELETE
!DEC$ ENDIF
      integer  nf_delete
      character(*) ncid
      END FUNCTION nf_delete

!
! general inquiry routines:
!

      FUNCTION NF_INQ(NCID, NDIMS, NVARS, NGATTS, UNLIMDIMID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ' :: NF_INQ
!DEC$ ENDIF
      integer  nf_inq
      integer  ncid, ndims, nvars, ngatts, unlimdimid
      END FUNCTION  nf_inq
      
      FUNCTION NF_INQ_NDIMS(NCID, NDIMS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_NDIMS' :: NF_INQ_NDIMS
!DEC$ ENDIF
      integer nf_inq_ndims
			integer ncid, ndims
      END FUNCTION nf_inq_ndims
      
      FUNCTION NF_INQ_NVARS(NCID, NVARS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_NVARS' :: NF_INQ_NVARS
!DEC$ ENDIF
      integer  nf_inq_nvars
			integer  ncid, nvars
      END FUNCTION nf_inq_nvars

      FUNCTION NF_INQ_NATTS(NCID, NGATTS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_NATTS' :: NF_INQ_NATTS
!DEC$ ENDIF
      integer   nf_inq_natts
      integer  ncid, ngatts
      END FUNCTION nf_inq_natts

      FUNCTION NF_INQ_UNLIMDIM(NCID, UNLIMDIMID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_UNLIMDIM' :: NF_INQ_UNLIMDIM
!DEC$ ENDIF
      integer   nf_inq_unlimdim
      integer   ncid, unlimdimid
      END FUNCTION nf_inq_unlimdim

!
! dimension routines:
!

      FUNCTION NF_DEF_DIM(NCID, NAME, LEN, DIMID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_DEF_DIM' :: NF_DEF_DIM
!DEC$ ENDIF
      integer   nf_def_dim
      integer   ncid
      character(*)  name
      integer  len, dimid
      END FUNCTION nf_def_dim

      FUNCTION NF_INQ_DIMID(NCID,NAME, DIMID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_DIMID' :: NF_INQ_DIMID
!DEC$ ENDIF
      integer   nf_inq_dimid
      integer   ncid
      character(*) name
      integer  dimid
      END FUNCTION nf_inq_dimid

      FUNCTION NF_INQ_DIM(NCID, DIMID, NAME, LEN)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_DIM' :: NF_INQ_DIM
!DEC$ ENDIF
      integer   nf_inq_dim
      integer   ncid, dimid
      character(*) name
      integer   len
      END FUNCTION nf_inq_dim

      FUNCTION NF_INQ_DIMNAME(NCID, DIMID, NAME)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_DIMNAME' :: NF_INQ_DIMNAME
!DEC$ ENDIF
      integer   nf_inq_dimname
      integer   ncid, dimid
      character(*) name
      END FUNCTION nf_inq_dimname

      FUNCTION NF_INQ_DIMLEN(NCID, DIMID, LEN)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_DIMLEN' :: NF_INQ_DIMLEN
!DEC$ ENDIF
      integer   nf_inq_dimlen
      integer   ncid, dimid,len
      END FUNCTION nf_inq_dimlen

      FUNCTION NF_RENAME_DIM(NCID, DIMID, NAME)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_RENAME_DIM' :: NF_RENAME_DIM
!DEC$ ENDIF
      integer   nf_rename_dim
      integer   ncid, dimid
      character(*) name
      END FUNCTION nf_rename_dim

!
! general attribute routines:
!

      FUNCTION NF_INQ_ATT(NCID, VARID, NAME, XTYPE, LEN)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_ATT' :: NF_INQ_ATT
!DEC$ ENDIF
      integer   nf_inq_att
      integer   ncid, varid
      character(*) name
      integer   xtype, len
      END FUNCTION nf_inq_att

      FUNCTION NF_INQ_ATTID(NCID, VARID, NAME, ATTNUM)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_ATTID' :: NF_INQ_ATTID
!DEC$ ENDIF
      integer   nf_inq_attid
      integer   ncid, varid
      character(*) name
      integer   attnum
      END FUNCTION nf_inq_attid

      FUNCTION NF_INQ_ATTTYPE(NCID, VARID, NAME, XTYPE)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_ATTTYPE' :: NF_INQ_ATTTYPE
!DEC$ ENDIF
      integer   nf_inq_atttype
      integer   ncid, varid
      character(*) name
      integer   xtype
      END FUNCTION nf_inq_atttype

      FUNCTION NF_INQ_ATTLEN(NCID, VARID, NAME, LEN)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_ATTLEN' :: NF_INQ_ATTLEN
!DEC$ ENDIF
      integer   nf_inq_attlen
      integer   ncid, varid
      character(*) name
      integer   len
      END FUNCTION nf_inq_attlen

      FUNCTION NF_INQ_ATTNAME(NCID, VARID, ATTNUM, NAME)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_ATTNAME' :: NF_INQ_ATTNAME
!DEC$ ENDIF
      integer   nf_inq_attname
      integer   ncid, varid, attnum
      character(*) name
      END FUNCTION nf_inq_attname

      FUNCTION NF_COPY_ATT(NCID_IN, VARID_IN, NAME, NCID_OUT, VARID_OUT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_COPY_ATT' :: NF_COPY_ATT
!DEC$ ENDIF
      integer   nf_copy_att
      integer   ncid_in, varid_in
      character(*) name
      integer   ncid_out, varid_out
      END FUNCTION nf_copy_att

      FUNCTION NF_RENAME_ATT(NCID, VARID, CURNAME, NEWNAME)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_RENAME_ATT' :: NF_RENAME_ATT
!DEC$ ENDIF
      integer   nf_rename_att
      integer   ncid,  varid
      character(*) curname, newname
      END FUNCTION nf_rename_att

      FUNCTION NF_DEL_ATT(NCID, VARID, NAME)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_DEL_ATT' :: NF_DEL_ATT
!DEC$ ENDIF
      integer   nf_del_att
      integer   ncid, varid
      character(*) name
      END FUNCTION nf_del_att

!
! attribute put/get routines:
!

      FUNCTION NF_PUT_ATT_TEXT(NCID, VARID, NAME, LEN, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_ATT_TEXT' :: NF_PUT_ATT_TEXT
!DEC$ ENDIF
      integer   nf_put_att_text
      integer   ncid, varid
      character(*) name
      integer   len
      character(*) text
      END FUNCTION nf_put_att_text

      FUNCTION NF_GET_ATT_TEXT(NCID, VARID, NAME, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_ATT_TEXT' :: NF_GET_ATT_TEXT
!DEC$ ENDIF
      integer   nf_get_att_text
      integer   ncid, varid
      character(*) name
      character(*) text
      END FUNCTION nf_get_att_text

      FUNCTION NF_PUT_ATT_INT1(NCID, VARID, NAME, XTYPE, LEN, i1vals )
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_ATT_INT1' :: NF_PUT_ATT_INT1
!DEC$ ENDIF
      integer   nf_put_att_int1
      integer   ncid,varid
      character(*) name
      integer   xtype,len
      integer*1  i1vals(1)
      END FUNCTION nf_put_att_int1

      FUNCTION NF_GET_ATT_INT1(NCID, VARID, NAME, I1VALS)

!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_ATT_INT1' :: NF_GET_ATT_INT1
!DEC$ ENDIF
      integer   nf_get_att_int1
      integer   ncid, varid
      character(*) name
      integer*1   i1vals(1)
      END FUNCTION nf_get_att_int1

      FUNCTION NF_PUT_ATT_INT2(NCID, VARID, NAME, XTYPE, LEN, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_ATT_INT2' :: NF_PUT_ATT_INT2
!DEC$ ENDIF
      integer   nf_put_att_int2
      integer   ncid, varid
      character(*) name
      integer   xtype, len
      integer*2   i2vals(1)
      END FUNCTION nf_put_att_int2

      FUNCTION NF_GET_ATT_INT2(NCID, VARID, NAME, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_ATT_INT2' :: NF_GET_ATT_INT2
!DEC$ ENDIF
      integer   nf_get_att_int2
      integer   ncid, varid
      character(*) name
      integer*2   i2vals(1)
      END FUNCTION nf_get_att_int2

      FUNCTION NF_PUT_ATT_INT(NCID, VARID, NAME, XTYPE, LEN, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_ATT_INT' :: NF_PUT_ATT_INT
!DEC$ ENDIF
      integer   nf_put_att_int
      integer   ncid, varid
      character(*) name
      integer   xtype,len
      integer   ivals(1)
      END FUNCTION nf_put_att_int

      FUNCTION NF_GET_ATT_INT(NCID, VARID, NAME, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_ATT_INT' :: NF_GET_ATT_INT
!DEC$ ENDIF
      integer   nf_get_att_int
      integer   ncid, varid
      character(*) name
      integer   ivals(1)
      END FUNCTION nf_get_att_int

      FUNCTION NF_PUT_ATT_REAL(NCID, VARID, NAME, XTYPE, LEN, RVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_ATT_REAL' :: NF_PUT_ATT_REAL
!DEC$ ENDIF
      integer   nf_put_att_real
      integer   ncid, varid
      character(*) name
      integer   xtype, len
      real   rvals(1)
      END FUNCTION nf_put_att_real

      FUNCTION NF_GET_ATT_REAL(NCID, VARID, NAME, RVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_ATT_REAL' :: NF_GET_ATT_REAL
!DEC$ ENDIF
      integer   nf_get_att_real
      integer   ncid, varid
      character(*) name
      real  rvals(1)
      END FUNCTION nf_get_att_real

      FUNCTION NF_PUT_ATT_DOUBLE(NCID, VARID, NAME, XTYPE, LEN, DVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_ATT_DOUBLE' :: NF_PUT_ATT_DOUBLE
!DEC$ ENDIF
      integer   nf_put_att_double
      integer   ncid, varid
      character(*) name
      integer   xtype, len
      doubleprecision  dvals(1)
      END FUNCTION nf_put_att_double

      FUNCTION NF_GET_ATT_DOUBLE(NCID, VARID, NAME, DVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_ATT_DOUBLE' :: NF_GET_ATT_DOUBLE
!DEC$ ENDIF
      integer   nf_get_att_double
      integer   ncid, varid
      character(*) name
      doubleprecision  dvals(1)
      END FUNCTION nf_get_att_double

!
! general variable routines:
!

      FUNCTION NF_DEF_VAR(NCID, NAME, DATATYPE, NDIMS, DIMIDS, VARID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_DEF_VAR' :: NF_DEF_VAR
!DEC$ ENDIF
      integer   nf_def_var
      integer   ncid
      character(*) name
      integer   datatype, ndims
      integer   dimids(1)
      integer   varid
      END FUNCTION nf_def_var

      FUNCTION NF_INQ_VAR(NCID, VARID, NAME, DATATYPE, NDIMS, DIMIDS, NATTS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_VAR' :: NF_INQ_VAR
!DEC$ ENDIF
      integer   nf_inq_var
      integer   ncid, varid
      character(*) name
      integer   datatype, ndims
      integer   dimids(1)
      integer   natts
      END FUNCTION nf_inq_var

      FUNCTION NF_INQ_VARID(NCID, NAME, VARID)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_VARID' :: NF_INQ_VARID
!DEC$ ENDIF
      integer   nf_inq_varid
      integer   ncid
      character(*) name
      integer   varid
      END FUNCTION nf_inq_varid

      FUNCTION NF_INQ_VARNAME(NCID, VARID, NAME)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_VARNAME' :: NF_INQ_VARNAME
!DEC$ ENDIF
      integer   nf_inq_varname
      integer   ncid,  varid
      character(*) name
      END FUNCTION nf_inq_varname

      FUNCTION NF_INQ_VARTYPE(NCID, VARID, XTYPE)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_VARTYPE' :: NF_INQ_VARTYPE
!DEC$ ENDIF
      integer   nf_inq_vartype
      integer   ncid, varid, xtype
      END FUNCTION nf_inq_vartype

      FUNCTION NF_INQ_VARNDIMS(NCID, VARID, NDIMS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_VARNDIMS' :: NF_INQ_VARNDIMS
!DEC$ ENDIF
      integer   nf_inq_varndims
      integer   ncid,  varid, ndims
      END FUNCTION nf_inq_varndims

      FUNCTION NF_INQ_VARDIMID(NCID, VARID, DIMIDS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_VARDIMID' :: NF_INQ_VARDIMID
!DEC$ ENDIF
      integer   nf_inq_vardimid
      integer   ncid, varid
      integer   dimids(1)
      END FUNCTION nf_inq_vardimid

      FUNCTION NF_INQ_VARNATTS(NCID, VARID, NATTS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_INQ_VARNATTS' :: NF_INQ_VARNATTS
!DEC$ ENDIF
      integer   nf_inq_varnatts
      integer   ncid, varid, natts
      END FUNCTION nf_inq_varnatts

      FUNCTION NF_RENAME_VAR(NCID, VARID, NAME)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_RENAME_VAR' :: NF_RENAME_VAR
!DEC$ ENDIF
      integer   nf_rename_var
      integer   ncid,  varid
      character(*) name
      END FUNCTION nf_rename_var

      FUNCTION NF_COPY_VAR(NCID_IN, VARID, NCID_OUT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_COPY_VAR' :: NF_COPY_VAR
!DEC$ ENDIF
      integer   nf_copy_var
      integer   ncid_in, varid
      integer   ncid_out
      END FUNCTION nf_copy_var

!
! entire variable put/get routines:
!

      FUNCTION NF_PUT_VAR_TEXT(NCID, VARID, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR_TEXT' :: NF_PUT_VAR_TEXT
!DEC$ ENDIF
      integer   nf_put_var_text
      integer   ncid, varid
      character(*) text
      END FUNCTION nf_put_var_text

      FUNCTION NF_GET_VAR_TEXT(NCID, VARID, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR_TEXT' :: NF_GET_VAR_TEXT
!DEC$ ENDIF
      integer   nf_get_var_text
      integer   ncid, varid
      character(*) text
      END FUNCTION nf_get_var_text

      FUNCTION NF_PUT_VAR_INT1(NCID, VARID, I1VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR_INT1' :: NF_PUT_VAR_INT1
!DEC$ ENDIF
      integer   nf_put_var_int1
      integer   ncid, varid
      integer*1   i1vals(1)
      END FUNCTION nf_put_var_int1

      FUNCTION NF_GET_VAR_INT1(NCID, VARID, I1VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR_INT1' :: NF_GET_VAR_INT1
!DEC$ ENDIF
      integer   nf_get_var_int1
      integer   ncid,  varid
      integer*1   i1vals(1)
      END FUNCTION nf_get_var_int1

      FUNCTION NF_PUT_VAR_INT2(NCID, VARID, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR_INT2' :: NF_PUT_VAR_INT2
!DEC$ ENDIF
      integer   nf_put_var_int2
      integer   ncid, varid
      integer*2   i2vals(1)
      END FUNCTION nf_put_var_int2

      FUNCTION NF_GET_VAR_INT2(NCID, VARID, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR_INT2' :: NF_GET_VAR_INT2
!DEC$ ENDIF
      integer   nf_get_var_int2
      integer   ncid, varid
      integer*2   i2vals(1)
      END FUNCTION nf_get_var_int2

      FUNCTION NF_PUT_VAR_INT(NCID, VARID, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR_INT' :: NF_PUT_VAR_INT
!DEC$ ENDIF
      integer   nf_put_var_int
      integer   ncid, varid
      integer   ivals(1)
      END FUNCTION nf_put_var_int

      FUNCTION NF_GET_VAR_INT(NCID, VARID, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR_INT' :: NF_GET_VAR_INT
!DEC$ ENDIF
      integer   nf_get_var_int
      integer   ncid, varid
      integer   ivals(1)
      END FUNCTION nf_get_var_int

      FUNCTION NF_PUT_VAR_REAL(NCID, VARID, RVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR_REAL' :: NF_PUT_VAR_REAL
!DEC$ ENDIF
      integer   nf_put_var_real
      integer   ncid, varid
      real  rvals(1)
      END FUNCTION nf_put_var_real

      FUNCTION NF_GET_VAR_REAL(NCID, VARID, RVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR_REAL' :: NF_GET_VAR_REAL
!DEC$ ENDIF
      integer   nf_get_var_real
      integer   ncid, varid
      real  rvals(1)
      END FUNCTION nf_get_var_real

      FUNCTION NF_PUT_VAR_DOUBLE(NCID, VARID, DVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR_DOUBLE' :: NF_PUT_VAR_DOUBLE
!DEC$ ENDIF
      integer   nf_put_var_double
      integer   ncid,  varid
      doubleprecision     dvals(1)
      END FUNCTION nf_put_var_double

      FUNCTION NF_GET_VAR_DOUBLE(NCID, VARID, DVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR_DOUBLE' :: NF_GET_VAR_DOUBLE
!DEC$ ENDIF
      integer   nf_get_var_double
      integer   ncid, varid
      doubleprecision     dvals(1)
      END FUNCTION nf_get_var_double

!
! single variable put/get routines:
!

      FUNCTION NF_PUT_VAR1_TEXT(NCID, VARID, INDEX, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR1_TEXT' :: NF_PUT_VAR1_TEXT
!DEC$ ENDIF
      integer   nf_put_var1_text
      integer   ncid, varid
      integer   index(1)
      character*1  text
      END FUNCTION nf_put_var1_text

      FUNCTION NF_GET_VAR1_TEXT(NCID, VARID, INDEX, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR1_TEXT' :: NF_GET_VAR1_TEXT
!DEC$ ENDIF
      integer   nf_get_var1_text
      integer   ncid,  varid
      integer   index(1)
      character*1     text
      END FUNCTION nf_get_var1_text

      FUNCTION NF_PUT_VAR1_INT1(NCID, VARID, INDEX, I1VAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR1_INT1' :: NF_PUT_VAR1_INT1
!DEC$ ENDIF
      integer   nf_put_var1_int1
      integer   ncid, varid
      integer   index(1)
      integer*1       i1val
      END FUNCTION nf_put_var1_int1

      FUNCTION NF_GET_VAR1_INT1(NCID, VARID, INDEX, I1VAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR1_INT1' :: NF_GET_VAR1_INT1
!DEC$ ENDIF
      integer   nf_get_var1_int1
      integer   ncid, varid
      integer   index(1)
      integer*1      i1val
      END FUNCTION nf_get_var1_int1

      FUNCTION NF_PUT_VAR1_INT2(NCID, VARID, INDEX, I2VAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR1_INT2' :: NF_PUT_VAR1_INT2
!DEC$ ENDIF
      integer   nf_put_var1_int2
      integer   ncid, varid
      integer   index(1)
      integer*2   i2val
      END FUNCTION nf_put_var1_int2

      FUNCTION NF_GET_VAR1_INT2(NCID, VARID, INDEX, I2VAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR1_INT2' :: NF_GET_VAR1_INT2
!DEC$ ENDIF
      integer   nf_get_var1_int2
      integer   ncid,  varid
      integer   index(1)
      integer*2   i2val
      END FUNCTION nf_get_var1_int2

      FUNCTION NF_PUT_VAR1_INT(NCID, VARID, INDEX, IVAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR1_INT' :: NF_PUT_VAR1_INT
!DEC$ ENDIF
      integer   nf_put_var1_int
      integer   ncid, varid
      integer   index(1)
      integer   ival
      END FUNCTION nf_put_var1_int

      FUNCTION NF_GET_VAR1_INT(NCID, VARID, INDEX, IVAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR1_INT' :: NF_GET_VAR1_INT
!DEC$ ENDIF
      integer   nf_get_var1_int
      integer   ncid, varid
      integer   index(1)
      integer   ival
      END FUNCTION nf_get_var1_int

      FUNCTION NF_PUT_VAR1_REAL(NCID, VARID, INDEX, RVAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR1_REAL' :: NF_PUT_VAR1_REAL
!DEC$ ENDIF
      integer   nf_put_var1_real
      integer   ncid, varid
      integer   index(1)
      real  rval
      END FUNCTION nf_put_var1_real

      FUNCTION NF_GET_VAR1_REAL(NCID, VARID, INDEX, RVAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR1_REAL' :: NF_GET_VAR1_REAL
!DEC$ ENDIF
      integer   nf_get_var1_real
      integer   ncid, varid
      integer   index(1)
      real  rval
      END FUNCTION nf_get_var1_real

      FUNCTION NF_PUT_VAR1_DOUBLE(NCID, VARID, INDEX, DVAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VAR1_DOUBLE' :: NF_PUT_VAR1_DOUBLE
!DEC$ ENDIF
      integer   nf_put_var1_double
      integer   ncid, varid
      integer   index(1)
      doubleprecision     dval
      END FUNCTION nf_put_var1_double

      FUNCTION NF_GET_VAR1_DOUBLE(NCID, VARID, INDEX, DVAL)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VAR1_DOUBLE' :: NF_GET_VAR1_DOUBLE
!DEC$ ENDIF
      integer   nf_get_var1_double
      integer   ncid, varid
      integer   index(1)
      doubleprecision     dval
      END FUNCTION nf_get_var1_double

!
! variable array put/get routines:
!

      FUNCTION NF_PUT_VARA_TEXT(NCID, VARID, START, COUNT, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARA_TEXT' :: NF_PUT_VARA_TEXT
!DEC$ ENDIF
      integer   nf_put_vara_text
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      character(*) text
      END FUNCTION nf_put_vara_text

      FUNCTION NF_GET_VARA_TEXT(NCID, VARID, START,COUNT, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARA_TEXT' :: NF_GET_VARA_TEXT
!DEC$ ENDIF
      integer   nf_get_vara_text
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      character(*) text
      END FUNCTION nf_get_vara_text

      FUNCTION NF_PUT_VARA_INT1(NCID, VARID, START, COUNT, I1VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARA_INT1' :: NF_PUT_VARA_INT1
!DEC$ ENDIF
      integer   nf_put_vara_int1
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer*1       i1vals(1)
      END FUNCTION nf_put_vara_int1

      FUNCTION NF_GET_VARA_INT1(NCID, VARID, START,COUNT, I1VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARA_INT1' :: NF_GET_VARA_INT1
!DEC$ ENDIF
      integer   nf_get_vara_int1
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer*1       i1vals(1)
      END FUNCTION nf_get_vara_int1

      FUNCTION NF_PUT_VARA_INT2(NCID, VARID, START, COUNT, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARA_INT2' :: NF_PUT_VARA_INT2
!DEC$ ENDIF
      integer   nf_put_vara_int2
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer*2   i2vals(1)
      END FUNCTION nf_put_vara_int2

      FUNCTION NF_GET_VARA_INT2(NCID, VARID, START, COUNT, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARA_INT2' :: NF_GET_VARA_INT2
!DEC$ ENDIF
      integer   nf_get_vara_int2
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer*2   i2vals(1)
      END FUNCTION nf_get_vara_int2

      FUNCTION NF_PUT_VARA_INT(NCID, VARID, START, COUNT, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARA_INT' :: NF_PUT_VARA_INT
!DEC$ ENDIF
      integer   nf_put_vara_int
      integer   ncid,  varid
      integer   start(1)
      integer   count(1)
      integer   ivals(1)
      END FUNCTION nf_put_vara_int

      FUNCTION NF_GET_VARA_INT(NCID, VARID, START, COUNT, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARA_INT' :: NF_GET_VARA_INT
!DEC$ ENDIF
      integer   nf_get_vara_int
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   ivals(1)
      END FUNCTION nf_get_vara_int

      FUNCTION NF_PUT_VARA_REAL(NCID, VARID, START, COUNT, RVALS) 
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARA_REAL' :: NF_PUT_VARA_REAL
!DEC$ ENDIF
      integer   nf_put_vara_real
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      real  rvals(1)
      END FUNCTION nf_put_vara_real

      FUNCTION NF_GET_VARA_REAL(NCID, VARID, START, COUNT, RVALS)   
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARA_REAL' :: NF_GET_VARA_REAL
!DEC$ ENDIF
      integer   nf_get_vara_real
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      real  rvals(1)
      END FUNCTION nf_get_vara_real

      FUNCTION NF_PUT_VARA_DOUBLE(NCID, VARID, START, COUNT, DVALS) 
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARA_DOUBLE' :: NF_PUT_VARA_DOUBLE
!DEC$ ENDIF
      integer   nf_put_vara_double
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      doubleprecision     dvals(1)
      END FUNCTION nf_put_vara_double

      FUNCTION NF_GET_VARA_DOUBLE(NCID, VARID, START, COUNT, DVALS) 
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARA_DOUBLE' :: NF_GET_VARA_DOUBLE
!DEC$ ENDIF
      integer   nf_get_vara_double
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      doubleprecision     dvals(1)
      END FUNCTION nf_get_vara_double

!
! strided variable put/get routines:
!

      FUNCTION NF_PUT_VARS_TEXT(NCID, VARID, START, COUNT, STRIDE, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARS_TEXT' :: NF_PUT_VARS_TEXT
!DEC$ ENDIF
      integer   nf_put_vars_text
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      character(*) text
      END FUNCTION nf_put_vars_text

      FUNCTION NF_GET_VARS_TEXT(NCID, VARID, START, COUNT, STRIDE, TEXT)
      integer   nf_get_vars_text
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARS_TEXT' :: NF_GET_VARS_TEXT
!DEC$ ENDIF
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      character(*) text
      END FUNCTION nf_get_vars_text

      FUNCTION NF_PUT_VARS_INT1(NCID, VARID, START, COUNT, STRIDE, I1VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARS_INT1' :: NF_PUT_VARS_INT1
!DEC$ ENDIF
			integer   nf_put_vars_int1
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer*1       i1vals(1)
      END FUNCTION nf_put_vars_int1

      FUNCTION NF_GET_VARS_INT1(NCID, VARID, START, COUNT, STRIDE, I1VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARS_INT1' :: NF_GET_VARS_INT1
!DEC$ ENDIF
      integer   nf_get_vars_int1
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer*1        i1vals(1)
      END FUNCTION nf_get_vars_int1

      FUNCTION NF_PUT_VARS_INT2(NCID, VARID, START, COUNT, STRIDE, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARS_INT2' :: NF_PUT_VARS_INT2
!DEC$ ENDIF
      integer   nf_put_vars_int2
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer*2   i2vals(1)
      END FUNCTION nf_put_vars_int2

      FUNCTION NF_GET_VARS_INT2(NCID, VARID, START, COUNT, STRIDE, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARS_INT2' :: NF_GET_VARS_INT2
!DEC$ ENDIF
      integer   nf_get_vars_int2
      integer   ncid,  varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer*2   i2vals(1)
      END FUNCTION nf_get_vars_int2

      FUNCTION NF_PUT_VARS_INT(NCID, VARID, START, COUNT, STRIDE, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARS_INT' :: NF_PUT_VARS_INT
!DEC$ ENDIF
      integer   nf_put_vars_int
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   ivals(1)
      END FUNCTION nf_put_vars_int

      FUNCTION NF_GET_VARS_INT(NCID, VARID, START, COUNT, STRIDE, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARS_INT2' :: NF_GET_VARS_INT
!DEC$ ENDIF
      integer   nf_get_vars_int
      integer   ncid,  varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   ivals(1)
      END FUNCTION nf_get_vars_int

      FUNCTION NF_PUT_VARS_REAL(NCID, VARID, START, COUNT, STRIDE, RVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARS_REAL' :: NF_PUT_VARS_REAL
!DEC$ ENDIF
      integer   nf_put_vars_real
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      real  rvals(1)
      END FUNCTION nf_put_vars_real

      FUNCTION NF_GET_VARS_REAL(NCID, VARID, START, COUNT, STRIDE, RVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARS_REAL' :: NF_GET_VARS_REAL
!DEC$ ENDIF
      integer   nf_get_vars_real
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      real  rvals(1)
      END FUNCTION nf_get_vars_real

      FUNCTION NF_PUT_VARS_DOUBLE(NCID, VARID, START, COUNT, STRIDE, DVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARS_DOUBLE' :: NF_PUT_VARS_DOUBLE
!DEC$ ENDIF
      integer   nf_put_vars_double
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      doubleprecision     dvals(1)
      END FUNCTION nf_put_vars_double

      FUNCTION NF_GET_VARS_DOUBLE(NCID, VARID, START, COUNT, STRIDE, DVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARS_DOUBLE' :: NF_GET_VARS_DOUBLE
!DEC$ ENDIF
      integer   nf_get_vars_double
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      doubleprecision     dvals(1)
      END FUNCTION nf_get_vars_double

!
! mapped variable put/get routines:
!

      FUNCTION NF_PUT_VARM_TEXT(NCID, VARID, START, COUNT, STRIDE, IMAP, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARM_TEXT' :: NF_PUT_VARM_TEXT
!DEC$ ENDIF
      integer   nf_put_varm_text
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      character(*) text
      END FUNCTION nf_put_varm_text

      FUNCTION NF_GET_VARM_TEXT(NCID, VARID, START, COUNT, STRIDE, IMAP, TEXT)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARM_TEXT' :: NF_GET_VARM_TEXT
!DEC$ ENDIF
      integer   nf_get_varm_text
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      character(*) text
      END FUNCTION nf_get_varm_text

      FUNCTION NF_PUT_VARM_INT1(NCID, VARID, START, COUNT, STRIDE, IMAP, I1VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARM_INT1' :: NF_PUT_VARM_INT1
!DEC$ ENDIF
      integer   nf_put_varm_int1
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      integer*1    i1vals(1)
      END FUNCTION nf_put_varm_int1

      FUNCTION NF_GET_VARM_INT1(NCID, VARID, START, COUNT, STRIDE, IMAP, I1VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARM_INT1' :: NF_GET_VARM_INT1
!DEC$ ENDIF
      integer   nf_get_varm_int1
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      integer*1   i1vals(1)
      END FUNCTION nf_get_varm_int1

      FUNCTION NF_PUT_VARM_INT2(NCID, VARID, START, COUNT, STRIDE, IMAP, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARM_INT2' :: NF_PUT_VARM_INT2
!DEC$ ENDIF
      integer   nf_put_varm_int2
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      integer*2   i2vals(1)
      END FUNCTION nf_put_varm_int2

      FUNCTION NF_GET_VARM_INT2(NCID, VARID, START, COUNT, STRIDE, IMAP, I2VALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARM_INT2' :: NF_GET_VARM_INT2
!DEC$ ENDIF
      integer   nf_get_varm_int2
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      integer*2   i2vals(1)
      END FUNCTION nf_get_varm_int2

      FUNCTION NF_PUT_VARM_INT(NCID, VARID, START, COUNT, STRIDE, IMAP, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARM_INT' :: NF_PUT_VARM_INT
!DEC$ ENDIF
      integer   nf_put_varm_int
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      integer   ivals(1)
      END FUNCTION nf_put_varm_int

      FUNCTION NF_GET_VARM_INT(NCID, VARID, START, COUNT, STRIDE, IMAP, IVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARM_INT1' :: NF_GET_VARM_INT
!DEC$ ENDIF
      integer   nf_get_varm_int
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      integer   ivals(1)
      END FUNCTION nf_get_varm_int

      FUNCTION NF_PUT_VARM_REAL(NCID, VARID, START, COUNT, STRIDE, IMAP, RVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARM_REAL' :: NF_PUT_VARM_REAL
!DEC$ ENDIF
      integer   nf_put_varm_real
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      real  rvals(1)
      END FUNCTION nf_put_varm_real

      FUNCTION NF_GET_VARM_REAL(NCID, VARID, START, COUNT, STRIDE, IMAP, RVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARM_REAL' :: NF_GET_VARM_REAL
!DEC$ ENDIF
      integer   nf_get_varm_real
      integer   ncid,  varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      real  rvals(1)
      END FUNCTION nf_get_varm_real

      FUNCTION NF_PUT_VARM_DOUBLE(NCID, VARID, START, COUNT, STRIDE, IMAP, DVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_PUT_VARM_DOUBLE' :: NF_PUT_VARM_DOUBLE
!DEC$ ENDIF
      integer   nf_put_varm_double
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      doubleprecision     dvals(1)
      END FUNCTION nf_put_varm_double

      FUNCTION NF_GET_VARM_DOUBLE(NCID, VARID, START, COUNT, STRIDE, IMAP, DVALS)
!DEC$ IF( _DF_VERSION_ > 650 )
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:'NF_GET_VARM_DOUBLE' :: NF_GET_VARM_DOUBLE
!DEC$ ENDIF
      integer   nf_get_varm_double
      integer   ncid, varid
      integer   start(1)
      integer   count(1)
      integer   stride(1)
      integer   imap(1)
      doubleprecision     dvals(1)
      END FUNCTION nf_get_varm_double
      
END INTERFACE      

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 32)
      parameter(maxncdim = 100)
      parameter(maxncatt = 2000)
      parameter(maxncvar = 2000)
!     not enforced 
      parameter(maxncnam = 128)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)
