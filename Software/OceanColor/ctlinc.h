/*
   include file for control point routines
  
   Created by Gary Fu, 7/27/95
*/

typedef struct CtlInfoStruc {
         int    totpixl;   /* total pixel # */
         int    totline;   /* total line # */
         int    ncpp;      /* total # of pixel control points */
         int    ncpl;      /* total # of line control points */
         int    *cppix;    /* ptr to array of pixel control points (1-rel) */
         int    *cplin;    /* ptr to array of line control points (1-rel) */
         float  *cplat;    /* pointer to 2D lat array of control points */
         float  *cplon;    /* pointer to 2D lon array of control points */
         float  limits[4];   /* latmin, lonmin, latmax, lonmax */
} CtlInfoType;
