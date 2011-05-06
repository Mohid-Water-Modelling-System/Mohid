/* Compile with: 
cc readL2scan.c -g -c readL2scan.o -I$HDFROOT/include 
*/

/*
    Modification history:
    Programmer       Organization      Date      Description of change
    --------------   ------------    --------    ---------------------
    Joel Gales       Futuretech      01/31/00    Original Development
    Bryan Franz      GSC             03/01/00    Add reading of LAC 
                                                 pix start & subsamp to 
                                                 readL2meta
    Joel Gales       Futuretech      03/03/00    Change nflag in l2_str
                                                 from byte to int32
                                                 Fix allocation
                                                 problem with prodlist

    Joel Gales       Futuretech      03/14/00    Add "getL3units" routine

    Joel Gales       Futuretech      05/25/00    Fix isalpha (Linux) bug
    Joel Gales       Futuretech      06/14/00    Fix case where there are
                                                 no L2 products (flags only)
    Joel Gales       Futuretech      06/15/00    Fix units problem for flag
                                                 products
    Joel Gales       Futuretech      06/20/00    Add read support for FLOAT32
                                                 data products
*/


#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "readL2scan.h"

//extern void __stdcall FGEONAV (float a, float b, float *c);
//extern void __stdcall FGEONAV (float pos[2], float rm[8], float coef[6], float sun[3],
//							   int *nsta, int *ninc, int *npix, float *latitude[1284], float *longitude[1284]);

static int32 n_files_open=0;
static int32 sd_id_file[MAXNFILES];
static int32 sds_id_prod[MAXNFILES][100];
static int32 sds_id_ll[MAXNFILES][3];
static int32 sds_id_date[MAXNFILES][3];
static int32 sds_id_geonav[MAXNFILES][6];
static int32 sds_id_l2_flags[MAXNFILES];
static int32 sds_id_eng_qual[MAXNFILES];
static int32 sds_id_s_flags[MAXNFILES];
static int32 sds_id_nflag[MAXNFILES];
static int32 nsta[MAXNFILES];
static int32 ninc[MAXNFILES];
static int32 n_cntl_pnts;
static int32 prev_n_cntl_pnts=-1;
static int32 prodtype[100];

static float32  slope[100];
static float32  intercept[100];
static float32  geonav[6][9];

static int32 one = 1;

static int32 l2_flags_type;

static unsigned char *databuf[MAXNFILES];
static char *prodlist[MAXNFILES];
int16 *t_ranges[2*20];



int32 openL2(char *fname, char *plist, l2_prod *l2_str)
{
  
  int32 i;
  int32 sd_id;
  int32 sds_id;
  int32 len;
  int32 dims[8];
  int32 rank;
  int32 dtype;
  int32 nattrs;
  int32 n_l2flags;
  int32 attr_indx;
  int32 status;
  int32 tilt_start[2] = {0,0};
  int32 tilt_edges[2] = {20,2};
  int px,ln;
 

  int32 HDFfid;
  int32 vg_ref;
  int32 vgid;
  int32 tag;
  int32 ref;
  int32 listlen=0;
  int32 fileindex;
  static int32 prev_listlen=-1;
  

  char  buffer[2048];
  char  *cptr1;
  char  *cptr2;
  char  *numstr[] = {"01","02","03","04","05","06","07","08","09","10",
		     "11","12","13","14","15","16","17","18","19","20",
                     "21","22","23","24","25","26","27","28","29","30",
                     "31","32","33"};

  static int32 first=1;


  /* Copy filename and product list into L2 structure */
  /* ------------------------------------------------ */
  if (l2_str->nrec == 0) {
    strcpy(l2_str->filename, fname);
    n_files_open++;
    fileindex = n_files_open-1;
    l2_str->fileindex = fileindex;
  } else {
    fileindex = l2_str->fileindex;
  }


  if (first) {
    for (i=0; i<MAXNFILES; i++) prodlist[i]=NULL;
    for (i=0; i<MAXNFILES; i++) databuf[i]=NULL;
  }
  /* pina

  /* Generate prodlist */
  /* ----------------- */
  if (plist != 0x0) {
    prodlist[fileindex] = (char *) calloc(strlen(plist) + 1, sizeof(char));
    strcpy(prodlist[fileindex], plist);
  }
  else {
    HDFfid = Hopen(fname, DFACC_READ, 0);

    if (HDFfid == -1) {
      printf("\n%s not found\n", fname);
      exit(1);
    }

    sd_id = SDstart(l2_str->filename, DFACC_RDONLY);
    if (sd_id == -1) {
      printf("Error opening (SDstart) %s\n", l2_str->filename);
      exit(-1);
    }

    status = Vstart(HDFfid);
    vg_ref = Vfind(HDFfid, "Geophysical Data"); 
    vgid = Vattach(HDFfid, vg_ref, "r");
    
    for (i=0; i<Vntagrefs(vgid); i++) {
      status = Vgettagref(vgid, i, &tag, &ref);
      sds_id = SDselect(sd_id, SDreftoindex(sd_id, ref));
      if (sds_id == -1) {
	printf("Error accessing SDS (reference #: %d) in: %s .\n", 
	       ref, l2_str->filename);
	exit(-1);
      }
      status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);
      if (strcmp(buffer, "l2_flags") != 0) { 
	listlen += strlen(buffer) + 1;
      }
      SDendaccess(sds_id);
    }


    prodlist[fileindex] = (char *) calloc(listlen + 1, sizeof(char));

    for (i=0; i<Vntagrefs(vgid); i++) {
      status = Vgettagref(vgid, i, &tag, &ref);
      sds_id = SDselect(sd_id, SDreftoindex(sd_id, ref));
      status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);
      if (strcmp(buffer, "l2_flags") != 0) {
	if (i == 0) strcpy(prodlist[fileindex], buffer); 
	else strcat(prodlist[fileindex], buffer);
	strcat(prodlist[fileindex], ":");
      }
      else {
	l2_flags_type = dtype;	
      }
      SDendaccess(sds_id);
    }
    prodlist[fileindex][listlen-1] = 0;

    SDend(sd_id);
    Vdetach(vgid);
    Vend(HDFfid);

    Hclose(HDFfid);
  }


  /* Parse Product list */
  /* ------------------ */
  l2_str->nprod = 1;
  l2_str->prodname[0] = &prodlist[fileindex][0];
  len = strlen(prodlist[fileindex]);
  for (i=0; i<len; i++) {
    if (prodlist[fileindex][i] == ':') {
      l2_str->prodname[l2_str->nprod] = prodlist[fileindex] + i + 1; 
      l2_str->nprod++;
      prodlist[fileindex][i] = 0;
    }
  }
  /* No L2 products in L2 file (flags only) */
  if (strlen(prodlist[fileindex]) == 0) l2_str->nprod = 0;

    

  /* Start SD interface */
  /* ------------------ */
  sd_id = SDstart(l2_str->filename, DFACC_RDONLY);


  /* Get # of scans and # of pixels */
  /* ------------------------------ */
  SDreadattr(sd_id, SDfindattr(sd_id,"Number of Scan Lines"),
	     (VOIDP) &l2_str->nrec);
  SDreadattr(sd_id, SDfindattr(sd_id,"Pixels per Scan Line"),
	     (VOIDP) &l2_str->nsamp);


  /* Get start & end times, orbit number and data type */
  /* ------------------------------------------------- */
  SDreadattr(sd_id, SDfindattr(sd_id,"Start Year"),(VOIDP) &l2_str->syear);
  SDreadattr(sd_id, SDfindattr(sd_id,"Start Day"),(VOIDP) &l2_str->sday);
  SDreadattr(sd_id, SDfindattr(sd_id,"Start Millisec"),(VOIDP) &l2_str->smsec);
  SDreadattr(sd_id, SDfindattr(sd_id,"End Year"),(VOIDP) &l2_str->eyear);
  SDreadattr(sd_id, SDfindattr(sd_id,"End Day"),(VOIDP) &l2_str->eday);
  SDreadattr(sd_id, SDfindattr(sd_id,"End Millisec"),(VOIDP) &l2_str->emsec);
  SDreadattr(sd_id, SDfindattr(sd_id,"Orbit Number"),(VOIDP) &l2_str->orbit);
  SDreadattr(sd_id, SDfindattr(sd_id,"Data Type"),(VOIDP) &l2_str->dtype);


  /* Allocate geoloc (lon/lat) arrays */
  /* -------------------------------- */
  l2_str->geoloc = (float32 *) calloc(2 * l2_str->nsamp, sizeof(float32));
  l2_str->latitude  = l2_str->geoloc;
  l2_str->longitude = l2_str->geoloc + l2_str->nsamp;


  /* Get longitude, latitude, & date sds ids */
  /* --------------------------------------- */
  sds_id_ll[fileindex][0] = SDselect(sd_id, SDnametoindex(sd_id,"longitude"));
  sds_id_ll[fileindex][1] = SDselect(sd_id, SDnametoindex(sd_id,"latitude"));
  sds_id_ll[fileindex][2] = -1;
  sds_id_date[fileindex][0]  = SDselect(sd_id, SDnametoindex(sd_id,"year"));
  sds_id_date[fileindex][1]  = SDselect(sd_id, SDnametoindex(sd_id,"day"));
  sds_id_date[fileindex][2]  = SDselect(sd_id, SDnametoindex(sd_id,"msec"));

  
   


  l2_str->geointerp = 0;


  /* Test if full-size lon/lat SDS */
  /* ----------------------------- */
  status = SDgetinfo(sds_id_ll[fileindex][0], buffer, &rank, dims, &dtype, &nattrs);
  if (dims[0] != l2_str->nrec || dims[1] != l2_str->nsamp || 
      sds_id_ll[fileindex][0] == -1) {

    l2_str->geointerp = 1;

    /* Test for geonav arrays (SeaWIFS) */
    /* -------------------------------- */
    sds_id_geonav[fileindex][0] = SDselect(sd_id, SDnametoindex(sd_id,"orb_vec"));
    sds_id_geonav[fileindex][1] = SDselect(sd_id, SDnametoindex(sd_id,"sen_mat"));
    sds_id_geonav[fileindex][2] = SDselect(sd_id, SDnametoindex(sd_id,"scan_ell"));
    sds_id_geonav[fileindex][3] = SDselect(sd_id, SDnametoindex(sd_id,"sun_ref"));
    sds_id_geonav[fileindex][4] = SDselect(sd_id, SDnametoindex(sd_id,"l_vert"));
    sds_id_geonav[fileindex][5] = SDselect(sd_id, SDnametoindex(sd_id,"att_ang"));

    nsta[fileindex] = -1;
    ninc[fileindex] = -1;
    SDreadattr(sd_id, SDfindattr(sd_id,"LAC Pixel Start Number"),
	       (VOIDP) &nsta[fileindex]);
    SDreadattr(sd_id, SDfindattr(sd_id,"LAC Pixel Subsampling"),
	       (VOIDP) &ninc[fileindex]);

    if (sds_id_geonav[fileindex][0] != -1 && sds_id_geonav[fileindex][1] != -1 &&
	sds_id_geonav[fileindex][2] != -1 && sds_id_geonav[fileindex][3] != -1 &&
	nsta[fileindex] != -1 && ninc[fileindex] != -1) {

      l2_str->geointerp = 2;

      for (i=0; i<6; i++) l2_str->geonav[i] = geonav[i];
    }
    else {

      /* Get # of control points */
      /* ----------------------- */
      SDreadattr(sd_id, SDfindattr(sd_id,"Number of Pixel Control Points"),
		 (VOIDP) &n_cntl_pnts);


      /* Check that all L2 files have same number of control points */
      /* ---------------------------------------------------------- */
      if (prev_n_cntl_pnts != -1 && prev_n_cntl_pnts != n_cntl_pnts) {
	printf("L2 file #:%4ld has %ld control points.\n", fileindex, prev_n_cntl_pnts);
	printf("L2 file #:%4ld has %ld control points.\n", fileindex+1, n_cntl_pnts);
	printf("These must be identical.\n");
	exit(-1);
      }
      prev_n_cntl_pnts = n_cntl_pnts;


      /* Allocate arrays needed for lon/lat interpolation */
      /* ------------------------------------------------ */
      l2_str->lon_cntl = (float32 *) calloc(n_cntl_pnts, sizeof(float32));
      l2_str->lat_cntl = (float32 *) calloc(n_cntl_pnts, sizeof(float32));
      l2_str->cntl_pnts = (float32 *) calloc(n_cntl_pnts, sizeof(float32));
      l2_str->spline_arr = (float32 *) calloc(n_cntl_pnts, sizeof(float32));

      /* Get control point sds id */
      /* ------------------------ */
      sds_id_ll[fileindex][2] = SDselect(sd_id, SDnametoindex(sd_id,"cntl_pt_cols"));
    }
  }



  /* Store data products info in L2 structure */
  /* ---------------------------------------- */
  for (i=0; i<l2_str->nprod; i++) {
    sds_id = SDselect(sd_id, SDnametoindex(sd_id,l2_str->prodname[i]));

    if (sds_id != -1) {
      status = SDgetinfo(sds_id, buffer, &rank, dims, &prodtype[i], &nattrs);

      /* Read scaling slope and intercept */
      /* -------------------------------- */
      SDreadattr(sds_id,SDfindattr(sds_id,"slope"), (VOIDP) &slope[i]);
      SDreadattr(sds_id,SDfindattr(sds_id,"intercept"), (VOIDP) &intercept[i]);

      sds_id_prod[fileindex][i] = sds_id;
    }
    else {
      printf("Data Product: \"%s\" not found.\n", l2_str->prodname[i]);
      exit(1);
    }
  }

  /* Allocate databuf data array */
  /* --------------------------- */
  databuf[fileindex] = (unsigned char *) calloc(l2_str->nsamp, 8);


  /* Allocate L2 data array */
  /* ---------------------- */
  l2_str->l2_data = (float32 *) calloc(l2_str->nprod * l2_str->nsamp, 
				       sizeof(float32));


  /* Read tilt data (if applicable) */
  /* ------------------------------ */
  if (SDselect(sd_id, SDnametoindex(sd_id,"ntilts")) != -1) {
    sds_id = SDselect(sd_id, SDnametoindex(sd_id,"ntilts"));
    status = SDreaddata(sds_id, tilt_start, NULL, &one, 
			(VOIDP) &l2_str->ntilts);
    sds_id = SDselect(sd_id, SDnametoindex(sd_id,"tilt_flags"));
    status = SDreaddata(sds_id, tilt_start, NULL, tilt_edges, 
			(VOIDP) l2_str->tilt_flags);
    sds_id = SDselect(sd_id, SDnametoindex(sd_id,"tilt_ranges"));
    status = SDreaddata(sds_id, tilt_start, NULL, tilt_edges, 
			(VOIDP) t_ranges);

    for (i=0; i<l2_str->ntilts; i++) {
      l2_str->tilt_ranges[0][i] = t_ranges[i*2];
      l2_str->tilt_ranges[1][i] = t_ranges[i*2+1];
    }
  }

  /* Open L2 flags SDS */
  /* ----------------- */
  sds_id_l2_flags[fileindex] = SDselect(sd_id, SDnametoindex(sd_id,"l2_flags"));
  if (sds_id_l2_flags[fileindex] != -1) {
    l2_str->l2_flags = (int32 *) calloc(l2_str->nsamp, sizeof(int32));

    /* Read L2 flagnames */
    /* ----------------- */
    n_l2flags = 0;
    listlen = 0;
    while(1) {
      sprintf(buffer, "f%s_name", numstr[n_l2flags]);
      attr_indx = SDfindattr(sds_id_l2_flags[fileindex],buffer);
      if (attr_indx != -1) {
	SDreadattr(sds_id_l2_flags[fileindex], attr_indx, (VOIDP) buffer);
	if (strcmp(buffer, "SPARE") == 0) break;
	listlen += strlen(buffer) + 1;
	n_l2flags++;
      } else break;
    }
    l2_str->flagnames = (char *) calloc(listlen, sizeof(char));
    for (i=0; i<n_l2flags; i++) {
      sprintf(buffer, "f%s_name", numstr[i]);
      attr_indx = SDfindattr(sds_id_l2_flags[fileindex],buffer);
      SDreadattr(sds_id_l2_flags[fileindex], attr_indx, (VOIDP) buffer);
      strcat(l2_str->flagnames, buffer);
      if (i < n_l2flags-1) strcat(l2_str->flagnames, ",");
    }
  } else {
    l2_str->l2_flags = 0x0;
  }


  /* Build flag mask */
  /* --------------- */
  SDreadattr(sd_id,SDfindattr(sd_id,"Processing Control"), (VOIDP) buffer);
  cptr1 = strstr(buffer, "|MSKFLG");
  if (cptr1 != 0x0) {
    strcpy(buffer, cptr1+7);
    for (i=0; i<strlen(buffer); i++)
      if (isalpha((int) buffer[i]) != 0) break;
    strcpy(buffer, &buffer[i]);
    cptr1 = strstr(buffer, "|");
    *cptr1 = 0;

    l2_str->flagmask = 0;
    strcat(l2_str->flagnames, ",");
    cptr1 = l2_str->flagnames;
    for (i=0; i<16; i++) {
      cptr2 = strstr(cptr1, ",");
      *cptr2 = 0;
      if (strstr(buffer, cptr1) != 0x0)	l2_str->flagmask += (one << i);
      *cptr2 = ',';
      cptr1 = cptr2 + 1;
    }
    *cptr2 = 0;
  }


  /* Open eng_qual SDS */
  /* ----------------- */
  sds_id_eng_qual[fileindex] = SDselect(sd_id, SDnametoindex(sd_id,"eng_qual"));


  /* Open s_flags SDS */
  /* ---------------- */
  sds_id_s_flags[fileindex] = SDselect(sd_id, SDnametoindex(sd_id,"s_flags"));


  /* Open nflag SDS */
  /* -------------- */
  sds_id_nflag[fileindex] = SDselect(sd_id, SDnametoindex(sd_id,"nflag"));


  sd_id_file[fileindex] = sd_id;
  /*  printf("sd_id(O): %ld %ld\n", fileindex, sd_id);*/

  first = 1;


  return 0;
}


int32 reopenL2(int32 fileindex, l2_prod *l2_str)
{
  int32 i;
  int32 sd_id;

  sd_id_file[fileindex] = SDstart(l2_str->filename, DFACC_RDONLY);

  sd_id = sd_id_file[fileindex];
  /*  printf("sd_id(R): %ld %ld\n", fileindex, sd_id);*/

  sds_id_ll[fileindex][0] = SDselect(sd_id, SDnametoindex(sd_id,"longitude"));
  sds_id_ll[fileindex][1] = SDselect(sd_id, SDnametoindex(sd_id,"latitude"));

  sds_id_date[fileindex][0]  = SDselect(sd_id, SDnametoindex(sd_id,"year"));
  sds_id_date[fileindex][1]  = SDselect(sd_id, SDnametoindex(sd_id,"day"));
  sds_id_date[fileindex][2]  = SDselect(sd_id, SDnametoindex(sd_id,"msec"));

  if (l2_str->geointerp == 1) {
    sds_id_ll[fileindex][2] = SDselect(sd_id, SDnametoindex(sd_id,"cntl_pt_cols"));
  }

  if (l2_str->geointerp == 2) {
    sds_id_geonav[fileindex][0] = SDselect(sd_id, SDnametoindex(sd_id,"orb_vec"));
    sds_id_geonav[fileindex][1] = SDselect(sd_id, SDnametoindex(sd_id,"sen_mat"));
    sds_id_geonav[fileindex][2] = SDselect(sd_id, SDnametoindex(sd_id,"scan_ell"));
    sds_id_geonav[fileindex][3] = SDselect(sd_id, SDnametoindex(sd_id,"sun_ref"));
    sds_id_geonav[fileindex][4] = SDselect(sd_id, SDnametoindex(sd_id,"l_vert"));
    sds_id_geonav[fileindex][5] = SDselect(sd_id, SDnametoindex(sd_id,"att_ang"));
  }

  for (i=0; i<l2_str->nprod; i++) {
    sds_id_prod[fileindex][i] = SDselect(sd_id, SDnametoindex(sd_id,l2_str->prodname[i]));
  }

  sds_id_l2_flags[fileindex] = SDselect(sd_id, SDnametoindex(sd_id,"l2_flags"));

  sds_id_eng_qual[fileindex] = SDselect(sd_id, SDnametoindex(sd_id,"eng_qual"));
  sds_id_s_flags[fileindex] = SDselect(sd_id, SDnametoindex(sd_id,"s_flags"));
  sds_id_nflag[fileindex] = SDselect(sd_id, SDnametoindex(sd_id,"nflag"));
}


int32 readL2(l2_prod *l2_str, int32 ifile, int32 recnum, int32 iprod)
{
  int32 i;
  int32 start[3] = {0,0,0};
  int32 edges[3];
  int32 ipix;
  int16 tempi16;
  int32 status;

  float32 slp;
  float32 itp;

  float32 tempf32;
  unsigned char tempi8;

  int32 flag_edges[2] = {1,4};
  int32 nflag_edges[2] = {1,8};
  int16 zero = 0;


  start[0] = recnum;
  start[1] = 0;
  edges[0] = 1;
  edges[1] = l2_str->nsamp;


  /* Main product loop */
  /* ----------------- */
  for (i=0; i<l2_str->nprod; i++) {

    if ((iprod != -1) && (i != iprod)) continue;

    slp = slope[i];
    itp = intercept[i];

    /* Read into data buffer */
    /* --------------------- */
    status = SDreaddata(sds_id_prod[ifile][i], start, NULL, edges, (VOIDP) databuf[ifile]);
    if (status != 0) {
      printf("Read Error: %ld %ld\n", ifile, i);
      exit(-1);
    }


    /* Convert to proper data type (unscale) */
    /* ------------------------------------- */
    switch (prodtype[i]) {

    case DFNT_UINT8:
      for (ipix=0; ipix<l2_str->nsamp; ipix++) {
	memcpy(&tempi8, &databuf[ifile][ipix], 1);
	tempf32 = tempi8 * slp + itp;
	memcpy(&l2_str->l2_data[i*l2_str->nsamp + ipix], &tempf32, 
	       sizeof(float32));
      }
      break;

    case DFNT_INT16:
      for (ipix=0; ipix<l2_str->nsamp; ipix++) {
	memcpy(&tempi16, &databuf[ifile][2*ipix], 2);
	tempf32 = tempi16 * slp + itp;
	memcpy(&l2_str->l2_data[i*l2_str->nsamp + ipix], &tempf32, 
	       sizeof(float32));
      }
      break;

    case DFNT_FLOAT32:
      for (ipix=0; ipix<l2_str->nsamp; ipix++) {

	memcpy(&tempf32, &databuf[ifile][4*ipix], 4);
	tempf32 = tempf32 * slp + itp;
	memcpy(&l2_str->l2_data[i*l2_str->nsamp + ipix], &tempf32, 
	       sizeof(float32));
      }
      break;

    }; /* end switch */
    
  } /* product loop */

 
  return 0;
}

int32 readInicL2(l2_prod *l2_str, int32 ifile, int32 recnum)
{
  int32 i;
  int32 start[3] = {0,0,0};
  int32 edges[3];
  int32 ipix;
  int16 tempi16;
  int32 status;

  float32 slp;
  float32 itp;

  float32 tempf32;
  unsigned char tempi8;

  int32 flag_edges[2] = {1,4};
  int32 nflag_edges[2] = {1,8};
  int16 zero = 0;


  start[0] = recnum;
  start[1] = 0;
  edges[0] = 1;
  edges[1] = l2_str->nsamp;


  /* Read L2 flags */
  /* ------------- */
  if (sds_id_l2_flags[ifile] != -1) { 
    status = SDreaddata(sds_id_l2_flags[ifile], start, NULL, edges, 
			(VOIDP) l2_str->l2_flags);

    /* If INT16 array then convert to INT32 */
    /* ------------------------------------ */
    if (l2_flags_type == DFNT_INT16) {
      for (i=edges[1]-1; i>=0; i--)
	memcpy((int16 *) l2_str->l2_flags + (2*i+1), 
	       (int16 *) l2_str->l2_flags + i, 2);
      for (i=0; i<edges[1]; i++)
	memcpy((int16 *) l2_str->l2_flags + (2*i), &zero, 2);
    }
  }

  /* Read eng_qual */
  /* ------------- */
  if (sds_id_eng_qual[ifile] != -1) 
    status = SDreaddata(sds_id_eng_qual[ifile], start, NULL, flag_edges, 
			(VOIDP) l2_str->eng_qual);

  /* Read s_flags */
  /* ------------ */
  if (sds_id_s_flags[ifile] != -1) 
    status = SDreaddata(sds_id_s_flags[ifile], start, NULL, flag_edges, 
			(VOIDP) l2_str->s_flags);

  /* Read nflag */
  /* ---------- */
  if (sds_id_nflag[ifile] != -1) 
    status = SDreaddata(sds_id_nflag[ifile], start, NULL, nflag_edges, 
			(VOIDP) l2_str->nflag);

  /* Read date fields */
  /* ---------------- */
  if (sds_id_date[ifile][0] != -1) 
    status = SDreaddata(sds_id_date[ifile][0], &recnum, NULL, &one, 
			(VOIDP) &l2_str->year);
  if (sds_id_date[ifile][1] != -1) 
    status = SDreaddata(sds_id_date[ifile][1], &recnum, NULL, &one, 
			(VOIDP) &l2_str->day);
  if (sds_id_date[ifile][2] != -1) 
    status = SDreaddata(sds_id_date[ifile][2], &recnum, NULL, &one, 
			(VOIDP) &l2_str->msec);


  /* Read lon/lat fields (Note: start & edges are changed) */
  /* ----------------------------------------------------- */
  readlonlat(l2_str, ifile, start, edges);


  /* Check whether lon/lat values are within range */
  /* --------------------------------------------- */
  for(i=0; i<l2_str->nsamp; i++){
    if ((l2_str->longitude[i] > 180 || l2_str->longitude[i] < -180) &&
	(l2_str->l2_flags[i] & 33554432 == 0)) {
      printf("Scheme: %d\n", l2_str->geointerp);
      printf("Pixel Longitude %d out of range (%f) for scan %d in %s.\n", 
	     i, l2_str->longitude[i], recnum, l2_str->filename);
      exit(-1);
    }


    if ((l2_str->latitude[i] > 180 || l2_str->latitude[i] < -180) &&
	(l2_str->l2_flags[i] & 33554432 == 0)) {
      printf("Scheme: %d\n", l2_str->geointerp);
      printf("Pixel Latitude %d out of range (%f) for scan %d in %s.\n", 
	     i, l2_str->latitude[i], recnum, l2_str->filename);
      exit(-1);
    }
  }

  return 0;
}


int32 readlonlat(l2_prod *l2_str, int32 ifile, int32 *start, int32 *edges)
{
  int32 i;
  int32 status;
  int32 tempi32;
  int32 geo_edge[6] = {3,3,6,3,3,3};
  
   float lat,lon;

   int    tot_pixl ;  
   int    tot_line  ; 
   int    lstart_pix ;
   int    lpix_sample;
   int    ntilts     ;
   short* tilt_flags ;

  float32 tempf32;
  float32 delta;
  long first;

  float geovec[2];
  float no[2];
  float up[2];
  float ea[2];
  float rmtq[2];
  float sina[1284];
  float cosa[1284];
  float *auxSmat;
  float *auxOrbVect;
  float *auxScanEll;

  float auxLat[1284];
  float auxLon[1284];

  int j,in,px, ln, ipair;

  char *sds_name;
  int32 *rank;
  int32 dimsizes[50];
  int32 *data_type;
  int32 *num_attrs;

  float32 pi,radeg,re,rem,f,omf2,omegae ;
  float32 sinc,elev,sinl,cosl,h;
  float32 a,b,c,r;
  float32 q,Qx,Qy,Qz;
  float32 tmp,temp,uxy,upxy,sv,sn,se,s;
  float32 sunv,sunn,sune;
       
  first=1;
  switch (l2_str->geointerp) {

  case 0:

    /* Read full-size lon/lat fields */
    /* ----------------------------- */
    status = SDreaddata(sds_id_ll[ifile][0], start, NULL, edges, 
			(VOIDP) l2_str->longitude);
    status = SDreaddata(sds_id_ll[ifile][1], start, NULL, edges, 
			(VOIDP) l2_str->latitude);
    break;


  case 1:

    /* Read subsampled lon/lat fields */
    /* ------------------------------ */
    edges[1] = n_cntl_pnts;
    status = SDreaddata(sds_id_ll[ifile][0], start, NULL, edges, 
			(VOIDP)	l2_str->lon_cntl);

    status = SDreaddata(sds_id_ll[ifile][1], start, NULL, edges, 
			(VOIDP) l2_str->lat_cntl);

    /* Read control points array if first time through */
    /* ----------------------------------------------- */
    status = SDreaddata(sds_id_ll[ifile][2], &start[1], NULL, &edges[1], 
			(VOIDP) databuf[ifile]);

    /* Convert cntl pnts from I32 to F32 */
    /* --------------------------------- */
    for (i=0; i<n_cntl_pnts; i++) {
      memcpy(&tempi32, &databuf[ifile][4*i], sizeof(int32));
      tempf32 = (float32) tempi32;
      memcpy(&l2_str->cntl_pnts[i], &tempf32, sizeof(float32));
    }


    /* Remove any dateline discontinuity in the longitudes */
    /* --------------------------------------------------- */
    for(i=1; i<n_cntl_pnts; i++){
      delta = l2_str->lon_cntl[i] - l2_str->lon_cntl[i-1];
      if (delta < -180) l2_str->lon_cntl[i] += 360; 
      else if (delta > 180) l2_str->lon_cntl[i] -= 360;
    }


    /* Interpolate Latitude */
    /* -------------------- */
    spline(l2_str->cntl_pnts,l2_str->lat_cntl,n_cntl_pnts,
	   1e30,1e30,l2_str->spline_arr);
    for(i=0; i<l2_str->nsamp; i++){
      splint(l2_str->cntl_pnts,l2_str->lat_cntl,l2_str->spline_arr,n_cntl_pnts,
	     i+1.0,&l2_str->latitude[i]);
    }


    /* Interpolate Longitude */
    /* --------------------- */
    spline(l2_str->cntl_pnts,l2_str->lon_cntl,n_cntl_pnts,
	   1e30,1e30,l2_str->spline_arr);

    for(i=0; i<l2_str->nsamp; i++){
      splint(l2_str->cntl_pnts,l2_str->lon_cntl,l2_str->spline_arr,n_cntl_pnts,
	     i+1.0,&l2_str->longitude[i]);

      /* Put the longitudes back in the [-180,180] range */
      /* ----------------------------------------------- */
      while(l2_str->longitude[i] >  180) 
	l2_str->longitude[i] -= 360;
      while(l2_str->longitude[i] < -180) 
	l2_str->longitude[i] += 360;
    }

    break;

  case 2:

    edges[2] = 3;
    for (i=0; i<4; i++) {
      edges[1] = geo_edge[i];
      if (sds_id_geonav[ifile][i] != -1) {
	status = SDreaddata(sds_id_geonav[ifile][i], start, NULL, edges, 
			    (VOIDP) geonav[i]);
      }
    }


	FGEONAV(geonav[0], geonav[1], geonav[2], geonav[3],(int *) &nsta[ifile], 
		    (int *) &ninc[ifile], (int *) &l2_str->nsamp, &auxLat, &auxLon );

	for (i=0; i<l2_str->nsamp; i++) {
		l2_str->longitude[i]=auxLon[i];
		l2_str->latitude[i]= auxLat[i];
	}

		break;
  }


//  tot_pixl    = l2_str->nsamp;
//  tot_line    = l2_str->nrec;
//  lstart_pix  = ninc[0];
//  lpix_sample = nsta[0];
//  ntilts      = l2_str->ntilts;
//  tilt_flags  = l2_str->tilt_flags;
  


//   start[0]=0;
//   start[1]=0;
//   start[2]=0;

//   edges[0]=l2_str->nrec;
//   edges[1]=3;
//   edges[2]=0;

//   if ( (auxOrbVect = (float *)malloc(edges[0]*edges[1]*sizeof(float))) == NULL ) {
//           	printf("-E- %s: Error allocating memory for L2 file index.\n");
//           	exit(1);
//    	    }
//   status = SDreaddata(sds_id_geonav[0][0], start, NULL, edges, 
//			    (VOIDP) auxOrbVect);

 //  start[0]=0;
 //  start[1]=0;
 //  start[2]=0;
//
 //  edges[0]=l2_str->nrec;
 //  edges[1]=3;
 //    edges[2]=3;

 //     if ( (auxSmat = (float *)malloc(edges[0]*edges[1]*edges[2]*sizeof(float))) == NULL ) {
 //          	printf("-E- %s: Error allocating memory for L2 file index.\n");
 //          	exit(1);
 //   	    }
 //  status = SDreaddata(sds_id_geonav[0][1], start, NULL, edges, 
//			    (VOIDP) auxSmat);
//
//
//   start[0]=0;
 //  start[1]=0;
 //  start[2]=0;
//
 //  edges[0]=l2_str->nrec;
 //  edges[1]=6;
 //  edges[2]=0;


 //  if ( (auxScanEll = (float *)malloc(edges[0]*edges[1]*sizeof(float))) == NULL ) {
 //          	printf("-E- %s: Error allocating memory for L2 file index.\n");
 //          	exit(1);
//    	    }
//   status = SDreaddata(sds_id_geonav[0][0], start, NULL, edges, 
//			    (VOIDP) auxScanEll);

//  lat= -12.2;
//  lon = 38.6;
//  px=-1;
//  ln=-1;
//  ipair =1;

 // IDLswfll2tv ( l2_str->nsamp, l2_str->nrec,
   //         ninc[0], nsta[0], l2_str->ntilts, l2_str->tilt_flags,
     //       t_ranges, auxOrbVect, auxSmat, auxScanEll,ipair);
//IDLswfll2tv (lat, lon, px, ln, l2_str->nsamp, l2_str->nrec,
  //          ninc[0], nsta[0], l2_str->ntilts, l2_str->tilt_flags,
    //        t_ranges, auxOrbVect, auxSmat, auxScanEll,ipair);
 
 //IDLswfll2tv', navdata.tot_pixel, navdata.tot_line, $
   //      navdata.start_pix, navdata.pix_sample, navdata.tilt_str.ntilts, $
     //    navdata.tilt_str.tilt_flags, navdata.tilt_str.tilt_ranges, $
       //  navdata.pos, navdata.smat, navdata.coef, ipair, latlon, pxln)


  return 0;
}


/*griflet: translated from fgeonav.f existing in the ConvertToHdf5\SeaDas folder of the source-safe project. */
int32 FGEONAV (float pos[3], float rm[8], float coef[6], float sun[3],
			   int *nsta, int *ninc, int *npix, float latitude[1284], float longitude[1284])
{

	float a,b,c,r,q,Qx,Qy,Qz,tmp;
	int i,j;
	int in;

    float geovec[3];
	float no[3];
	float up[3];
	float ea[3];
	float rmtq[3];
    double pi,radeg,re,rem,f,omf2,omegae;
    double sinc,elev,sinl,cosl,h;
	int first;
	
	double sina[1285];
	double cosa[1285];

	first = 1;

	sinc = 0.0015911e0;
    ea[0] = 0.0;
	ea[1] = 0.0;
	ea[2] = 0.0;

	if (first == 1)
	{
		first = 0;

		/*cdata*/
		pi = acos(-1.0e0);
		re = 6378.137e0;
		rem = 6371.e0;
		radeg = 180.e0/pi;
		f = 1.e0/298.257e0;
		omf2 = pow(1.e0-f,2);
	    omegae = 7.29211585494e-5;
		/*end of cdata*/

		elev = sinc*1.2;
		sinl = sin(elev);
		cosl = cos(elev);
		for (i=0;i<1285;i++)
		{
		    sina[i] = sin((i+1-643)*sinc)*cosl;
			cosa[i] = cos((i+1-643)*sinc)*cosl;
		}
	}

	h = (rm[1,0]*pos[1]+rm[1,1]*pos[1]+rm[1,2]*pos[2]/omf2)*2.e0;

	for (i=0; i<npix; i++)
	{
	  in = (*ninc)*i + (*nsta);
	  a = coef[0]*cosa[in]*cosa[in] + coef[1]*cosa[in]*sina[in]*coef[2]*sina[in]*sina[in];
	  b = coef[3]*cosa[in] + coef[4]*sina[in];
	  c = coef[5];
	  r = b * b-4.e0 * c*a;

	  if (r < 0.)
	  {

	    latitude[i] = (float)999.;
		longitude[i] = (float)999.;

	  } else {

	    q = (-b-sqrt(r))/(2.e0*a);
	    q = q*(1.e0 + sinl*h/sqrt(r));
	    
	    Qx = q*cosa[in];
	    Qy = q*sinl;
	    Qz = q*sina[in];

	    for (j=0;j<3;j++)
		{
	      rmtq[j] = Qx*rm[0,j] + Qy*rm[1,j] + Qz*rm[2,j];
	      geovec[j] = rmtq[j] + pos[j];
		}

		tmp = sqrt(geovec[0]*geovec[0]+geovec[1]*geovec[1])*omf2;
	    latitude[i] = radeg*atan2(geovec[2],tmp);
	    longitude[i] = radeg*atan2(geovec[1],geovec[0]);
	  }
	}

	return 0;

}


int32 closeL2(l2_prod *l2_str, int32 ifile)
{
  int32 i;
  int32 status;

  for (i=0; i<l2_str->nprod; i++) {
    if (sds_id_prod[ifile][i] != -1) {
      status = SDendaccess(sds_id_prod[ifile][i]);
      if (status != 0) {
	printf("Error ending access to product sds: %ld for file: %ld\n", i, ifile);
	exit(-1);
      }
    }
  }

  for (i=0; i<3; i++) {
    if (sds_id_ll[ifile][i] != -1) {
      status = SDendaccess(sds_id_ll[ifile][i]);
      if (status != 0) {
	printf("Error ending access to ll sds: %ld for file: %ld\n", i, ifile);
	exit(-1);
      }
    }
  }

  for (i=0; i<3; i++) {
    if (sds_id_date[ifile][i] != -1) {
      status = SDendaccess(sds_id_date[ifile][i]);
      if (status != 0) {
	printf("Error ending access to date sds: %ld for file: %ld\n", i, ifile);
	exit(-1);
      }
    }
  }


  if (l2_str->geointerp == 2) {
    for (i=0; i<6; i++) {
      if (sds_id_geonav[ifile][i] != -1)
	status = SDendaccess(sds_id_geonav[ifile][i]);
      if (status != 0) {
	printf("Error ending access to geonav sds: %ld for file: %ld\n", i, ifile);
	exit(-1);
      }
    }
  }

  if (sds_id_l2_flags[ifile] != -1) {
    status = SDendaccess(sds_id_l2_flags[ifile]);
    if (status != 0) {
      printf("Error ending access to l2_flags sds for file: %ld\n", ifile);
      exit(-1);
    }
  }


  if (sds_id_eng_qual[ifile] != -1) {
    status = SDendaccess(sds_id_eng_qual[ifile]);
    if (status != 0) {
      printf("Error ending access to eng_qual sds for file: %ld\n", ifile);
      exit(-1);
    }
  }


  if (sds_id_s_flags[ifile]  != -1) {
    status = SDendaccess(sds_id_s_flags[ifile]);
    if (status != 0) {
      printf("Error ending access to s_flags sds for file: %ld\n", ifile);
      exit(-1);
    }
  }

  if (sds_id_nflag[ifile]    != -1) {
    status = SDendaccess(sds_id_nflag[ifile]);
    if (status != 0) {
      printf("Error ending access to n_flag sds for file: %ld\n", ifile);
      exit(-1);
    }
  }

  status = SDend(sd_id_file[ifile]);
  /*  printf("sd_id(C): %ld %ld\n", ifile, sd_id_file[ifile]);*/
  if (status != 0) {
    printf("Error ending access to file: %ld\n", ifile);
    exit(-1);
  }

  return 0;
}




int32 freeL2(l2_prod *l2_str)
{
  int32 i;
  int32 status;

  if (l2_str == NULL) {

    for (i=0; i<MAXNFILES; i++) {
      if (databuf[i] != NULL) free(databuf[i]);
      if (prodlist[i] != NULL) free(prodlist[i]);
    }

  } else {

    if (l2_str->geointerp == 1) {
      free(l2_str->lon_cntl);
      free(l2_str->lat_cntl);
      free(l2_str->spline_arr);
      free(l2_str->cntl_pnts);
    }

    free(l2_str->l2_data);
    free(l2_str->geoloc);

    if (l2_str->l2_flags != NULL) free(l2_str->l2_flags);
    if (l2_str->flagnames != NULL) free(l2_str->flagnames);
  }



  return 0;
}





int32 findprod(l2_prod *l2_str, char* prodname)
{
  int32 i;

  for (i=0; i<l2_str->nprod; i++) {
    if (strcmp(l2_str->prodname[i], prodname) == 0) return i;
  }

  return -1;
}


int32 readL2meta(meta_l2Type *meta_l2, int32 ifile)
{
  int32 attr_index;
  int32 dtype;
  int32 count;
  int32 sd_id;
  char buf[256];

  sd_id = sd_id_file[ifile];

  if ((attr_index = SDfindattr(sd_id, TITLE)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->title = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->title);
  }

  if ((attr_index = SDfindattr(sd_id, STATION)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->station = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->station);
  }

  if ((attr_index = SDfindattr(sd_id, INFILES)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->infiles = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->infiles);
  }

  if ((attr_index = SDfindattr(sd_id, SENNME)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->sensor_name = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->sensor_name); 
  }

  if ((attr_index = SDfindattr(sd_id, DCENTER)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->data_center = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->data_center); 
  }


  SDreadattr(sd_id, SDfindattr(sd_id, STLAT), (VOIDP)&meta_l2->station_lat);
  SDreadattr(sd_id, SDfindattr(sd_id, STLON), (VOIDP)&meta_l2->station_lon);
  SDreadattr(sd_id, SDfindattr(sd_id, SNCNTR), (VOIDP)&meta_l2->ncrec); 
  SDreadattr(sd_id, SDfindattr(sd_id, NFREC), (VOIDP)&meta_l2->nfrec);
  SDreadattr(sd_id, SDfindattr(sd_id, PCTFLAG), (VOIDP)meta_l2->flags_pc);

  if ((attr_index = SDfindattr(sd_id, CTIME)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->ctime = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->ctime);
  }

  if ((attr_index = SDfindattr(sd_id, NTIME)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->ntime = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->ntime);
  }

  if ((attr_index = SDfindattr(sd_id, SNODE)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->snode = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->snode);
  }

  if ((attr_index = SDfindattr(sd_id, ENODE)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->enode = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->enode); 
  }

  if ((attr_index = SDfindattr(sd_id, MISSION)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->mission = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->mission); 
  }

  if ((attr_index = SDfindattr(sd_id, MSNCHAR)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->mission_char = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->mission_char); 
  }

  if ((attr_index = SDfindattr(sd_id, SENSOR)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->sensor = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->sensor);
  }

  if ((attr_index = SDfindattr(sd_id, SNSCHAR)) != -1) {
    SDattrinfo(sd_id, attr_index, buf, &dtype, &count);
    meta_l2->sensor_char = (char *) calloc(count, sizeof(char));
    SDreadattr(sd_id, attr_index, (VOIDP)meta_l2->sensor_char);
  }

  SDreadattr(sd_id, SDfindattr(sd_id, ORBNUM), (VOIDP)&meta_l2->orbnum);
  SDreadattr(sd_id, SDfindattr(sd_id, CLAT),  (VOIDP)&meta_l2->sclat);
  SDreadattr(sd_id, SDfindattr(sd_id, CLON),  (VOIDP)&meta_l2->sclon);
  SDreadattr(sd_id, SDfindattr(sd_id, SCSOL_Z), (VOIDP)&meta_l2->scsol_z);
  SDreadattr(sd_id, SDfindattr(sd_id, ULLAT), (VOIDP)&meta_l2->ullat);
  SDreadattr(sd_id, SDfindattr(sd_id, ULLON), (VOIDP)&meta_l2->ullon);
  SDreadattr(sd_id, SDfindattr(sd_id, URLAT), (VOIDP)&meta_l2->urlat);
  SDreadattr(sd_id, SDfindattr(sd_id, URLON), (VOIDP)&meta_l2->urlon);
  SDreadattr(sd_id, SDfindattr(sd_id, LLLAT), (VOIDP)&meta_l2->lllat);
  SDreadattr(sd_id, SDfindattr(sd_id, LLLON), (VOIDP)&meta_l2->lllon);
  SDreadattr(sd_id, SDfindattr(sd_id, LRLAT), (VOIDP)&meta_l2->lrlat);
  SDreadattr(sd_id, SDfindattr(sd_id, LRLON), (VOIDP)&meta_l2->lrlon);
  SDreadattr(sd_id, SDfindattr(sd_id, NLAT),  (VOIDP)&meta_l2->northlat);
  SDreadattr(sd_id, SDfindattr(sd_id, SLAT),  (VOIDP)&meta_l2->southlat);
  SDreadattr(sd_id, SDfindattr(sd_id, WLON),  (VOIDP)&meta_l2->westlon);
  SDreadattr(sd_id, SDfindattr(sd_id, ELON),  (VOIDP)&meta_l2->eastlon);
  SDreadattr(sd_id, SDfindattr(sd_id, STCLAT), (VOIDP)&meta_l2->startclat);
  SDreadattr(sd_id, SDfindattr(sd_id, STCLON), (VOIDP)&meta_l2->startclon);
  SDreadattr(sd_id, SDfindattr(sd_id, ENDCLAT), (VOIDP)&meta_l2->endclat);
  SDreadattr(sd_id, SDfindattr(sd_id, ENDCLON), (VOIDP)&meta_l2->endclon);
  SDreadattr(sd_id, SDfindattr(sd_id, NODEL),  (VOIDP)&meta_l2->nodel);
  SDreadattr(sd_id, SDfindattr(sd_id, LAC_PX_ST),  (VOIDP)&meta_l2->pix_start);
  SDreadattr(sd_id, SDfindattr(sd_id, LAC_PX_SUBSAMP),  (VOIDP)&meta_l2->pix_sub);

  return 0;
}


int32 freeL2meta(meta_l2Type *meta_l2)
{
#define FREE(ptr) if((ptr) != 0x0) free(ptr);

  FREE(meta_l2->title);
  FREE(meta_l2->station);
  FREE(meta_l2->infiles);
  FREE(meta_l2->sensor_name); 
  FREE(meta_l2->data_center); 
  FREE(meta_l2->ctime);
  FREE(meta_l2->ntime);
  FREE(meta_l2->snode);
  FREE(meta_l2->enode);
  FREE(meta_l2->mission);
  FREE(meta_l2->mission_char); 
  FREE(meta_l2->sensor);
  FREE(meta_l2->sensor_char); 

  return 0;
}


int32 getL3units(l2_prod *l2_str, int32 ifile, char *l3b_prodname, char *units)
{
  intn i;
  int32 sds_id;

  char bufnum[128];
  char bufden[128];
  char* char_ptr;

  char_ptr = strchr(l3b_prodname, '/');
  if (char_ptr != NULL) *char_ptr = 0;

  for (i=0; i<l2_str[ifile].nprod; i++) {
    if (strcmp(l3b_prodname, l2_str[ifile].prodname[i]) == 0) {
      sds_id = sds_id_prod[ifile][i];
      SDreadattr(sds_id,SDfindattr(sds_id,"units"), (VOIDP) bufnum);
      break;
    }
  }

  if (char_ptr != NULL) {
    for (i=0; i<l2_str[ifile].nprod; i++) {
      if (strcmp(char_ptr+1, l2_str[ifile].prodname[i]) == 0) {
	sds_id = sds_id_prod[ifile][i];
	SDreadattr(sds_id,SDfindattr(sds_id,"units"), (VOIDP) bufden);
	break;
      }
    }

    if (strcmp(bufnum, bufden) == 0) {
      strcpy(units, "dimensionless");
    }
    else if (strcmp(bufnum, "dimensionless") == 0) {
      strcpy(units, "1 / ");
      strcat(units, bufden);
    }
    else if (strcmp(bufden, "dimensionless") == 0) {
      strcpy(units, bufnum);
    }
    else {
      strcpy(units, bufnum);
      strcat(units, " / ");
      strcat(units, bufden);
    }
  }
  else {
    strcpy(units, bufnum);
  }

  return 0;
}


/*

This code is based on the cubic spline interpolation code presented in:
Numerical Recipes in C: The Art of Scientific Computing
by
William H. Press,
Brian P. Flannery,
Saul A. Teukolsky, and
William T. Vetterling .
Copyright 1988 (and 1992 for the 2nd edition)

I am assuming zero-offset arrays instead of the unit-offset arrays
suggested by the authors.  You may style me rebel or conformist
depending on your point of view.

Norman Kuring	31-Mar-1999

*/

#include <stdio.h>
#include <stdlib.h>

#define MALLOC(ptr,typ,num) {                                           \
  (ptr) = (typ *)malloc((num) * sizeof(typ));                           \
  if((ptr) == NULL){                                                    \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n",      \
    __FILE__,__LINE__);                                                 \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}

void
spline(
float	x[],
float	y[],
int	n,
float	yp1,
float	ypn,
float	y2[]
){

  int	i,k;
  float	p,qn,sig,un,*u;

  MALLOC(u,float,n-1);

  if(yp1 > 0.99e30)
    y2[0] = u[0] = 0.0;
  else{
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for(i = 1; i < n-1; i++){
    sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig - 1.0)/p;
    u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
  }
  if(ypn > 0.99e30)
    qn = un = 0.0;
  else{
    qn = 0.5;
    un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
  }
  y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);
  for(k = n-2; k >= 0; k--){
    y2[k] = y2[k]*y2[k+1] + u[k];
  }

  free(u);
}

void
splint(
float	xa[],
float	ya[],
float	y2a[],
int	n,
float	x,
float	*y
){

  int		klo,khi,k;
  float		h,b,a;
  static int	pklo=0,pkhi=1;

  /*
  Based on the assumption that sequential calls to this function are made
  with closely-spaced, steadily-increasing values of x, I first try using
  the same values of klo and khi as were used in the previous invocation.
  If that interval is no longer correct, I do a binary search for the
  correct interval.
  */
  if(xa[pklo] <= x && xa[pkhi] > x){
    klo = pklo;
    khi = pkhi;
  }
  else{
    klo = 0;
    khi = n - 1;
    while(khi - klo > 1){
      k = (khi + klo) >> 1;
      if(xa[k] > x) khi = k;
      else          klo = k;
    }
  }

  h = xa[khi] - xa[klo];
  if(h == 0){
    fprintf(stderr,"-E- %s line %d: Bad xa input to function splint()\n",
            __FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;
  *y = a*ya[klo] + b*ya[khi] +
       ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
}


