#include "hdf.h"
#include "mfhdf.h"

#define byte unsigned char

#define MAXNFILES 525

typedef struct l2prod_struct {

  int32    fileindex;

  char     filename[128];

  int32    nrec;
  int32    nsamp;

  int16    syear;
  int16    sday;
  int32    smsec;
  int16    eyear;
  int16    eday;
  int32    emsec;
  int32    orbit;
  char     dtype[8];

  int32    ntilts;
  int16    tilt_flags[20];
  int16    tilt_ranges[2][20];

  char     *flagnames;
  int32    flagmask;

  int32    year;
  int32    day;
  int32    msec;

  float32  *geoloc;
  float32  *longitude;
  float32  *latitude;
  float32  *geonav[6];

  float32  *lon_cntl;
  float32  *lat_cntl;
  float32  *cntl_pnts;
  float32  *spline_arr;

  int32    nprod;
  char     *prodname[100];

  float32  *l2_data;
  int32    *l2_flags;

  byte     eng_qual[4];
  byte     s_flags[4];
  int32    nflag[8];

  int32    geointerp;

} l2_prod;


typedef struct meta_l2Struct {
        char    *product_name;  /* ATTR Product name(file name)         */
        char    *title;
        char    *data_center;   /* ATTR data_center, processing center  */
        char    *station;       /* ATTR station                         */
        float   station_lat;    /* ATTR station latitude                */
        float   station_lon;    /* ATTR station longitude               */
        char    *mission;       /* ATTR mission                         */
        char    *mission_char;  /* ATTR Mission Characteristics         */
        char    *sensor_name;   /* ATTR sensor name                     */
        char    *sensor;        /* ATTR sensor                          */
        char    *sensor_char;   /* ATTR Sensor Characteristics          */
        char    *sw_id;         /* ATTR Software ID                     */
        char    *infiles;       /* ATTR Input files                     */
        char    *stime;         /* ATTR Start time                      */
        char    *etime;         /* ATTR End time                        */
        char    *ctime;         /* ATTR scene center time               */
        char    *ntime;         /* ATTR Node crossing time              */
        char    *snode;         /* ATTR Start Node                      */
        char    *enode;         /* ATTR End Node                        */
        int     orbnum;         /* ATTR orbit number                    */
        char    *norad1;        /* ATTR NORAD elements, first line      */
        char    *norad2;        /* ATTR NORAD elements, second line     */
        int     pix_start;      /* ATTR LAC Pixel Start Number          */
        int     pix_sub;        /* ATTR LAC Pixel Subsampling           */
        int     ncrec;          /* ATTR scene center scan line          */
        int     nfrec;          /* ATTR number of filled scan line      */
        byte    ff_mis;         /* ATTR FF missing frames               */
        byte    sd_mis;         /* ATTR SDPS missing frames             */
        float   flags_pc[32];   /* MFSD % data for each quality flag    */
        char    *lat_units;     /* ATTR Latitude units                  */
        char    *lon_units;     /* ATTR Longitude units                 */
        float   sclat;          /* ATTR scene center latitude           */
        float   sclon;          /* ATTR scene center longitude          */
        float   scsol_z;        /* ATTR scene center solar zenith angle */
        float   ullat;          /* ATTR scene upper left latitude       */
        float   ullon;          /* ATTR scene upper left longitude      */
        float   urlat;          /* ATTR scene upper right latitude      */
        float   urlon;          /* ATTR scene upper right longitude     */
        float   lllat;          /* ATTR scene lower left latitude       */
        float   lllon;          /* ATTR scene lower left longitude      */
        float   lrlat;          /* ATTR scene lower right latitude      */
        float   lrlon;          /* ATTR scene lower right longitude     */
        float   northlat;       /* ATTR Northernmost latitude           */
        float   southlat;       /* ATTR Southernmost latitude           */
        float   westlon;        /* ATTR Westernmost longitude           */
        float   eastlon;        /* ATTR Easternmost longitude           */
        float   startclat;      /* ATTR Start Center Latitude           */
        float   startclon;      /* ATTR Start Center Longitude          */
        float   endclat;        /* ATTR End Center Latitude             */
        float   endclon;        /* ATTR End Center Longitude            */
        float   nodel;          /* ATTR Orbit node longitude            */
        int     ntilts;         /* MFSD Sensor Tilt                     */
        /* Calibration Vgroup */
        short   entry_year;
        short   entry_day;
        short   ref_year;
        short   ref_day;
        short   ref_minute;
} meta_l2Type;


#define PNAME           "Product Name"
#define TITLE           "Title"
#define LEGEND          "Legend"
#define DCENTER         "Data Center"
#define MISSION         "Mission"
#define MSNCHAR         "Mission Characteristics"
#define SENNME          "Sensor Name"
#define SENSOR          "Sensor"
#define SNSCHAR         "Sensor Characteristics"
#define REPLACES        "Replacement Flag"
#define SOFTID          "Software ID"
#define PTIME           "Processing Time" 
#define INFILES         "Input Files"
#define PROC_CON        "Processing Control"
#define PROC_LOG        "Processing Log"

/*  Constants defined for attributes defined in sec.5.3.2 of product specs */
  
#define STATION    "Station Name"
#define STLAT      "Station Latitude"
#define STLON      "Station Longitude"
#define DTYPE      "Data Type"
#define NPIX       "Pixels per Scan Line"
#define NSAMP      "Parent Pixels per Scan Line"
#define NSCAN      "Number of Scan Lines"
#define NREC       "Parent Number of Scan Lines"
#define SNCNTR     "Scene Center Scan Line"
#define NFREC      "Filled Scan Lines"
#define FFLAG1     "FF Missing Frames"
#define FFLAG2     "SDPS Missing Frames"
#define MASKNAMES  "Mask Names"
#define PCTFLAG    "Flag Percentages"
#define FLAGNAMES  "EPSILON1,LAND1,ANCIL1,SUNGLINT1,HIGHLT1,SATZEN1,COASTZ1,NEGLW1,STRAYLIGHT1,CLDICE1,COCCOLITH1,TRUBIDW1,SOLZEN1,HIGHTAU1,LOWLW1,CHLOR1"

/*  Constants defined for attributes defined in sec.5.3.3 of product specs */
#define STIME       "Start Time"
#define ETIME       "End Time"
#define CTIME       "Scene Center Time"
#define NTIME       "Node Crossing Time"
#define SYEAR       "Start Year"
#define SDAY        "Start Day"
#define SMSEC       "Start Millisec"
#define EYEAR       "End Year"
#define EDAY        "End Day"
#define EMSEC       "End Millisec"
#define SNODE       "Start Node"
#define ENODE       "End Node"
#define ORBNUM      "Orbit Number"
#define NORAD1      "NORAD Line 1"
#define NORAD2      "NORAD Line 2"

/*  Constants defined for attributes defined in sec.5.3.4 of product specs */
#define LATUNITS                "Latitude Units"
#define LONUNITS                "Longitude Units"
#define ULLAT                   "Upper Left Latitude"
#define ULLON                   "Upper Left Longitude"
#define URLAT                   "Upper Right Latitude"
#define URLON                   "Upper Right Longitude"
#define LLLAT                   "Lower Left Latitude"
#define LLLON                   "Lower Left Longitude"
#define LRLAT                   "Lower Right Latitude"
#define LRLON                   "Lower Right Longitude"
#define CLAT                    "Scene Center Latitude"
#define CLON                    "Scene Center Longitude"
#define SCSOL_Z                 "Scene Center Solar Zenith"
#define NLAT                    "Northernmost Latitude"
#define SLAT                    "Southernmost Latitude"
#define WLON                    "Westernmost Longitude"
#define ELON                    "Easternmost Longitude"
#define STCLAT                  "Start Center Latitude"
#define STCLON                  "Start Center Longitude"
#define ENDCLAT                 "End Center Latitude"
#define ENDCLON                 "End Center Longitude"
#define NODEL                   "Orbit Node Longitude"

/*  Constants defined for attributes defined in sec.5.3.4 of product specs */
#define PARAM                   "Parameter"
#define UNITS                   "Units"
#define PX_START                "Start Pixel"
#define LAC_PX_ST               "LAC Pixel Start Number"
#define PX_SUBSAMP              "Pixel Subsampling Rate"
#define LAC_PX_SUBSAMP          "LAC Pixel Subsampling" 
#define PX_NUM                  "Pixels per Scan Line"
#define SC_START                "Start Scan" 
#define SC_SUBSAMP              "Scan Subsampling Rate"
#define SC_NUM                  "Number of Scan Lines"
#define PX_LL_NUM               "Pixel Coordinates"
#define SC_LL_NUM               "Scan Coordinates"
#define SC_TYPE                 "Scaling"
#define SC_EQN                  "Scaling Equation"
#define BASE                    "Base"
#define SLOPE                   "Slope"
#define INTERCEPT               "Intercept"
#define SCALE_OFF               "Scale Offset"
  
/*  Constants defined for attributes defined in sec.5.4.1 of product specs */
#define PX_LL_FIRST             "px_ll_first"
#define PX_LL_FST_ATTR          "Lat/lon of pixels along first scan line"
#define PX_LL_LAST              "px_ll_last"
#define PX_LL_LST_ATTR          "Lat/lon of pixels along last scan line"
#define SC_LL_FIRST             "sc_ll_first"
#define SC_LL_FST_ATTR          "Lat/lon of starts of scan lines"
#define SC_LL_LAST              "sc_ll_last"
#define SC_LL_LST_ATTR          "Lat/lon of ends of scan lines"

/*  Constants defined for attributes defined in sec.5.4.2 of product specs */
#define NTILTS                  "ntilts"
#define NTILTS_NAME             "Number of scene tilt states"
#define T_FLAGS                 "tilt_flags"
#define T_FLAGS_NAME            "Tilt indicators"
#define T_RANGES                "tilt_ranges"
#define T_RANGES_NAME           "Scan-line number ranges of scene tilt states"
#define T_LATS                  "tilt_lats"
#define T_LATS_NAME             "Latitudes of tilt-range scan line end points"
#define T_LONS                  "tilt_lons"
#define T_LONS_NAME             "Longitudes of tilt-range scan line end points"

/* Constants defined for navigation data */
#define ORBVEC           "orb_vec"
#define ORBVEC_NAME      "Orbit position vector at scan line time"
#define ORBVEC_UNITS     "kilometers"
#define LVERT            "l_vert"
#define LVERT_NAME       "Local vertical vector in ECEF frame"
#define SUNREF           "sun_ref"
#define SUNREF_NAME      "Reference Sun vector in ECEF frame"
#define ATTANG           "att_ang"
#define ATTANG_NAME      "Computed yaw, roll, pitch"
#define SENMAT           "sen_mat"
#define SENMAT_NAME      "ECEF-to-sensor-frame matrix"
#define SCANELL          "scan_ell"
#define SCANELL_NAME     "Scan-track ellipse coefficients"
#define NFLAG            "nflag"
#define NFLAG_NAME       "Navigation flags"


#define NTILTS          "ntilts"
#define TILT_FLAGS      "tilt_flags"
#define TILT_RANGES     "tilt_ranges"
#define TILT_LATS       "tilt_lats"
#define TILT_LONS       "tilt_lons"
#define MSEC            "msec"
#define ENG_QUAL        "eng_qual"
#define S_FLAGS         "s_flags"
#define NDVI_DATA       "NDVI"
#define L2_FLAGS        "l2_flags"
#define ORBVEC          "orb_vec"
#define LVERT           "l_vert"
#define SUNREF          "sun_ref"
#define ATTANG          "att_ang"
#define SENMAT          "sen_mat"
#define SCANELL         "scan_ell"
#define NFLAG           "nflag"


/* Prototypes */
int32 openL2(char *, char *, l2_prod *);
int32 reopenL2(int32, l2_prod *);
int32 readL2(l2_prod *, int32, int32, int32);
int32 readlonlat(l2_prod *, int32, int32 *, int32 *);
int32 closeL2(l2_prod *, int32);
int32 freeL2(l2_prod *);
int32 findprod(l2_prod *, char* );
int32 readL2meta(meta_l2Type *, int32);
int32 freeL2meta(meta_l2Type *);
int32 getL3units(l2_prod *, int32, char *, char *);

void spline(float [], float [], int, float, float, float []);
void splint(float [], float [], float [], int, float, float*);

int32 FGEONAV (float pos[3], float rm[8], float coef[6], float sun[3],
		   int *nsta, int *ninc, int *npix, float latitude[1284], float longitude[1284]);


