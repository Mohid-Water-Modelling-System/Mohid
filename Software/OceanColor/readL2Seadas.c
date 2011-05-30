
#include "readL2scan.h"
#include <ctype.h>
#include <string.h>



meta_l2Type *meta_l2;      
l2_prod     *l2_str;

float       *Data2D;
float       *Lat2D;
float       *Lon2D;
long        *inx;
int32       **l2_flags;
float       ***l2vec;

//void __stdcall Ucase2( char *string, int length )
//{
//    char *ptr;

//    for (ptr = string; *ptr; ptr++)
//        *ptr = toupper( *ptr );
//}


void __stdcall OpenL2Seadas(char *filename, int length) {

 
   long        nl2files,i,iscan,n, ipix;
   float       lat,lon;
   char ptr;

   void lonlat2geovec( float lon, float lat, float *v );
   
   int32       nprod,flags,nrec,nsamp;
   float32     longitude,latitude,data;
   float  auxz[10];
     
 


        nl2files = 1;

		if ( (inx = (long *)malloc(nl2files*sizeof(long))) == NULL ) {
           	printf("-E- %s: Error allocating memory for L2 file index.\n");
           	exit(1);
    	    }
	    
		
		if ( (l2vec = (float ***)malloc(nl2files*sizeof(float **))) == NULL ) {
           printf("-E- %s: Error allocating memory for L2 geocentric pixel location vector.\n");
           exit(1);
    	} 
    	if ( (l2_flags = (int32 **)malloc(nl2files*sizeof(int32 *))) == NULL ) {
           printf("-E- %s: Error allocating memory for L2 flags.\n");
           exit(1);
    	}
        
		
		
		
		for (i=0; i<nl2files; i++) inx[i] = i;

   
      	if ( (l2_str = (l2_prod *)malloc(nl2files*sizeof(l2_prod))) == NULL ) {
            printf("-E- %s: Error allocating memory for reading L2 data.\n");
            exit(1);
      	}
      	
		

		for (i=0; i<nl2files; i++) {
    	    l2_str[i].nrec = 0;
    	    l2_str[i].nsamp = 0;
	    if ( openL2(filename, 0x0, &l2_str[i]) ) {
        	printf("-E- %s: readL2Seadas.c err01 - reading L2 data %s.\n");
        	exit(1);
            }
		
		    n = l2_str[inx[i]].nsamp*l2_str[inx[i]].nrec;
           if ( (l2vec[i] = (float **)malloc(n*sizeof(float *))) == NULL ) {
         	printf("-E- %s: Error allocating memory for L2 geocentric pixel location.\n");
         	exit(1);
           }
           if ( (l2_flags[i] = (int32 *)malloc(n*sizeof(int32))) == NULL ) {
         	printf("-E- %s: Error allocating memory for L2 flags.\n");
         	exit(1);
           }

		   if ( (Data2D = (float *)malloc(n*sizeof(float))) == NULL ) {
         	printf("-E- %s: Error allocating memory for Data2D.\n");
         	exit(1);
           }

		   if ( (Lon2D = (float *)malloc(n*sizeof(float))) == NULL ) {
         	printf("-E- %s: Error allocating memory for Lon2D.\n");
         	exit(1);
           }

		   if ( (Lat2D = (float *)malloc(n*sizeof(float))) == NULL ) {
         	printf("-E- %s: Error allocating memory for Lat2D.\n");
         	exit(1);
           }

		   
           if ( (meta_l2 = (meta_l2Type *)malloc(1*sizeof(meta_l2Type))) == NULL ) {
              printf("-E- %s: Error allocating memory for L2 metadata.\n");
              exit(0);
		   }  
 
	       if ( readL2meta(&(meta_l2[i]), i) ) {
                printf("-E- %s: Error reading L2 metadata %s.\n");
                exit(0);
		   }

 }
}


void ReadDataL2(int Mini, int Maxi, int ProdId) 
{

int iscan,i,ipix;
float lat,lon;		
	

		for (iscan = Mini-1; iscan <= Maxi-1; iscan++) {
      
		  if (readL2(&(l2_str[0]), 0, iscan, ProdId)) {
                   printf("%s: Error: Cannot read L2 data file %s.\n", 
                        l2_str[inx[i]].filename);
                   exit(1);
      	  	
		  }
		
      		for (ipix = 0; ipix < l2_str[0].nsamp; ipix++) {
		
			   Data2D[iscan*l2_str[0].nsamp+ipix] = l2_str[0].l2_data[ProdId*l2_str->nsamp + ipix];

      	}
    }

    
	 
		  
}

void ReadFirstL2() 
{

int iscan,i,ipix;
float lat,lon;		
	
	i=0;

		for (iscan = 0; iscan < l2_str[inx[i]].nrec; iscan++) {
      
		  if (readInicL2(&(l2_str[inx[i]]), inx[i], iscan)) {
                   printf("%s: Error: Cannot read L2 data file %s.\n", 
                        l2_str[inx[i]].filename);
                   
				   exit(1);
				   
      	  	
		  }
	
		
      		for (ipix = 0; ipix < l2_str[inx[i]].nsamp; ipix++) {
		
	           lat = l2_str[inx[i]].latitude[ipix];
	           lon = l2_str[inx[i]].longitude[ipix];
			  
		   	         
		   /*  get the L2 flags                   */
	           l2_flags[i][iscan*l2_str[inx[i]].nsamp+ipix] = l2_str[inx[i]].l2_flags[ipix];
		      
		       Lon2D[iscan*l2_str[inx[i]].nsamp+ipix]  = l2_str[inx[i]].longitude[ipix];
               Lat2D[iscan*l2_str[inx[i]].nsamp+ipix]  = l2_str[inx[i]].latitude[ipix];
			   
			   /*  get the L2 geocentric coordinates  */
	           //if ( (l2vec[i][iscan*l2_str[inx[i]].nsamp+ipix] = (float *)malloc(3*sizeof(float))) == NULL ) {
         	//	printf("-E- %s: Error allocating memory for L2 geocentric pixel location vector.\n" );
         	//	exit(1);
    	      //     }
	           
			    //lonlat2geovec( lon, lat, l2vec[i][iscan*l2_str[inx[i]].nsamp+ipix] );

      	}
    }

    
	 
		  
}

void __stdcall GetProdName(char *ProductName, int length) {

	int i;
	i=0;
    ProductName=l2_str->prodname[i];

}

void GetLatLon (int j, int i, float *lat, float *lon) {
int geointerp;

   	geointerp = l2_str->geointerp;
	
	if (geointerp==1 || geointerp==0) {   // modis
		*lon    = Lon2D[(i-1)*l2_str[0].nsamp+j-1];
        *lat    = Lat2D[(i-1)*l2_str[0].nsamp+j-1];
    } else if (geointerp==2)  { // seawifs
		*lon    = Lon2D[(i-1)*l2_str[0].nsamp+j-1+3];
        *lat    = Lat2D[(i-1)*l2_str[0].nsamp+j-1];
    } else {

		printf("Error! GeoInterp Not Defined");
                   exit(1);
	}
	
}


void GetData (int i, int j, int idprod, float *data, int32 *flag) {
    *data   = Data2D[(j-1)*l2_str[0].nsamp+i-1];
//    *flag   = l2_str[0].l2_flags[(j-1)*l2_str[0].nsamp+i-1];
}

void GetInfo (int *i, int *j) {
  *i = l2_str[0].nrec;
  *j = l2_str[0].nsamp;
}


void GetMetaData (float *ullat,   
float *ullon,   	
float *urlat,   	
float *urlon,   	
float *lllat,   	
float *lllon,   	
float *lrlat,   	
float *lrlon,	
float *northlat,	
float *southlat,	
float *westlon, 	
float *eastlon,
int   *columns,
int   *lines,
int   *nprod,
int   *geointerp){ 	
	  


	   *ullat = meta_l2->ullat;
       *ullon    = meta_l2->ullon;
       *urlat    = meta_l2->urlat;
       *urlon    = meta_l2->urlon;
       *lllat    = meta_l2->lllat;
       *lllon    = meta_l2->lllon;
       *lrlat    = meta_l2->lrlat;
       *lrlon    = meta_l2->lrlon;
       *northlat = meta_l2->northlat;
       *southlat = meta_l2->southlat;
       *westlon  = meta_l2->westlon;
       *eastlon  = meta_l2->eastlon;
	   *lines    = l2_str[0].nrec;
	   *columns  = l2_str[0].nsamp;
	   *nprod    = l2_str->nprod;
	   *geointerp= l2_str->geointerp;
}
void closeFile()
{

		closeL2(&l2_str[0], 0);
	   	freeL2(&l2_str[0]);
		

}

