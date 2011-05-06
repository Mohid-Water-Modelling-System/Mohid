#include <math.h>
#include <malloc.h>
#include <float.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>

void lonlat2geovec( float lon, float lat, float *v )
{
#define  PI      3.141592653589793
    double Re = 6378.137;                          /* Earth radius in km */
    double f  = 1/298.257223563;                   /* Earth flattening factor */
    double rlon, rlat, phi, R;

    rlon = (double)(lon*PI/180.0);                 /* geodetic lon in radians */
    rlat = (double)(lat*PI/180.0);                 /* geodetic lat in radians */
    phi  = atan(tan(rlat)*(1-f)*(1-f));            /* geocentric latitude     */
    R    = Re*(1.0-f)/sqrt(1.0-(2.0-f)*f*(cos(phi)*cos(phi)));  /* dist to Earth surface  */

    *v = (float)(R*cos(phi)*cos(rlon));
    *(v+1) = (float)(R*cos(phi)*sin(rlon));
    *(v+2) = (float)(R*sin(phi));

}