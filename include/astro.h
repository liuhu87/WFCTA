/* $Source: /hires_soft/uvm2k/uti/astro.h,v $
 * $Log: astro.h,v $
 * Revision 1.1  2000/09/14 21:13:16  jeremy
 * Initial revision
 *
 *
 **************************************************************
 * These C routines are adapted from the FORTRAN routines     *
 * from the astro program, which were significantly adaped    *
 * from the book "astronomical formulae for calculators"      *
 * by jean meeus.                                             *
 **************************************************************
*/

#ifndef M_2PI
#define M_2PI (M_PI + M_PI)
#endif

#ifndef NULL
#define NULL ((void *)0)
#endif

/* Location of Big-H according to CT GPS clock */
#define DUGLON -1.7470048 /* 100.096d East */
#define DUGLAT 0.51281264 /* 29.382d North */
#define DUGALT 4410.0     /* meters */

/* Number of planets and stars we know about */
#define N_PLANETS 4
#define N_STARS 52

/* pi */
#define DPI 3.1415926535897932384626433832795028841971693993751
/* 2pi */
#define D2PI 6.2831853071795864769252867665590057683943387987502

/**************************************************************/
/* Convert date/time to/from Julian day                              */
/**************************************************************/
int astro_local(int local);

double astroTojday(int year, int month, int day, int hour, int min, int sec);

void astroFromjday(double jd,
  int *year, int *month, int *day, int *hour, int *min, int *sec, char *zone);

/**************************************************************/
/* Convert ecliptic coordinates to equatorial coordinates     */
/**************************************************************/
void astroEctoeq(float lam, float bet, float *ra, float *dec);

/**************************************************************/
/* Convert Julian day to local sidereal time (rads) at Dugway */
/**************************************************************/
double astroSidtim(double jd);
double astroLST(double jd);//Given by llma.

/**************************************************************/
/* Convert equatorial coordinates to horizon coordinates      */
/**************************************************************/
void slaDe2h (double ha, double dec, double phi, double *az, double *el);

/**************************************************************/
/* Computes the ecliptic longitude, latitude and time         */
/* derivatives (radians/day) of the sun (if moon = 0)         */
/*                             and moon (if moon = 1)         */
/**************************************************************/
void astroSunmon(double jd, int moon,
	    float *lam, float *bet, float *lamdot, float *betdot,
	    float *rsun);

/**************************************************************/
/* Return the time (in days) from jd until:                   */
/*   next astro-twilight ends   (stcode = 0)                  */
/*   next astro-twilight begins (stcode = 1)                  */
/*   next sunset                (stcode = 2)                  */
/*   next sunrise ???           (stcode = 3)                  */
/**************************************************************/
double astroSuntwi(double jd, int stcode);

/**************************************************************/
/* Return the time (in days) from jd until:                   */
/*   next moonset               (rscode = 0)                  */
/*   next moonrise              (rscode = 1)                  */
/**************************************************************/
double astroMoonrs(double jd, int rscode);

/**************************************************************/
/* Compute Run period, Start and Stop times for this/next run */
/*   twilight-ends to twilight-begins -> moon = 0             */
/*   twilight-ends to moon-rise       -> moon = 1             */
/*   moon-set to twilight-begins      -> moon = 2             */
/*   moon-set to moon-rise            -> moon = 3             */
/*   no dark period this night        -> moon = -1            */
/* Returns dark period (in days)                              */
/**************************************************************/
float astroRuntime(double jday, double *start, double *stop, int *moon);

/**************************************************************/
/* Return illuminated fraction of the moon's disk at time jd  */
/**************************************************************/
double astroIllfra(double jd);

/**************************************************************/
/* Compute direction cosines toward an object at equatorial   */
/* coordinates ra, dec (radians, equinox of date) at local    */
/* sidereal time st (radians).                                */
/*           x: east;   y: north;   z: upwards                */
/**************************************************************/
void objectCosines(float ra, float dec, float st,
		float *x, float *y, float *z);

/**************************************************************/
/* Adjust local coordinates (x, y, z) of object being viewed  */
/* to account approximatly for the effects of refraction in   */
/* Earths atmosphere. This is not particularly accurate       */
/**************************************************************/
void refractionAdjust(float *x, float *y, float *z);

/**************************************************************/
/* Account approximatly for the effects of parallax on the    */
/* position of the moon. It is assumed that the horizontal    */
/* parallax is 57.4 arcmin.                                   */
/**************************************************************/
void parallaxAdjust(float *x, float *y, float *z);

/**************************************************************/
/* Finds tube(s) that contain the coordinates (x, y, z) and   */
/* returns the number of tube found (<= maxt).                */
/**************************************************************/
int findTubes(float x, float y, float z, int maxt,
	       int mir[], int sclu[], int stub[]);

/**************************************************************/
/* Computes right ascension and declination of a planet at    */
/* time jd (julian day).                                      */
/**************************************************************/
void planetDirection(int p_code, double jd, float *ra, float *dec);

/**************************************************************/
/* Computes right ascension and declination of a fixed star   */
/* at time jd (julian day).                                   */
/**************************************************************/
void starDirection(int s_code, double jd, float *ra, float *dec);

/**************************************************************/
/* Planet and star names. starName also returns a flag        */
/* indicating whether the star is visable in the UV.          */
/**************************************************************/
char *planetName(int p_code);
char *starName(int s_code, int *uvob);
