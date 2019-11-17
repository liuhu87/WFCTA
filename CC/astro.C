/* $Source: /hires_soft/uvm2k/uti/astro.c,v $
 * $Log: astro.c,v $
 * Revision 1.1  2000/09/14 21:13:16  jeremy
 * Initial revision
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "astro.h"

/* All routines in this library use MODIFIED julian day */

/* modified julian day for 00:00:00 Jan 1, 1970 */
#define JD0 2440587.5
/* convert modified julian day to time_t */
#define timejd(jd) ((time_t)((jd - JD0) * 86.4e3))

static int local_time = 0;

/**************************************************************/
/* Set time conversion mode (0 = UTC, non-zero = LOCAL)       */
/**************************************************************/
int astroLocal(int local)
{
  int old = local_time;
  local_time = local;
  return old;
}

int astro_local_(int *local)
{
  return astroLocal(*local);
}

/**************************************************************/
/* Convert date/time to Julian day                              */
/**************************************************************/
double astroTojday(int year, int month, int day, int hour, int min, int sec)
{
  double jd; int a, b;
  if (month <= 2) year -= 1, month += 12;
  a = year / 100;
  b = 2 - a + (a / 4);
  jd = (double)((int)(365.25 * year) + (int)(30.6001 * (month + 1)) + day + b)
    + 1720994.5 + (((sec / 60.0) + min) / 60.0 + hour) / 24.0;
  if (local_time)
  { /* correct julian day if local time */
    time_t s = timejd(jd);
    struct tm *t = gmtime(&s);
    min = t->tm_hour * 60 + t->tm_min;
    t = localtime(&s);
    min -= t->tm_hour * 60 + t->tm_min;
    if (min > 720) min -= 720;
    else if (min <= -720) min += 720;
    jd += min / 1440.0;
  }
  return jd;
}

double astro_tojday_(int *year, int *month, int *day,
		     int *hour, int *min, int *sec)
{
  return astroTojday(*year, *month, *day, *hour, *min, *sec);
}

/**************************************************************/
/* Convert Julian day to date/time                            */
/**************************************************************/
void astroFromjday(double jd,
  int *year, int *month, int *day, int *hour, int *min, int *sec, char *zone)
{
  time_t s = timejd(jd);
  struct tm *t = (local_time ? localtime : gmtime)(&s);
  if (year != NULL) *year = t->tm_year + 1900;
  if (month != NULL) *month = t->tm_mon + 1;
  if (day != NULL) *day = t->tm_mday;
  if (hour != NULL) *hour = t->tm_hour;
  if (min != NULL) *min = t->tm_min;
  if (sec != NULL) *sec = t->tm_sec;
  if (zone != NULL) strncpy(zone, local_time ? tzname[(t->tm_isdst > 0)] : "UTC", 4);
}

void astro_fromjday_(double *jd, int *year, int *month, int *day,
		     int *hour, int *min, int *sec, char *zone)
{
  astroFromjday(*jd, year, month, day, hour, min, sec, zone);
}

/**************************************************************/
/* Convert ecliptic coordinates to equatorial coordinates     */
/**************************************************************/
void astroEctoeq(float lam, float bet, float *ra, float *dec)
{
  double local_ra, sin_dec;
  double sin_lam = sin(lam), cos_lam = cos(lam);
  double sin_bet = sin(bet), cos_bet = cos(bet);
  double sin_eps = 0.397784, cos_eps = 0.917479;
  local_ra = atan2((sin_lam * cos_eps - sin_bet * sin_eps / cos_bet), cos_lam);
  if (ra != NULL)
    *ra = (local_ra < 0.0) ? local_ra + M_2PI : local_ra;
  sin_dec = sin_bet * cos_eps + cos_bet * sin_eps * sin_lam;
  if (dec != NULL)
    *dec = atan(sin_dec / sqrt(1.0 - sin_dec * sin_dec));
}

void astro_ectoeq_(float *lam, float *bet, float *ra, float *dec)
{
  astroEctoeq(*lam, *bet, ra, dec);
}

/**************************************************************/
/* Convert Julian day to local sidereal time (rads) at Dugway */
/**************************************************************/
double astroSidtim(double jd)
{
  double jd0, t, theta;
  jd -= 2415020.0;
  jd0 = (int)(jd + 0.5) - 0.5;
  t = jd0 / 36525.0;
  theta = (jd - jd0) * 1.0027379 - DUGLON / M_2PI;
  theta += 0.276919398 + 100.0021359 * t + 1.075E-6 * t * t;
  theta -= (int)theta;
  return theta * M_2PI;
}

double astro_sidtim_(double *jd)
{
  return astroSidtim(*jd);
}

//////////////////////////Given by llma/////////////////////////
double astroLST(double jd)
{
  double lst;
  double mjd;
  mjd = jd - 2400000.5;
  lst = 0.671262 + 1.002737909*(mjd - 40000)+(-DUGLON)/D2PI;
  lst = lst-(int)(lst);
  lst = D2PI*lst;
  return lst;
}

double astro_lst_(double *jd)
{
  return astroLST(*jd);
}

/**************************************************************/
/* Computes the ecliptic longitude, latitude and time         */
/* derivatives (radians/day) of the sun (if moon = 0)         */
/*                             and moon (if moon = 1)         */
/**************************************************************/
void astroSunmon(double jd, int moon,
  float *lam, float *bet, float *lamdot, float *betdot, float *rsun)
{
  float local_lam, t, m, e;
  t = (jd - 2415020) / 36525.0;
  m = 0.99576619 + 99.997360 * t;
  m = M_2PI * (m - (int)m);
  if (moon == 0)
  { /* coordinates of the sun */
    float ecc, l, nu;
    int i;
    l = 0.77693522 + 100.00214 * t;
    l = M_2PI * (l - (int)l);
    ecc = 0.01675104 - 0.0000418 * t;
    e = m;
    for (i = 0; i < 5; ++i)
      e = m + ecc * sin(e);
    nu = 2.0 * atan((sin(0.5 * e) / cos(0.5 * e)) *
		    sqrt((1.0 + ecc) / (1.0 - ecc)));
    local_lam = l + nu - m;
    if (lamdot != NULL)
      *lamdot = 0.017202 * sqrt((1.0 + ecc) / (1.0 - ecc)) *
	pow(cos(0.5 * nu) / cos(0.5 * e), 2.0) / (1.0 - ecc * cos(e));
    if (bet != NULL) *bet = 0.0;
    if (betdot != NULL) *betdot = 0.0;
    /* rsun is the earth-sun distance */
    if (rsun != NULL) *rsun = 1.0000002 * (1.0 - ecc * cos(e));
  }
  else
  { /* coordinates of the moon */
    float ll, mm, d, f;
    e = 1.0 - 0.002495 * t;
    ll = 0.75120601 + 1336.855231 * t;
    ll = M_2PI * (ll - (int)ll);
    mm = 0.82251280 + 1325.552359 * t;
    mm = M_2PI * (mm - (int)mm);
    d = 0.97427079 + 1236.853095 * t;
    d = M_2PI * (d - (int)d);
    f = 0.03125247 + 1342.227848 * t;
    f = M_2PI * (f - (int)f);
    local_lam = ll + 0.109759 * sin(mm) + 0.022236 * sin(2.0 * d - mm)
      + 0.011490 * sin(2.0 * d) + 0.003728 * sin(2.0 * mm)
      - 0.003239 * e * sin(m) - 0.001996 * sin(2.0 * f)
      + 0.001026 * sin(2.0 * (d - mm)) + 0.000999 * e * sin(2.0 * d - m - mm)
      + 0.000931 * sin(2.0 * d + mm) + 0.000801 * e * sin(2.0 * d - m)
      + 0.000716 * e * sin(mm - m);
    if (bet != NULL)
      *bet = 0.089504 * sin(f) + 0.004897 * sin(mm - f)
	+ 0.004847 * sin(mm - f) + 0.003024 * sin(2.0 * d - f)
	+ 0.000967 * sin(2.0 * d + f - mm);
    if (lamdot != NULL)
      *lamdot = 0.229972 + 0.025028 * cos(mm) + 0.004392 * cos(2.0 * d - mm)
	+ 0.004889 * cos(2.0 * d) + 0.001700 * cos(2.0 * mm)
	- 0.000922 * cos(2.0 * f);
    if (betdot != NULL)
      *betdot = 0.020666 * cos(f) + 0.002248 * cos(mm + f)
	+ 0.002248 * cos(mm + f);
  }
  if (lam != NULL)
  {
    local_lam = local_lam / M_2PI;
    local_lam = local_lam - (int)local_lam + 2.0;
    *lam = M_2PI * (local_lam - (int)local_lam);
  }
}

void astro_sunmon_(double *jd, int *moon,
     float *lam, float *bet, float *lamdot, float *betdot, float *rsun)
{
  astroSunmon(*jd, *moon, lam, bet, lamdot, betdot, rsun);
}

/**************************************************************/
/* Return the time (in days) from jd until:                   */
/*   next astro-twilight ends   (stcode = 0)                  */
/*   next astro-twilight begins (stcode = 1)                  */
/*   next sunset                (stcode = 2)                  */
/*   next sunrise ???           (stcode = 3)                  */
/**************************************************************/
double astroSuntwi(double jd, int stcode)
{
  float lam, lamdot, bet, betdot;
  float ra, dec, sinh0, cosh0, h;
  astroSunmon(jd, 0, &lam, &bet, &lamdot, &betdot, NULL);
  astroEctoeq(lam, bet, &ra, &dec);
  cosh0 = -0.309 - sin(DUGLAT) * sin(dec);
  if (stcode == 2 || stcode == 3) cosh0 += 0.295;
  cosh0 /= cos(DUGLAT) * cos(dec);
  sinh0 = sqrt(1.0 - cosh0 * cosh0);
  if (stcode == 1 || stcode == 3) sinh0 = -sinh0;
  h = (atan2(sinh0, cosh0) - (astroSidtim(jd) - ra)) / M_2PI + 3.0;
  return h - (int)h;
}

double astro_suntwi_(double *jd, int *stcode)
{
  return astroSuntwi(*jd, *stcode);
}

/**************************************************************/
/* Return the time (in days) from jd until:                   */
/*   next moonset               (rscode = 0)                  */
/*   next moonrise              (rscode = 1)                  */
/**************************************************************/
double astroMoonrs(double jd, int rscode)
{
  float lam, lamdot, bet, betdot, st0;
  float moon, ra, dec, sinh0, cosh0;
  int i, flag;
  astroSunmon(jd, 1, &lam, &bet, &lamdot, &betdot, NULL);
  st0 = astroSidtim(jd);
  moon = 0.0;
  flag = 0;
  for (i = 0; i < 10; ++i)
  {
    astroEctoeq(lam + lamdot * moon, bet + betdot * moon, &ra, &dec);
    cosh0 = (0.0068 - sin(DUGLAT) * sin(dec)) / (cos(DUGLAT) * cos(dec));
    sinh0 = sqrt(1.0 - cosh0 * cosh0);
    if (rscode == 1) sinh0 *= -1.0;
    moon = (atan2(sinh0, cosh0) - (st0 - ra)) / M_2PI + 5.0;
    moon = (moon - (int)moon) / 1.0027379;
    if (i == 0 && moon > 0.8) flag = 1;
    else if (flag && moon < 0.2) moon += 0.99726958;
  }
  return moon;
}

double astro_moonrs_(double *jd, int *rscode)
{
  return astroMoonrs(*jd, *rscode);
}

/**************************************************************/
/* Compute Run period, Start and Stop times for this/next run */
/* if Illfra is more than 50%                                 */
/*   twilight-ends to twilight-begins -> moon = 0             */
/*   twilight-ends to moon-rise       -> moon = 1             */
/*   moon-set to twilight-begins      -> moon = 2             */
/*   moon-set to moon-rise            -> moon = 3             */
/*   no dark period this night        -> moon = -1            */
/* if Illfra is less than 50%                                 */
/*   twilight-ends to twilight-begins -> moon = 4             */
/* Returns dark period (in days)                              */
/**************************************************************/
float astroRuntime(double jday, double *start, double *stop, int *moon)
{
  double twil_end, twil_beg, moon_set, moon_ris, period;
  double lstart, lstop;
  int Illfra;
  int lmoon;
  twil_end = jday + astroSuntwi(jday, 0);
  twil_beg = jday + astroSuntwi(jday, 1);
  /* make sure twilight-ends is before next twilight-begins */
  if (twil_beg < twil_end)
  {
    jday = twil_end - 1.25; /* 30 hours before next twilight-ends */
    twil_end = jday + astroSuntwi(jday, 0);
  }

  Illfra = astroIllfra(twil_end);
  if(Illfra>-1)
  {
    moon_set = twil_end + astroMoonrs(twil_end, 0);
    moon_ris = twil_end + astroMoonrs(twil_end, 1);

    /* sequence of events: twilight-ends is always first */
    if (twil_beg < moon_set && moon_set < moon_ris)
    { /* no dark period */
      lmoon = -1;
      lstart = lstop = 0.0;
      period = 0.0;
    }
    else
    {
      if (twil_beg < moon_ris && moon_ris < moon_set)
        lmoon = 0;     /* twilight-ends -> twilight-begins */
      else if (moon_ris < twil_beg && twil_beg < moon_set)
        lmoon = 1;     /* twilight-ends -> moon-rise */
      else if (moon_ris < moon_set && moon_set < twil_beg)
        /* Two dark periods this night; choose the longest */
        lmoon = ((moon_ris - twil_end) > (twil_beg - moon_set)) ? 1 : 2;
      else if (moon_set < twil_beg && twil_beg < moon_ris)
        lmoon = 2;     /* moon-set -> twilight-begins */
      else
        lmoon = 3;     /* moon-set -> moon-rise */
      lstart = (lmoon & 2) ? moon_set : twil_end;
      lstop = (lmoon & 1) ? moon_ris : twil_beg;
    }
  }
  else
  {/* if Illfra is less than 50%, SiPM can work at moon night*/
   lstart = twil_end;
   lstop = twil_beg;
   lmoon = 4;
  }
  period = lstop - lstart;

  if (start != NULL) *start = lstart;
  if (stop != NULL) *stop = lstop;
  if (moon != NULL) *moon = lmoon;
  return period;
}

float astro_runtime_(double *jday, double *start, double *stop, int *moon)
{
  return astroRuntime(*jday, start, stop, moon);
}

/**************************************************************/
/* Return illuminated fraction of the moon's disk at jd       */
/**************************************************************/
double astroIllfra(double jd)
{
  double t;
  float m, m1, d, i;
  t = (jd - 2415020.0) / 36525.0;
  m = 6.25658358 + 628.301947 * t;
  m1 = 5.16800034 + 8328.69110 * t;
  d = 6.12152394 + 7771.37772 * t;
  i = M_PI - d - 0.1098 * sin(m1) + 0.03665 * sin(m) -
    0.0222 * sin(2.0 * d - m1) - 0.0115 * sin(2.0 * d) -
    0.00374 * sin(2.0 * m1) - 0.00195 * sin(d);
  return (50.0 * (1.0 + cos(i)) + 0.5);
}

double astro_illfra_(double *jd)
{
  return astroIllfra(*jd);
}

/**************************************************************/
/* Convert equatorial coordinates to horizon coordinates      */
/* Copyright P.T.Wallace.  All rights reserved.               */
/* Copied from slalib.                                        */
/**************************************************************/
void slaDe2h ( double ha, double dec, double phi, double *az, double *el )
{
   double sh, ch, sd, cd, sp, cp, x, y, z, r, a;
   
/* Useful trig functions */
   sh = sin ( ha );
   ch = cos ( ha );
   sd = sin ( dec );
   cd = cos ( dec );
   sp = sin ( phi );
   cp = cos ( phi );

/* Az,El as x,y,z */
   x = - ch * cd * sp + sd * cp;
   y = - sh * cd;
   z = ch * cd * cp + sd * sp;

/* To spherical */
   r = sqrt ( x * x + y * y );
   a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
   *az = ( a < 0.0 ) ? a + D2PI : a;
   *el = atan2 ( z, r );
}
