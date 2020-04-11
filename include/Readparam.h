#include <iostream>
#include <fstream>
#include <string.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PM_ITEM_LIST   /* LIST OF ITEMS IN THE PARAMETERS FILE */  \
T(SIMMODE),           /*The simulation mode*/\
T(lhaaso_position),    /*lhaaso center coordinates*/\
T(ct_number),          /* number of CT, which is number of reflector files */ \
T(ct_position),          /* number of CT, which is number of reflector files */ \
T(mirror_size),    /*total size of the reflector mirror*/\
T(mirror_spot),   /*  imga of spot size caused by the roughness of the mirror (mm)*/ \
T(mirror_geomtry),  /* reflector geometry used */ \
T(mirror_pointerror), /* point error of each mirror */ \
T(wavelength), /*will the wavelength of photons be simulated */\
T(thin), /*will the thin option be used in the corsika*/\
T(filter), /*will the fiter be used*/\
T(fadc_bin),            /* total number of fadc bins */\
T(fadc_bin_length),     /* the bin length of each bin*/\
T(nsb),                 /* night sky background*/\
T(Aod_air),         /*extinction coefficient of air*/\
T(Aod_aerosol),         /*extinction coefficient of aerosol*/\
T(Scat_air),         /*scattering coefficient of air*/\
T(Scat_aerosol),         /*scattering coefficient of aerosol*/\
T(Laser_coo),            /*laser position*/\
T(Laser_dir),            /*laser generate direction*/\
T(Laser_intensity),      /*laser intensity*/\
T(Laser_wavelength),     /*laser wavelength*/\
T(Laser_frequency),      /*laser frequency (in Hz)*/\
T(Laser_pulsetime),      /*laser pulsetime*/\
T(Laser_spotrange),      /*laser spotrange (in mm)*/\
T(Laser_divergence),     /*laser divergence (in mrad)*/\
T(end_file)         /* end of the parameters file */
#define T(x)  x             // define T() as the name as it is

//#define T(x)  #x              // define T() as the string of x

enum ITEM_TYPE {
  PM_ITEM_LIST
};
#undef T

#define T(x)  #x              // define T() as the string of x
const char *const ITEM_NAMES[] = {
  PM_ITEM_LIST
};

#undef T

#define LINE_MAX_LENGTH  400
#define ITEM_MAX_LENGTH  40
#define PATH_MAX_LENGTH  256
#define MAX_NUMBER_OF_CTS 100
#define FALSE 0
#define TRUE 1

class WReadConfig
{
 
public:

    WReadConfig();
   ~WReadConfig();
   void readparam(char * filename);
   double GetLHAASOCoo(int i);
   int GetCTNumber();
   float GetCTPosition(int ict, int i);
   int GetMirrorGeometry();
   void GetMirrorSize(float *x , float *y);
   float GetMirrorSpot();
   int GetMirrorPointError();
   int GetMirrorPointErrorFlag();
   int GetFilterFlag();
   int GetThinFlag();
   int GetWLFlag();
   int GetFadcBins();
   int GetFadcLength();
   float GetNSB();

   float Getaod_air();
   float Getaod_aerosol();
   float Getscat_air();
   float Getscat_aerosol();
   float GetLaserCoo(int i);
   float GetLaserDir(int i);
   float GetLaserIntensity();
   float GetLaserIntensityErr();
   float GetLaserWavelength();
   float GetLaserWavelengthErr();
   float GetLaserFrequency();
   float GetLaserPulsetime();
   float GetLaserSpotrange(int i=0);
   float GetLaserDivergence(int i=0);
};
