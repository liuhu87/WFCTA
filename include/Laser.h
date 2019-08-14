#ifndef __Laser_H__
#define __Laser_H__
#include "common.h"
#include "WFTelescope.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TFile.h"

class WFCTAEvent;

/*
The class for absorbtion and scatter process calculation
*/
class Atmosphere {
   public:
      ///extinction coefficient
      static double aod_air;
      static double aod_aerosol;
      static double scat_air;
      static double scat_aerosol;
      static TGraph* gRayScatAngle;
      static TGraph* gMieScatAngle;
      static double scale;
   public:
      void Init(int seed=0);
      void Release();
      Atmosphere(int seed=0) {Init(seed);}
      ~Atmosphere() {Release();}
      static void SetParameters(char* filename=0);
      static double ZDependence(double z,int type=0);
      static double DeltaZ(double z);
      static double ProbTransform(double xx,double yy[2],double &weight,bool IsCenter);
      static bool RayScatterAngle(double wavelength, double &theta, double &phi,double anglerange[2],double &weight);
      static bool MieScatterAngle(double wavelength, double &theta, double &phi,double anglerange[2],double &weight);
      static double FreeIntgLength(double lengthrange[2],double &weight);
      static double FreePathLength(double z0,double dir0[3],double lengthrange[2],double &weight);
      static int IsScattering(double z0);
};

class Laser {
   public:
      static int jdebug;
      static TRandom3* prandom;
      static double TelSimDist;
      static double TelSimAngl;
      static double scale;
      static double unittime;
      static double intensity;
      static double intensity_err;
      static double wavelength0;
      static double wavelength0_err;
      static double frequency;
      static double pulsetime;
      static double spotrange;
      static double  divergence;

      ///position of the laser generator
      double lasercoo[3];	//in cm
      ///pointing direction of the laser generator
      double laserdir[2];
      double count_gen;
      int Time_gen;
      double time_gen;
      int ievent_gen;
      double wavelength_gen;	//in nm
      double coor_gen[3];	//in cm
      double dir_gen[3];
      vector<double> vgwav;
      vector<double> vgcoo[3];
      vector<double> vgdir[3];

      ///reweight during propagation in the atmosphere
      vector<double> vowei;

      ///coordinate and direction at telescope plane
      int Telindex;
      double coor_out[3];
      double dir_out[3];
      vector<double> votim;
      vector<int> votel;
      vector<double> vocoo[3];
      vector<double> vodir[3];

      WFCTAEvent* pwfc; //!

   public:
      void Init(int seed=0);
      void Release();
      void Reset();
      Laser(int seed=0) {Init(seed);}
      ~Laser() { Release(); }
      void SetParameters(char* filename=0);
      static void cross(double dir1[3],double dir2[3],double *dir3);
      static bool CartesianFrame(double zero[3],double coor_in[3],double dir_in[3],double *xdir,double *ydir,double *zdir);
      static double mindist(double zero[3],double coor_in[3],double dir_in[3],double *coor_min,bool &decrease);
      static double mindist(double coor_in[3],double dir_in[3],int &whichtel,double *coor_min,bool &decrease);
      static void PositionDis(double &xx,double &yy);
      static void DirectionDis(double &theta,double &phi);
      double WaveLengthGen();
      bool InitialGen();
      long int EventGen(int &Time,double &time,bool SimPulse=false);
      int Propagate(double &distance,double &weight);
      bool DoWFCTASim();
};

#endif
